/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_DOCKER_IMAGE_H
#define LINE_IO_DOCKER_IMAGE_H

/**
 * Docker primitives for the backends that legitimately ship an image.
 *
 * Port of jar/src/main/java/jline/io/DockerImage.java. In the JAR these are
 * shared by the JMT backend (imperialqore/jmt-rest) and the Sage symbolic
 * engine (imperialqore/line-sage-rest); this port has only the second so far.
 * Each backend owns its image name; this file only answers whether the daemon
 * is up, whether an image is present, whether there is room to pull it, and
 * performs the pull.
 *
 * LQNS, lqsim and qnsolver are deliberately absent, here as in the JAR: their
 * licence is an evaluation agreement that forbids redistribution, so LINE runs
 * them only from a binary the user installed.
 *
 * THE STORAGE GUARD IS NOT DECORATION. A pull that fills the filesystem backing
 * the Docker root breaks every container on the machine, not just LINE's, so a
 * pull is refused unless the free space clears the larger of 2 GiB and three
 * times the compressed manifest size. When the free space cannot be determined
 * the pull is ALLOWED: an unknown is not evidence of a full disk, and refusing
 * would strand users whose daemon does not report a root dir.
 */

#include <cstddef>
#include <cstdlib>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/statvfs.h>

#include "line/util/subprocess.h"

namespace line {
namespace io {

/** Conservative free-space floor required before a pull (2 GiB). */
static const long long DOCKER_DEFAULT_MIN_FREE_BYTES = 2LL * 1024 * 1024 * 1024;

/** True if the Docker daemon is reachable. */
inline bool docker_daemon_available() {
    const util::ProcResult r = util::capture({"docker", "info"}, 15);
    return r.exitCode == 0;
}

/**
 * @param image the image tag
 * @return true if it is already present in the local Docker store
 */
inline bool docker_has_local_image(const std::string& image) {
    if (image.empty()) return false;
    const util::ProcResult r = util::capture({"docker", "images", "-q", image}, 15);
    return r.exitCode == 0 && !util::trim(r.out).empty();
}

/**
 * Pulls an image, streaming Docker's progress to stdout and stderr. No timeout.
 *
 * @param image the image tag
 * @return true on success
 */
inline bool docker_pull(const std::string& image) {
    if (image.empty()) return false;
    return util::run_inherit({"docker", "pull", image}) == 0;
}

namespace detail {

/** Usable bytes on the filesystem backing the Docker root dir; -1 if unknown. */
inline long long docker_free_bytes() {
    const util::ProcResult r =
        util::capture({"docker", "info", "--format", "{{.DockerRootDir}}"}, 15);
    std::string root = r.exitCode == 0 ? util::trim(r.out) : std::string();
    if (root.empty()) root = "/var/lib/docker";

    // The root dir may not exist for, or be readable by, this user: walk up to
    // an existing ancestor so the statvfs reports a real filesystem.
    struct stat st;
    while (!root.empty() && ::stat(root.c_str(), &st) != 0) {
        const std::size_t slash = root.find_last_of('/');
        if (slash == std::string::npos) break;
        root = slash == 0 ? std::string("/") : root.substr(0, slash);
        if (root == "/") break;
    }
    if (root.empty()) root = "/";

    struct statvfs vfs;
    if (::statvfs(root.c_str(), &vfs) != 0) return -1;
    const long long usable =
        static_cast<long long>(vfs.f_bavail) * static_cast<long long>(vfs.f_frsize);
    return usable > 0 ? usable : -1;
}

/** On-disk estimate (compressed layer sizes x3), or 0 if it cannot be determined. */
inline long long docker_estimate_image_bytes(const std::string& image) {
    const util::ProcResult r = util::capture({"docker", "manifest", "inspect", image}, 30);
    if (r.exitCode != 0 || r.out.empty()) return 0;
    long long sum = 0;
    const std::string& json = r.out;
    const std::string key = "\"size\"";
    std::size_t at = 0;
    while ((at = json.find(key, at)) != std::string::npos) {
        std::size_t i = at + key.size();
        while (i < json.size() && (json[i] == ' ' || json[i] == ':' || json[i] == '\t')) ++i;
        if (i < json.size() && json[i] >= '0' && json[i] <= '9')
            sum += std::atoll(json.c_str() + i);
        at = i;
    }
    return sum > 0 ? sum * 3 : 0;
}

/** LINE_DOCKER_MIN_FREE_BYTES, or -1 when unset or unparsable. */
inline long long docker_override_min_free_bytes() {
    const char* v = std::getenv("LINE_DOCKER_MIN_FREE_BYTES");
    if (v == nullptr || *v == '\0') return -1;
    const long long parsed = std::atoll(v);
    return parsed > 0 ? parsed : -1;
}

}  // namespace detail

/**
 * @param image the image tag
 * @return true if the Docker storage location has room for it
 */
inline bool docker_has_storage_for(const std::string& image) {
    const long long free = detail::docker_free_bytes();
    if (free < 0) return true;  // could not determine; do not block the pull
    const long long override = detail::docker_override_min_free_bytes();
    if (override > 0) return free >= override;
    const long long est = detail::docker_estimate_image_bytes(image);
    const long long required = est > DOCKER_DEFAULT_MIN_FREE_BYTES ? est
                                                                   : DOCKER_DEFAULT_MIN_FREE_BYTES;
    return free >= required;
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_DOCKER_IMAGE_H
