/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SYM_SYM_ENGINES_H
#define LINE_API_SYM_SYM_ENGINES_H

/**
 * Resolves the symbolic backend to use, and owns the container that serves it.
 *
 * Port of jline.api.sym.SymEngines. Resolution order, the same in MATLAB
 * (SAGE.m), the JAR and Python (line_solver.api.sym):
 *   1. an explicit URL, from solver options or the `requested` argument;
 *   2. the LINE_SAGE_URL environment variable;
 *   3. a line-sage-rest service already listening on a conventional port;
 *   4. a container started here from a locally present image;
 *   5. nothing, in which case the caller keeps whatever native algebra it has,
 *      or reports that no backend is configured.
 *
 * STEP 3 VERIFIES IDENTITY THROUGH /api/v1/info rather than trusting the port:
 * every imperialqore line-*-rest service listens on 8080 by convention, so a
 * health probe alone would happily accept the LQNS service and then fail on the
 * first symbolic request with an unrecognizable error.
 *
 * The container started in step 4 is reused for the life of the process and
 * stopped by an atexit handler. It is bound to an ephemeral host port, so
 * several processes, or a process alongside a hand-started service, do not
 * collide. An atexit handler does NOT run on a signal or on _exit, so a
 * container may outlive a killed process; `docker ps` shows it under the name
 * line-sage-rest-<port> and it was started with --rm, so stopping it removes it.
 *
 * A PULL HAPPENS ONLY ON EXPLICIT OPT-IN, i.e. the "sage" keyword or a named
 * image. Bare "auto"/"true"/"" keep the native backend unless the image is
 * already local, so leaving the symbolic option on auto never triggers a
 * multi-gigabyte download, and the pull itself is refused when the Docker
 * storage location is short of space (see line/io/docker_image.h).
 *
 * DIVERGENCE FROM THE JAR: an https:// URL is refused by name rather than used,
 * because this port's HTTP client has no TLS (see line/util/http.h). Refusing
 * loudly is the point -- reporting "no backend" for a service that is up and
 * merely unreachable over plaintext would send the caller looking in the wrong
 * place.
 */

#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

#include "line/api/sym/sage_rest_engine.h"
#include "line/api/sym/sym_engine.h"
#include "line/io/docker_image.h"
#include "line/util/error.h"
#include "line/util/subprocess.h"

namespace line {
namespace sym {

/** Image serving the symbolic REST API. */
inline const char* const SYM_DOCKER_IMAGE = "imperialqore/line-sage-rest:latest";
/** Environment variable naming a service to use. */
inline const char* const SYM_URL_ENV = "LINE_SAGE_URL";
/** Seconds to wait for a container to report healthy. */
inline constexpr int SYM_STARTUP_TIMEOUT_SECONDS = 120;

/** Fallback tags, tried in order after SYM_DOCKER_IMAGE. */
inline std::vector<std::string> sym_docker_image_candidates() {
    return std::vector<std::string>{"imperialqore/line-sage-rest:latest",
                                    "imperialqore/line-sage-rest"};
}

/** Ports probed for an already running service, in order. */
inline std::vector<int> sym_probe_ports() { return std::vector<int>{8085, 8080}; }

namespace detail {

/** Process-wide record of the container this process started, if any. */
struct SymState {
    std::mutex mutex;
    std::shared_ptr<SageRestEngine> started;
    std::string container;
    bool atexitRegistered = false;
};

inline SymState& sym_state() {
    static SymState state;
    return state;
}

inline std::string lower(const std::string& s) {
    std::string out(s);
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = static_cast<char>(std::tolower(out[i]));
    return out;
}

/** An unused local TCP port, obtained the way the JAR does: bind port 0. */
inline int free_port() {
    const int fd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) throw SymEngineError("SymEngines: cannot open a socket to pick a free port");
    struct sockaddr_in addr;
    std::memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    addr.sin_port = 0;
    if (::bind(fd, reinterpret_cast<struct sockaddr*>(&addr), sizeof(addr)) != 0) {
        ::close(fd);
        throw SymEngineError("SymEngines: cannot bind a free port");
    }
    socklen_t len = sizeof(addr);
    if (::getsockname(fd, reinterpret_cast<struct sockaddr*>(&addr), &len) != 0) {
        ::close(fd);
        throw SymEngineError("SymEngines: cannot read the bound port");
    }
    const int port = static_cast<int>(ntohs(addr.sin_port));
    ::close(fd);
    return port;
}

/** Checks that a service is line-sage-rest and not another line-*-rest one. */
inline bool is_sage_service(const SageRestEngine& engine) {
    try {
        return engine.info().contains("sage_version");
    } catch (const Error&) {
        return false;
    }
}

inline void sleep_millis(long millis) {
    struct timespec ts;
    ts.tv_sec = millis / 1000;
    ts.tv_nsec = (millis % 1000) * 1000000L;
    ::nanosleep(&ts, nullptr);
}

}  // namespace detail

/** Stops the container started by this process, if any. */
inline void sym_stop_container() {
    detail::SymState& st = detail::sym_state();
    std::lock_guard<std::mutex> guard(st.mutex);
    if (st.container.empty()) return;
    util::capture({"docker", "stop", "-t", "1", st.container}, 30);
    st.container.clear();
    st.started.reset();
}

/**
 * @return the first locally present image tag, or the empty string if none is
 */
inline std::string sym_find_image() {
    const std::vector<std::string> candidates = sym_docker_image_candidates();
    for (std::size_t i = 0; i < candidates.size(); ++i)
        if (io::docker_has_local_image(candidates[i])) return candidates[i];
    return std::string();
}

namespace detail {

/**
 * Pulls `target` if the Docker storage location has room; returns the tag on
 * success, else the empty string. On refusal the caller keeps its native
 * algebra rather than failing.
 */
inline std::string pull_image(const std::string& target) {
    if (!io::docker_has_storage_for(target)) {
        std::cerr << "[LINE] Skipping docker pull of " << target
                  << ": insufficient free space at the Docker storage location; "
                  << "keeping the native symbolic backend." << std::endl;
        return std::string();
    }
    std::cout << "[LINE] Pulling Docker image " << target << " (this may take a while)..."
              << std::endl;
    if (io::docker_pull(target) && io::docker_has_local_image(target)) return target;
    return std::string();
}

/** Starts the service in a container and waits for it to report healthy. */
inline std::shared_ptr<SageRestEngine> start_container(const std::string& image) {
    const int port = free_port();
    const std::string name = "line-sage-rest-" + std::to_string(port);
    const util::ProcResult run =
        util::capture({"docker", "run", "-d", "--rm", "--name", name, "-p",
                       std::to_string(port) + ":8080", image},
                      120);
    if (run.exitCode != 0 || util::trim(run.out).empty())
        throw SymEngineError("SymEngines: could not start " + image);

    SymState& st = sym_state();
    {
        std::lock_guard<std::mutex> guard(st.mutex);
        st.container = name;
        if (!st.atexitRegistered) {
            std::atexit(&sym_stop_container);
            st.atexitRegistered = true;
        }
    }

    std::shared_ptr<SageRestEngine> engine =
        std::make_shared<SageRestEngine>("http://localhost:" + std::to_string(port));
    for (int waited = 0; waited < SYM_STARTUP_TIMEOUT_SECONDS * 1000; waited += 500) {
        if (engine->isAvailable()) {
            if (engine->isUsable()) {
                std::lock_guard<std::mutex> guard(st.mutex);
                st.started = engine;
                return engine;
            }
            // Booted, but its arithmetic dies on this CPU. Keeping it running
            // would only cost memory, and returning it would hand the caller a
            // backend that kills every request.
            sym_stop_container();
            throw SymEngineError("SymEngines: container " + name +
                                 " answers but cannot evaluate on this CPU");
        }
        sleep_millis(500);
    }
    sym_stop_container();
    throw SymEngineError("SymEngines: container " + name + " did not become healthy within " +
                         std::to_string(SYM_STARTUP_TIMEOUT_SECONDS) + " s");
}

}  // namespace detail

/**
 * Resolves an engine.
 *
 * @param requested "" or "auto" to search, a URL to use a specific service,
 *                  "none" to disable the backend, or an image name to start
 * @return an engine, or a null pointer if no backend could be resolved
 */
inline std::shared_ptr<SymEngine> sym_resolve(const std::string& requested = "auto") {
    const std::string req = util::trim(requested);
    const std::string reqLower = detail::lower(req);
    if (reqLower == "none" || reqLower == "off") return std::shared_ptr<SymEngine>();

    if (req.compare(0, 8, "https://") == 0)
        throw UnsupportedError(
            "SymEngines: this port's HTTP client has no TLS, so the symbolic service must be "
            "reached over http://; terminate TLS in front of it or use a local container");
    if (req.compare(0, 7, "http://") == 0) {
        std::shared_ptr<SageRestEngine> engine = std::make_shared<SageRestEngine>(req);
        return engine->isAvailable() && engine->isUsable() ? engine : std::shared_ptr<SymEngine>();
    }

    const char* env = std::getenv(SYM_URL_ENV);
    if (env != nullptr && !util::trim(env).empty()) {
        const std::string url = util::trim(env);
        if (url.compare(0, 8, "https://") == 0) {
            std::cerr << "[LINE] Ignoring " << SYM_URL_ENV
                      << ": this port's HTTP client has no TLS." << std::endl;
        } else {
            std::shared_ptr<SageRestEngine> engine = std::make_shared<SageRestEngine>(url);
            if (engine->isAvailable() && engine->isUsable()) return engine;
        }
    }

    {
        detail::SymState& st = detail::sym_state();
        std::shared_ptr<SageRestEngine> cached;
        {
            std::lock_guard<std::mutex> guard(st.mutex);
            cached = st.started;
        }
        if (cached && cached->isAvailable() && cached->isUsable()) return cached;
    }

    const std::vector<int> ports = sym_probe_ports();
    for (std::size_t i = 0; i < ports.size(); ++i) {
        std::shared_ptr<SageRestEngine> engine =
            std::make_shared<SageRestEngine>("http://localhost:" + std::to_string(ports[i]));
        if (detail::is_sage_service(*engine) && engine->isUsable()) return engine;
    }

    const bool search =
        req.empty() || reqLower == "auto" || reqLower == "true" || reqLower == "sage";
    std::string image = search ? sym_find_image() : req;
    if (image.empty() && reqLower == "sage")
        image = detail::pull_image(SYM_DOCKER_IMAGE);
    else if (!search && !image.empty() && !io::docker_has_local_image(image))
        image = detail::pull_image(image);
    if (image.empty()) return std::shared_ptr<SymEngine>();

    try {
        return detail::start_container(image);
    } catch (const Error&) {
        return std::shared_ptr<SymEngine>();
    }
}

}  // namespace sym
}  // namespace line

#endif  // LINE_API_SYM_SYM_ENGINES_H
