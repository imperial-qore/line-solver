/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_TEMPDIR_H
#define LINE_UTIL_TEMPDIR_H

/**
 * A scratch directory for the subprocess wrappers, the port's `lineTempName`.
 *
 * The wrappers hand an external binary a FILE, not a stream, so they need a
 * private directory to put it in and to collect the tool's output from. The
 * name comes from mkdtemp rather than from a counter or a PID: two solves of the
 * same model running side by side (a parity sweep does exactly that) would
 * otherwise write the same `model.jmva` and read each other's results.
 *
 * The directory is REMOVED IN THE DESTRUCTOR, including on the exception path,
 * because the wrappers throw on a tool failure and a leaked scratch directory
 * per failed solve accumulates silently. `keep()` suppresses the removal for
 * the case where the caller wants to inspect what was sent.
 */

#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <string>
#include <vector>

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "line/util/error.h"

namespace line {
namespace util {

namespace detail {

/** The string without its surrounding whitespace. */
inline std::string trim_ws(const std::string& s) {
    std::size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

/**
 * Create a directory and every missing parent, as MATLAB's `mkdir` does.
 *
 * "It is already there" is success: two solves staging side by side race on the
 * shared workspace directory, and the loser of that race has what it asked for.
 */
inline void mkdir_p(const std::string& dir) {
    for (std::size_t at = dir.find('/', 1); ; at = dir.find('/', at + 1)) {
        const std::string part = (at == std::string::npos) ? dir : dir.substr(0, at);
        if (::mkdir(part.c_str(), 0700) != 0 && errno != EEXIST)
            throw InputError("cannot create the staging directory '" + part + "'");
        if (at == std::string::npos) break;
    }
}

}  // namespace detail

/**
 * Create a private scratch directory named after its caller.
 *
 * LINE_WORKSPACE_ROOT RELOCATES EVERY STAGED MODEL, exactly as it does in
 * `lineTempName.m`, `SysUtils.lineTempName` and the native-Python wrappers.
 * run-tests.sh sets it when it wraps lqns, lqsim or qnsolver in a container: the
 * shim bind-mounts that root and nothing else, so a model left under TMPDIR is
 * invisible to the binary, which then writes no result and is reported as having
 * rejected the model. No wrapper here needs to know a container is involved.
 *
 * `mountable` is the second argument of `lineTempName`: a snap-confined Docker
 * daemon bind-mounts neither the system temp dir nor a dot-directory, so a
 * caller that hands the directory to a container asks for a HOME-rooted one.
 *
 * @param prefix a short tag naming the caller, as `lineTempName('qns')` does
 * @param mountable root under $HOME/.line rather than TMPDIR
 */
inline std::string make_temp_dir(const std::string& prefix, bool mountable = false) {
    std::string root;
    const char* ws = std::getenv("LINE_WORKSPACE_ROOT");
    if (ws != nullptr) root = detail::trim_ws(std::string(ws));
    if (root.empty() && mountable) {
        const char* home = std::getenv("HOME");
        if (home != nullptr && *home) root = std::string(home) + "/.line";
    }
    if (!root.empty()) {
        root += "/line_workspace/" + prefix;
        detail::mkdir_p(root);
    } else {
        const char* base = std::getenv("TMPDIR");
        root = (base && *base ? std::string(base) : std::string("/tmp"));
    }
    const std::string tmpl = root + "/line_" + prefix + "_XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    if (::mkdtemp(&buf[0]) == nullptr)
        throw InputError("cannot create a temporary directory under '" + tmpl + "'");
    return std::string(&buf[0]);
}

class TempDir {
public:
    /**
     * @param prefix a short tag naming the caller, as `lineTempName('qns')` does
     * @param mountable root it where a confined Docker daemon can bind-mount it
     */
    explicit TempDir(const std::string& prefix, bool mountable = false)
        : path_(make_temp_dir(prefix, mountable)) {}

    ~TempDir() {
        if (!keep_ && !path_.empty()) remove_tree(path_);
    }

    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;

    /** The directory itself, with no trailing separator. */
    const std::string& path() const { return path_; }

    /** A file inside it. */
    std::string file(const std::string& name) const { return path_ + "/" + name; }

    /** Leave the directory in place, for a caller that wants to inspect it. */
    void keep() { keep_ = true; }

private:
    static void remove_tree(const std::string& dir) {
        DIR* d = ::opendir(dir.c_str());
        if (d != nullptr) {
            struct dirent* e;
            while ((e = ::readdir(d)) != nullptr) {
                const std::string n(e->d_name);
                if (n == "." || n == "..") continue;
                const std::string child = dir + "/" + n;
                struct stat st;
                if (::lstat(child.c_str(), &st) == 0 && S_ISDIR(st.st_mode))
                    remove_tree(child);
                else
                    ::unlink(child.c_str());
            }
            ::closedir(d);
        }
        ::rmdir(dir.c_str());
    }

    std::string path_;
    bool keep_ = false;
};

}  // namespace util
}  // namespace line

#endif  // LINE_UTIL_TEMPDIR_H
