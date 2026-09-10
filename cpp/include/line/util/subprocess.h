/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_SUBPROCESS_H
#define LINE_UTIL_SUBPROCESS_H

/**
 * Running an external command and capturing its output, with a deadline.
 *
 * The port needs this for the Docker orchestration behind the symbolic backend,
 * where the JAR uses ProcessBuilder. popen() is not enough: it offers no timeout
 * and runs the command through a shell, so an image tag would be word-split and
 * glob-expanded. This forks and execs an argv vector directly, so no shell ever
 * sees the arguments, and enforces the deadline by polling the pipe and killing
 * the child when it expires.
 *
 * A TIMED-OUT CHILD IS KILLED AND REAPED, never merely abandoned: an orphaned
 * `docker run` would keep a container alive that nobody holds the name of.
 */

#include <cerrno>
#include <csignal>
#include <cstddef>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

#include <fcntl.h>
#include <poll.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace line {
namespace util {

/** Outcome of a captured command. */
struct ProcResult {
    int exitCode = -1;     ///< Exit status, or -1 when the command could not run
    std::string out;       ///< Everything the command wrote to stdout
    bool timedOut = false; ///< True when the deadline expired and the child was killed
};

namespace detail {

/** Milliseconds on the monotonic clock, for deadlines that survive a clock step. */
inline long monotonic_millis() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return static_cast<long>(ts.tv_sec) * 1000L + ts.tv_nsec / 1000000L;
}

inline void exec_never_returns(const std::vector<std::string>& argv) {
    std::vector<char*> raw;
    raw.reserve(argv.size() + 1);
    for (std::size_t i = 0; i < argv.size(); ++i)
        raw.push_back(const_cast<char*>(argv[i].c_str()));
    raw.push_back(nullptr);
    ::execvp(raw[0], &raw[0]);
    ::_exit(127);
}

}  // namespace detail

/**
 * Runs a command, capturing stdout and discarding stderr.
 *
 * `mergeStderr` sends the child's stderr down the SAME pipe instead, which is
 * what MATLAB's `[status, cmdout] = system(cmd)` does and what a wrapper needs
 * on the failure path: an external engine reports why it refused on stderr, and
 * dropping it leaves the wrapper reporting an exit code and nothing else.
 *
 * @param argv           the command and its arguments, argv[0] resolved on PATH
 * @param timeoutSeconds deadline; not positive waits indefinitely
 * @param mergeStderr    capture stderr too, interleaved with stdout
 * @return the exit code, the captured output and whether the deadline expired
 */
inline ProcResult capture(const std::vector<std::string>& argv, int timeoutSeconds,
                          bool mergeStderr = false) {
    ProcResult r;
    if (argv.empty()) return r;

    int pipefd[2];
    if (::pipe(pipefd) != 0) return r;

    const pid_t pid = ::fork();
    if (pid < 0) {
        ::close(pipefd[0]);
        ::close(pipefd[1]);
        return r;
    }
    if (pid == 0) {
        ::close(pipefd[0]);
        ::dup2(pipefd[1], STDOUT_FILENO);
        if (mergeStderr) {
            ::dup2(pipefd[1], STDERR_FILENO);
        } else {
            const int devnull = ::open("/dev/null", O_WRONLY);
            if (devnull >= 0) {
                ::dup2(devnull, STDERR_FILENO);
                ::close(devnull);
            }
        }
        ::close(pipefd[1]);
        detail::exec_never_returns(argv);
    }

    ::close(pipefd[1]);
    const long deadline =
        timeoutSeconds > 0 ? detail::monotonic_millis() + timeoutSeconds * 1000L : -1;
    char buf[4096];
    while (true) {
        int wait = -1;
        if (deadline >= 0) {
            const long left = deadline - detail::monotonic_millis();
            if (left <= 0) {
                r.timedOut = true;
                break;
            }
            wait = static_cast<int>(left);
        }
        struct pollfd pfd;
        pfd.fd = pipefd[0];
        pfd.events = POLLIN;
        pfd.revents = 0;
        const int pr = ::poll(&pfd, 1, wait);
        if (pr == 0) {
            r.timedOut = true;
            break;
        }
        if (pr < 0) {
            if (errno == EINTR) continue;
            break;
        }
        const ssize_t k = ::read(pipefd[0], buf, sizeof(buf));
        if (k > 0) {
            r.out.append(buf, static_cast<std::size_t>(k));
            continue;
        }
        if (k < 0 && errno == EINTR) continue;
        break;  // end of stream
    }
    ::close(pipefd[0]);

    if (r.timedOut) {
        ::kill(pid, SIGKILL);
        int status = 0;
        ::waitpid(pid, &status, 0);
        r.exitCode = -1;
        return r;
    }
    int status = 0;
    if (::waitpid(pid, &status, 0) == pid && WIFEXITED(status))
        r.exitCode = WEXITSTATUS(status);
    return r;
}

/**
 * Runs a command with the parent's stdout and stderr, e.g. so a `docker pull`
 * streams its progress. There is no deadline: a pull is long by nature and the
 * user can see it working.
 *
 * @param argv the command and its arguments
 * @return the exit code, or -1 if the command could not run
 */
inline int run_inherit(const std::vector<std::string>& argv) {
    if (argv.empty()) return -1;
    const pid_t pid = ::fork();
    if (pid < 0) return -1;
    if (pid == 0) detail::exec_never_returns(argv);
    int status = 0;
    if (::waitpid(pid, &status, 0) == pid && WIFEXITED(status)) return WEXITSTATUS(status);
    return -1;
}

/** Trims ASCII whitespace from both ends, as Java's String.trim() does. */
inline std::string trim(const std::string& s) {
    const std::size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return std::string();
    const std::size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

}  // namespace util
}  // namespace line

#endif  // LINE_UTIL_SUBPROCESS_H
