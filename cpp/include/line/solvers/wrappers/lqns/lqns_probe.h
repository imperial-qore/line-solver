/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_WRAPPERS_LQNS_LQNS_PROBE_H
#define LINE_SOLVERS_WRAPPERS_LQNS_LQNS_PROBE_H

/**
 * Is a usable lqns installed on this machine?
 *
 * Its own header, small enough for SolverAUTO to include: the layered arm of
 * `chooseAvgSolverHeur` gates its LQNS candidate on `SolverLQNS.isAvailable()`,
 * so the chooser needs the answer without taking on the wrapper, the .lqnx
 * writer and the XML DOM behind it.
 *
 * THE ANSWER IS CACHED for the life of the process. A binary does not appear or
 * vanish mid-run, and the probe is a fork and an exec: SolverAUTO would
 * otherwise pay for one on every layered choice it makes.
 */

#include <cstdlib>
#include <string>
#include <vector>

#include "line/util/subprocess.h"

namespace line {
namespace lqns {

namespace detail {

/** The first line of `lqns -V -H`, empty when there is no binary to ask. */
inline std::string probe_version_banner() {
    std::vector<std::string> argv;
    argv.push_back("lqns");
    argv.push_back("-V");
    argv.push_back("-H");
    const util::ProcResult r = util::capture(argv, 20);
    if (r.exitCode < 0 || r.timedOut) return std::string();
    const std::size_t nl = r.out.find('\n');
    return util::trim(nl == std::string::npos ? r.out : r.out.substr(0, nl));
}

}  // namespace detail

/** The version banner of the local lqns, empty when there is none. */
inline const std::string& lqns_version() {
    static const std::string banner = detail::probe_version_banner();
    return banner;
}

/** The major release number, 0 when no binary answered. */
inline int lqns_major_version() {
    const std::string& b = lqns_version();
    const std::size_t at = b.find("Version ");
    if (at == std::string::npos) return 0;
    return std::atoi(b.c_str() + at + 8);
}

/**
 * True when lqns is installed AND is a release this port speaks.
 *
 * LINE requires 6.0 or greater: the 5.x .lqxo grammar spells several result
 * attributes differently, so an older binary would parse into a table of NaNs
 * rather than fail, and the caller would read "not reported" as an answer.
 */
inline bool lqns_is_available() { return lqns_major_version() >= 6; }

}  // namespace lqns
}  // namespace line

#endif  // LINE_SOLVERS_WRAPPERS_LQNS_LQNS_PROBE_H
