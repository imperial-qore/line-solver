/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_REG_API_DISPATCH_H
#define LINE_REG_API_DISPATCH_H

/**
 * Direct invocation of a single API function from named JSON arguments.
 *
 * This is the entry point behind the CLI's --api flag and the one a pybind11
 * or MEX gateway should call, so that all hosts share one dispatch table and
 * one conversion policy (include/line/reg/api_json.h).
 *
 * Three refusals, all explicit and all by name, because the port is
 * incremental and a caller must never be able to mistake "not wired" for a
 * result:
 *   - a name absent from the registry is not ported to C++ at all;
 *   - a name in the registry but not in the dispatch table is ported and
 *     tested but not yet exposed over this boundary;
 *   - an arithmetic the registry does not list for that function is refused
 *     naming the modes it does list.
 * No path returns an empty result, a zero, or a value from a different
 * arithmetic than the caller asked for.
 */

#include <string>
#include <vector>

#include "line/reg/api_json.h"
#include "line/reg/registry.h"

namespace line {
namespace reg {

/** The arithmetic a call runs at: Real also carries the precision tier. */
struct ArithSpec {
    Arith mode = Arith::Double;
    unsigned digits = 0;  ///< significant decimal digits, Real only

    /** Canonical text, e.g. "double", "exact", "real:50". */
    std::string str() const;
};

/**
 * Parse --arith. Accepts "double", "exact" and "real:<digits>"; anything else
 * is an InputError listing the accepted forms. The port instantiates the
 * high-precision backend at three fixed tiers (50, 100, 200 digits), so a
 * request between tiers is rounded UP to the next one -- never down, so a
 * caller never silently receives less precision than asked -- and str()
 * reports the tier actually used. Above 200 digits the call is refused.
 */
ArithSpec parse_arith(const std::string& text);

/** Names exposed over this boundary, sorted; a subset of the registry. */
std::vector<std::string> api_exposed_functions();

/** True when the named function has a dispatch entry. */
bool api_is_exposed(const std::string& name);

/**
 * Invoke one API function.
 *
 * @param name  MATLAB function name, e.g. "pfqn_ca"
 * @param arith text of --arith, e.g. "exact" or "real:50"
 * @param args  object keyed by the MATLAB parameter names
 * @return {"function", "arith", "results"}, results keyed by the MATLAB output
 *         names and encoded per the policy in api_json.h
 * @throws UnsupportedError, InputError, NumericError
 */
Json api_invoke(const std::string& name, const std::string& arith, const Json& args);

/** Human-readable rendering of the object api_invoke returns, for -o readable. */
std::string api_render_readable(const Json& result);

}  // namespace reg
}  // namespace line

#endif  // LINE_REG_API_DISPATCH_H
