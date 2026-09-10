/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_BICOV_H
#define LINE_API_TRACE_TRACE_BICOV_H

/**
 * Bicovariance of a trace on a lag grid.
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_bicov.m, cross-checked
 * against `jar/src/main/java/jline/api/trace/Trace_var.java#trace_bicov` (same
 * grid construction; the JAR values differ because its trace_joint drops the
 * cumulative sum, see trace_joint.h).
 *
 * For every ordered pair (i,j) drawn from the grid the third-order joint
 * moment E[X_t X_{t+i} X_{t+i+j}] is evaluated. Note that both references
 * return the raw joint moment, not a centred cumulant, despite the name.
 *
 * ARITHMETIC: a vector of joint moments, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_joint.h"
#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** Return value of trace_bicov, mirroring [BiCov, BiCovLags]. */
template <class T>
struct TraceBicovResult {
    std::vector<T> bicov;                ///< one moment per grid pair
    std::vector<std::vector<int>> lags;  ///< the [1,i,j] increment triples
};

/**
 * @param S    the trace
 * @param grid the lag values swept in both positions
 */
template <class T>
TraceBicovResult<T> trace_bicov(const std::vector<T>& S, const std::vector<int>& grid) {
    detail::require_nonempty(S, "trace_bicov");
    TraceBicovResult<T> out;
    const std::vector<unsigned> order(3, 1u);
    for (std::size_t a = 0; a < grid.size(); ++a) {
        for (std::size_t b = 0; b < grid.size(); ++b) {
            std::vector<int> lag(3);
            lag[0] = 1;
            lag[1] = grid[a];
            lag[2] = grid[b];
            out.lags.push_back(lag);
            out.bicov.push_back(trace_joint(S, lag, order));
        }
    }
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_BICOV_H
