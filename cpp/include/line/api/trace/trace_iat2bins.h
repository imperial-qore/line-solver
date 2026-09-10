/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_IAT2BINS_H
#define LINE_API_TRACE_TRACE_IAT2BINS_H

/**
 * Bins a trace on a fixed time grid: the number of arrivals falling in each
 * interval ((i-1)*scale, i*scale], and the bin index of each arrival.
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_iat2bins.m,
 * cross-checked against
 * `jar/src/main/java/jline/api/trace/Trace_var.java#trace_iat2bins` (same
 * algorithm; the JAR stops one bin earlier because its outer loop runs to
 * `bins` while MATLAB runs to `bins+1`, so the JAR can drop the arrivals of
 * the final, partially filled bin).
 *
 * Unlike trace_iat2counts, the windows here are non-overlapping and anchored
 * at the origin, so sum(C) is the number of non-censored arrivals.
 *
 * ARITHMETIC: additions and comparisons plus one ceiling division for the bin
 * count, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** Return value of trace_iat2bins, mirroring [C, bC]. */
struct TraceBinsResult {
    std::vector<long> counts;      ///< arrivals per bin
    std::vector<long> membership;  ///< 1-based bin index of each binned arrival
};

/**
 * @param S     inter-arrival times
 * @param scale bin width
 */
template <class T>
TraceBinsResult trace_iat2bins(const std::vector<T>& S, const T& scale) {
    detail::require_nonempty(S, "trace_iat2bins");
    if (scale <= num_traits<T>::from_int(0))
        throw InputError("trace_iat2bins: the bin width must be positive");
    const long n = static_cast<long>(S.size());
    const std::vector<T> cs = detail::cumsum0(S);  // cs[k] = MATLAB CS(k)

    // bins = ceil((CS(end) - CS(1)) / scale), computed without a real ceil so
    // that the exact backends stay exact.
    const T span = cs[static_cast<std::size_t>(n)] - cs[1];
    long bins = 0;
    while (num_traits<T>::from_int(bins) * scale < span) ++bins;

    TraceBinsResult out;
    long cur = 1, last = 0;
    for (long i = 1; i <= bins + 1; ++i) {
        if (cur == n) break;
        while (cs[static_cast<std::size_t>(cur + 1)] <= num_traits<T>::from_int(i) * scale) {
            ++cur;
            if (cur == n) break;
        }
        if (static_cast<long>(out.counts.size()) < i) out.counts.resize(static_cast<std::size_t>(i), 0);
        out.counts[static_cast<std::size_t>(i - 1)] = cur - last;
        for (long k = 0; k < cur - last; ++k) out.membership.push_back(i);
        last = cur;
    }
    if (static_cast<long>(out.counts.size()) < bins) out.counts.resize(static_cast<std::size_t>(bins), 0);
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_IAT2BINS_H
