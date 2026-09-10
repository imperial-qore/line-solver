/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_IAT2COUNTS_H
#define LINE_API_TRACE_TRACE_IAT2COUNTS_H

/**
 * Counting process of a trace: the number of arrivals in the window of length
 * `scale` that starts at each arrival epoch.
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_iat2counts.m,
 * cross-checked against
 * `jar/src/main/java/jline/api/trace/Trace_var.java#trace_iat2counts`.
 *
 * With CS the cumulative arrival epochs, MATLAB advances `cur` while
 * CS(cur+1) - CS(i) <= scale, i.e. it measures the window from the epoch of
 * the i-th arrival, and returns cur-i arrivals. The series is truncated at
 * the first window that reaches the end of the trace, because from there on
 * the count is censored.
 *
 * DIVERGENCE, MATLAB vs JAR: the JAR compares CS[cur+1]-CS[i] with its CS
 * indexed from 0, where CS[i] is the epoch of arrival i-1, not i. Its window
 * therefore starts one arrival too early and its counts are shifted by one
 * index relative to MATLAB. MATLAB is the reference here.
 *
 * ARITHMETIC: only additions and comparisons, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/**
 * @param S     inter-arrival times
 * @param scale window length
 * @return counts, one per arrival, truncated at the first censored window
 */
template <class T>
std::vector<long> trace_iat2counts(const std::vector<T>& S, const T& scale) {
    detail::require_nonempty(S, "trace_iat2counts");
    if (scale <= num_traits<T>::from_int(0))
        throw InputError("trace_iat2counts: the time scale must be positive");
    const long n = static_cast<long>(S.size());
    const std::vector<T> cs = detail::cumsum0(S);  // cs[k] = CS(k) of MATLAB

    std::vector<long> C;
    C.reserve(static_cast<std::size_t>(n > 0 ? n - 1 : 0));
    for (long i = 1; i <= n - 1; ++i) {
        long cur = i;
        while (cs[static_cast<std::size_t>(cur + 1)] - cs[static_cast<std::size_t>(i)] <= scale) {
            ++cur;
            if (cur == n) {  // the window has reached the end of the trace
                C.push_back(cur - i);
                return C;
            }
        }
        C.push_back(cur - i);
    }
    return C;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_IAT2COUNTS_H
