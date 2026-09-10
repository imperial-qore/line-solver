/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_JOINT_H
#define LINE_API_TRACE_TRACE_JOINT_H

/**
 * Joint moments of a trace, E[X_i^{k_1} X_{i+l_2}^{k_2} ...].
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_joint.m, cross-checked
 * against `jar/src/main/java/jline/api/trace/Trace_var.java#trace_joint`.
 *
 * The lag argument is a vector of INCREMENTS: MATLAB forms
 * `lag = sort(cumsum(lag))` and then shifts it so that its first entry is 0,
 * so trace_bicov's [1,i,j] means the triple (X_t, X_{t+i}, X_{t+i+j}).
 *
 * REFERENCE DEFECT (JAR): Trace_var.trace_joint sorts the lag vector but
 * never takes the cumulative sum, and then indexes it with
 * `adjustedLag[min(j, adjustedLag.length-1)]`, which also silently reuses the
 * last lag when order is longer than lag. It therefore computes a different
 * joint moment from MATLAB for every lag vector that is not already
 * cumulative -- including the [1,i,j] grid that Trace_var.trace_bicov feeds
 * it, so the JAR bicovariance is wrong wherever i or j differs from 1.
 * MATLAB is the reference and is what this port implements; a lag and order
 * of different lengths is rejected instead of being padded.
 *
 * ARITHMETIC: products of integer powers of the samples and one division,
 * exact in Rational.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/**
 * @param S     the trace
 * @param lag   lag increments; cumulated and shifted to start at 0
 * @param order the exponent of each factor, same length as lag
 */
template <class T>
T trace_joint(const std::vector<T>& S, const std::vector<int>& lag,
              const std::vector<unsigned>& order) {
    detail::require_nonempty(S, "trace_joint");
    if (lag.size() != order.size())
        throw InputError("trace_joint: lag and order must have the same length");
    if (lag.empty()) throw InputError("trace_joint: an empty lag vector");

    std::vector<long> L(lag.size());
    long acc = 0;
    for (std::size_t i = 0; i < lag.size(); ++i) {
        acc += lag[i];
        L[i] = acc;
    }
    std::sort(L.begin(), L.end());
    const long base = L[0];
    for (std::size_t i = 0; i < L.size(); ++i) L[i] -= base;
    const long maxlag = L.back();
    if (maxlag < 0) throw InputError("trace_joint: the cumulated lags are not increasing");

    const long n = static_cast<long>(S.size());
    const long len = n - maxlag;
    if (len <= 0) throw InputError("trace_joint: the lag span exceeds the trace length");

    T sum = num_traits<T>::from_int(0);
    for (long t = 0; t < len; ++t) {
        T prod = num_traits<T>::from_int(1);
        for (std::size_t i = 0; i < L.size(); ++i)
            prod *= num_pow_int(S[static_cast<std::size_t>(t + L[i])], order[i]);
        sum += prod;
    }
    return sum / num_traits<T>::from_int(len);
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_JOINT_H
