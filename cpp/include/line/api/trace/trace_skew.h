/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_SKEW_H
#define LINE_API_TRACE_TRACE_SKEW_H

/**
 * Bias-corrected sample skewness (MATLAB's skewness(S,0), equivalently the
 * G1 estimator).
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_skew.m, cross-checked
 * against jar/src/main/java/jline/api/trace/TraceSkew.java, which delegates
 * to Apache Commons Math's Skewness. The two are algebraically IDENTICAL:
 * Apache computes n/((n-1)(n-2)) * sum d^3 / s^3 with s the n-1 standard
 * deviation, MATLAB computes m3/m2^{3/2} * sqrt((n-1)/n) * n/(n-2), and both
 * reduce to
 *
 *   G1 = sqrt(n(n-1)) / (n-2) * m3 / m2^{3/2},   m_k = (1/n) sum (x-xbar)^k.
 *
 * That closed form is what is evaluated here, so the result is symmetric in
 * the two references rather than favouring one rounding order.
 *
 * ARITHMETIC: a 3/2 power of the second central moment.
 *   static_assert(num_traits<T>::has_transcendental)
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_mean.h"
#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** @param S the trace, at least 3 samples. */
template <class T>
T trace_skew(const std::vector<T>& S) {
    static_assert(num_traits<T>::has_transcendental,
                  "trace_skew requires transcendental arithmetic");
    detail::require_nonempty(S, "trace_skew");
    const long n = static_cast<long>(S.size());
    if (n < 3) throw InputError("trace_skew: the skewness needs at least three samples");
    const T mu = trace_mean(S);
    T m2 = num_traits<T>::from_int(0), m3 = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < S.size(); ++i) {
        const T d = S[i] - mu;
        m2 += d * d;
        m3 += d * d * d;
    }
    const T nt = num_traits<T>::from_int(n);
    m2 /= nt;
    m3 /= nt;
    if (m2 == num_traits<T>::from_int(0))
        throw NumericError("trace_skew: the trace is constant, the skewness is undefined");
    const T corr = detail::num_sqrt(T(nt * num_traits<T>::from_int(n - 1))) / num_traits<T>::from_int(n - 2);
    return corr * m3 / (m2 * detail::num_sqrt(m2));
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_SKEW_H
