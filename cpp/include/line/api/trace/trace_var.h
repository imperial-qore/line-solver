/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_VAR_H
#define LINE_API_TRACE_TRACE_VAR_H

/**
 * Sample variance of a trace.
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_var.m (which is
 * MATLAB's var, denominator n-1), cross-checked against
 * jar/src/main/java/jline/api/trace/Trace_var.java.
 *
 * DIVERGENCE, MATLAB vs JAR: the JAR computes E[X^2] - E[X]^2, i.e. the
 * POPULATION variance with denominator n, while MATLAB's var uses the
 * unbiased denominator n-1. The two differ by the factor n/(n-1), which is
 * not negligible on the short traces used in tests and propagates into
 * trace_scv, trace_idi and trace_idc. MATLAB is the reference, so `unbiased`
 * defaults to true; pass false to reproduce the JAR.
 *
 * ARITHMETIC: sums and one division, exact in Rational.
 */

#include <vector>

#include "line/api/trace/trace_mean.h"
#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/**
 * @param S        the trace
 * @param unbiased true for the n-1 denominator (MATLAB var), false for the
 *                 n denominator (the JAR, and the lag-0 autocovariance)
 */
template <class T>
T trace_var(const std::vector<T>& S, bool unbiased = true) {
    detail::require_nonempty(S, "trace_var");
    const long n = static_cast<long>(S.size());
    if (unbiased && n < 2) throw InputError("trace_var: the unbiased variance needs n >= 2");
    const T mu = trace_mean(S);
    T s2 = num_traits<T>::from_int(0);
    for (const T& v : S) {
        const T d = v - mu;
        s2 += d * d;
    }
    return s2 / num_traits<T>::from_int(unbiased ? n - 1 : n);
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_VAR_H
