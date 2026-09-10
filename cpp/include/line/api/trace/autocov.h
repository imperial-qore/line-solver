/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_AUTOCOV_H
#define LINE_API_TRACE_AUTOCOV_H

/**
 * Sample autocovariance sequence of a trace, lags 0 .. n-2.
 *
 * Templated port of matlab/lib/kpctoolbox/contrib/autocov.m, cross-checked
 * against the private autocov() of
 * jar/src/main/java/jline/api/trace/Trace_var.java (identical formula; the
 * MATLAB implementation evaluates it by FFT, which is the same quantity up to
 * rounding but is not usable in exact arithmetic, so this port evaluates the
 * defining sum directly).
 *
 *   acv(p) = 1/(n-p) sum_{i=1}^{n-p} (S_i - Sbar)(S_{i+p} - Sbar)
 *
 * The lag-0 term is therefore the POPULATION variance (denominator n), not
 * MATLAB's var: the header comment of autocov.m claiming acv(1) = var(X) is
 * wrong by the factor n/(n-1). See trace_var.h.
 *
 * ARITHMETIC: sums, products and divisions, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/trace_mean.h"
#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/** @return acv[0..n-2], acv[p] the lag-p sample autocovariance. */
template <class T>
std::vector<T> autocov(const std::vector<T>& S) {
    detail::require_nonempty(S, "autocov");
    const std::size_t n = S.size();
    if (n < 2) throw InputError("autocov: the trace must have at least two samples");
    const T mu = trace_mean(S);
    std::vector<T> X(n);
    for (std::size_t i = 0; i < n; ++i) X[i] = S[i] - mu;

    std::vector<T> acv(n - 1, num_traits<T>::from_int(0));
    for (std::size_t p = 0; p + 1 < n; ++p) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i + p < n; ++i) s += X[i] * X[i + p];
        acv[p] = s / num_traits<T>::from_int(static_cast<long>(n - p));
    }
    return acv;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_AUTOCOV_H
