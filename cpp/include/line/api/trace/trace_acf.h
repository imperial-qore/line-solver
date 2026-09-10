/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_ACF_H
#define LINE_API_TRACE_TRACE_ACF_H

/**
 * Autocorrelation coefficients of a trace at the requested lags.
 *
 *   rho(k) = acv(k) / acv(0)
 *
 * Templated port of matlab/lib/kpctoolbox/trace/trace_acf.m, cross-checked
 * against `jar/src/main/java/jline/api/trace/Trace_var.java#trace_acf`.
 *
 * DIVERGENCE INSIDE MATLAB: trace_acf.m has two branches. When the Signal
 * Processing Toolbox is present it uses xcorr(...,'coeff'), whose estimator
 * divides every lag by the same n (the biased estimator); otherwise it falls
 * back to autocov, which divides lag p by n-p. The two therefore differ by
 * the factor n/(n-p) at lag p -- the same trace gives different acf values
 * depending on which toolboxes are installed. Measured on the trace
 * 1,2,3,4,5: MATLAB with the toolbox returns 0.4 and -0.1 at lags 1 and 2,
 * while the JAR and the MATLAB fallback return 1/2 and -1/6. The JAR ports
 * the fallback branch, and so does this header -- it is the branch the two
 * codebases share, and the one whose lag-0 value is the sample variance.
 *
 * OUT-OF-RANGE LAGS: MATLAB clamps them to n-2 and then deletes the clamped
 * entries, the JAR filters out every lag outside (0, n-2]. Both silently
 * return fewer values than lags requested; this port follows the JAR filter
 * and documents it rather than throwing, so that trace_summary and
 * trace_gamma keep working on short traces.
 *
 * ARITHMETIC: a ratio of autocovariances, exact in Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/trace/autocov.h"
#include "line/api/trace/trace_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

/**
 * @param S    the trace
 * @param lags the lags to evaluate; entries outside (0, n-2] are dropped
 * @return one coefficient per surviving lag, in the order given
 */
template <class T>
std::vector<T> trace_acf(const std::vector<T>& S, const std::vector<int>& lags) {
    detail::require_nonempty(S, "trace_acf");
    const long n = static_cast<long>(S.size());
    std::vector<int> kept;
    for (int l : lags)
        if (l > 0 && l <= n - 2) kept.push_back(l);
    if (kept.empty()) return std::vector<T>();

    const std::vector<T> acv = autocov(S);
    if (acv[0] == num_traits<T>::from_int(0))
        throw NumericError("trace_acf: the trace is constant, the acf is undefined");
    std::vector<T> rho;
    rho.reserve(kept.size());
    for (int l : kept) rho.push_back(acv[static_cast<std::size_t>(l)] / acv[0]);
    return rho;
}

/** Lag-1 autocorrelation, the MATLAB and JAR default. */
template <class T>
std::vector<T> trace_acf(const std::vector<T>& S) {
    return trace_acf(S, std::vector<int>(1, 1));
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_ACF_H
