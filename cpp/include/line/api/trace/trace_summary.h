/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_TRACE_TRACE_SUMMARY_H
#define LINE_API_TRACE_TRACE_SUMMARY_H

/**
 * Descriptive summary of a trace: moments, shape, order statistics,
 * autocorrelation and burstiness.
 *
 * Templated port of `jar/src/main/java/jline/api/trace/Trace_var.java#trace_summary`,
 * cross-checked against matlab/lib/kpctoolbox/trace/trace_summary.m. The
 * MATLAB version prints its result to a file id and returns the same
 * quantities as separate outputs; the JAR returns them packed in a vector.
 * This port returns a struct, so nothing is printed and nothing is positional.
 *
 * DIVERGENCES, MATLAB vs JAR, resolved as follows:
 *  - MAD: MATLAB uses mad(m,1), the median absolute deviation about the
 *    median. The JAR sorts the absolute deviations and takes element n/2,
 *    which is the upper median for even n and is not the median at all when
 *    n is even. MATLAB is the reference; the true median is used here.
 *  - KURT: the JAR subtracts 3, MATLAB's kurtosis does not. The field is
 *    named kurt_excess to make the convention explicit; add 3 for MATLAB.
 *  - SCV / IDC: the JAR uses the population variance everywhere, MATLAB the
 *    unbiased one; MATLAB is the reference (see trace_var.h).
 *  - percentiles: the linear-interpolation ("type 7") rule of the JAR is
 *    used. MATLAB's prctile interpolates on the (i-0.5)/n grid ("type 5")
 *    and gives different values on short traces; both agree in the limit.
 *
 * ARITHMETIC: the skewness and the standard deviation take square roots.
 *   static_assert(num_traits<T>::has_transcendental)
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/trace/trace_acf.h"
#include "line/api/trace/trace_idc.h"
#include "line/api/trace/trace_mean.h"
#include "line/api/trace/trace_scv.h"
#include "line/api/trace/trace_skew.h"
#include "line/api/trace/trace_types.h"
#include "line/api/trace/trace_var.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace trace {

namespace detail {

/** Median of an already sorted sample. */
template <class T>
inline T sorted_median(const std::vector<T>& x) {
    const std::size_t n = x.size();
    if (n % 2 == 1) return x[n / 2];
    return (x[n / 2 - 1] + x[n / 2]) / num_traits<T>::from_int(2);
}

/** Linear-interpolation percentile of an already sorted sample. */
template <class T>
inline T sorted_percentile(const std::vector<T>& x, long p_num, long p_den) {
    const long n = static_cast<long>(x.size());
    if (n == 1) return x[0];
    // index = (p/100) * (n-1), kept as an exact rational position
    const long num = p_num * (n - 1);
    const long lower = num / (p_den * 100);
    const long rem = num - lower * p_den * 100;
    if (rem == 0) return x[static_cast<std::size_t>(lower)];
    const T w = num_traits<T>::from_rational(rem, p_den * 100);
    return x[static_cast<std::size_t>(lower)] * (num_traits<T>::from_int(1) - w) +
           x[static_cast<std::size_t>(lower + 1)] * w;
}

}  // namespace detail

/** Return value of trace_summary. */
template <class T>
struct TraceSummary {
    T mean;
    T scv;
    T mad;           ///< median absolute deviation about the median
    T skew;          ///< bias-corrected skewness, see trace_skew
    T kurt_excess;   ///< population kurtosis minus 3
    T q25, q50, q75, p95;
    T min, max, iqr;
    std::vector<T> acf;  ///< lags 1..4
    T idc;
    T idc_scv_ratio;
};

/** @param S the trace, at least 6 samples so that the lag-4 acf exists. */
template <class T>
TraceSummary<T> trace_summary(const std::vector<T>& S) {
    static_assert(num_traits<T>::has_transcendental,
                  "trace_summary requires transcendental arithmetic");
    detail::require_nonempty(S, "trace_summary");
    if (S.size() < 6) throw InputError("trace_summary: at least six samples are required");

    TraceSummary<T> out;
    out.mean = trace_mean(S);
    out.scv = trace_scv(S);

    std::vector<T> x = S;
    std::sort(x.begin(), x.end());
    out.min = x.front();
    out.max = x.back();
    out.q25 = detail::sorted_percentile(x, 25, 1);
    out.q50 = detail::sorted_percentile(x, 50, 1);
    out.q75 = detail::sorted_percentile(x, 75, 1);
    out.p95 = detail::sorted_percentile(x, 95, 1);
    out.iqr = out.q75 - out.q25;

    std::vector<T> dev(S.size());
    for (std::size_t i = 0; i < S.size(); ++i) dev[i] = num_abs(T(S[i] - out.q50));
    std::sort(dev.begin(), dev.end());
    out.mad = detail::sorted_median(dev);

    out.skew = trace_skew(S);

    const T varp = trace_var(S, false);
    if (varp == num_traits<T>::from_int(0))
        throw NumericError("trace_summary: the trace is constant");
    T k4 = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < S.size(); ++i) {
        const T d = S[i] - out.mean;
        k4 += d * d * d * d;
    }
    k4 /= num_traits<T>::from_int(static_cast<long>(S.size()));
    out.kurt_excess = k4 / (varp * varp) - num_traits<T>::from_int(3);

    std::vector<int> lags;
    for (int l = 1; l <= 4; ++l) lags.push_back(l);
    out.acf = trace_acf(S, lags);

    out.idc = trace_idc(S);
    out.idc_scv_ratio = out.idc / out.scv;
    return out;
}

}  // namespace trace
}  // namespace line

#endif  // LINE_API_TRACE_TRACE_SUMMARY_H
