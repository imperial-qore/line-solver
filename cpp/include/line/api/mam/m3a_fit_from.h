/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_M3A_FIT_FROM_H
#define LINE_API_MAM_M3A_FIT_FROM_H

/**
 * The m3a fitters driven from a process or from a trace rather than from
 * moments.
 *
 * Templated port of matlab/lib/m3a/m3a/aph2/{aph2_fit_map,aph2_fit_trace}.m and
 * amap2/{amap2_fit_gamma_map,amap2_fit_gamma_trace}.m. Each is a four-line
 * wrapper that measures the descriptors and hands them to the moment fitter, so
 * what they really contribute is WHICH descriptors a given fitter consumes and
 * how they are measured -- and that is worth having in one place, because the
 * two families disagree:
 *
 *   aph2_*        (M1, M2, M3), no autocorrelation: an APH(2) is renewal.
 *   amap2_gamma_* (M1, M2, M3, gamma), where gamma is the autocorrelation decay
 *                 rate, measured by `map_gamma` on a process and by
 *                 `trace_gamma` on a trace.
 *
 * THE TWO GAMMAS ARE NOT THE SAME MEASUREMENT. `map_gamma` reads the decay rate
 * off the process's own autocorrelation; `trace_gamma` FITS one on a grid and
 * returns the residual of that fit alongside it. The trace entry point here
 * discards the residual, exactly as the reference does, but a caller that wants
 * to know how well the geometric decay described the trace should call
 * `trace_gamma` directly rather than infer it from the fit.
 *
 * The trace moments are the RAW sample moments mean(T), mean(T^2), mean(T^3),
 * not the unbiased estimators; that is the reference's choice and it matters at
 * small sample sizes.
 *
 * ARITHMETIC: transcendental, inherited from the fitters.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/amap2_fit_gamma.h"
#include "line/api/mam/aph2_fit.h"
#include "line/api/mam/map_gamma.h"
#include "line/api/mam/map_moment.h"
#include "line/api/trace/trace_gamma.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

namespace m3adetail {

/** The raw sample moments 1..3 of a trace, as the reference measures them. */
template <class T>
void trace_moments3(const std::vector<T>& S, T* m1, T* m2, T* m3) {
    if (S.empty()) throw InputError("m3a fit: the trace is empty");
    const T zero = num_traits<T>::from_int(0);
    T a = zero, b = zero, c = zero;
    for (std::size_t i = 0; i < S.size(); ++i) {
        const T x = S[i];
        a += x;
        b += x * x;
        c += x * x * x;
    }
    const T n = num_traits<T>::from_int(static_cast<long>(S.size()));
    *m1 = T(a / n);
    *m2 = T(b / n);
    *m3 = T(c / n);
}

}  // namespace m3adetail

/** Fit an APH(2) to the first three moments of a MAP. */
template <class T>
Aph2FitResult<T> aph2_fit_map(const Map<T>& m) {
    return aph2_fit(map_mean(m), map_moment(m, 2), map_moment(m, 3));
}

/** Fit an APH(2) to the first three sample moments of a trace. */
template <class T>
Aph2FitResult<T> aph2_fit_trace(const std::vector<T>& S) {
    T m1, m2, m3;
    m3adetail::trace_moments3(S, &m1, &m2, &m3);
    return aph2_fit(m1, m2, m3);
}

/** Fit an AMAP(2) to the three moments and the decay rate of a MAP. */
template <class T>
Amap2FitGammaResult<T> amap2_fit_gamma_map(const Map<T>& m) {
    return amap2_fit_gamma(map_mean(m), map_moment(m, 2), map_moment(m, 3), map_gamma(m));
}

/**
 * Fit an AMAP(2) to the three sample moments and the fitted decay rate of a
 * trace. The residual of the gamma fit is discarded, as in the reference.
 */
template <class T>
Amap2FitGammaResult<T> amap2_fit_gamma_trace(const std::vector<T>& S) {
    T m1, m2, m3;
    m3adetail::trace_moments3(S, &m1, &m2, &m3);
    return amap2_fit_gamma(m1, m2, m3, line::trace::trace_gamma(S).gamma);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_M3A_FIT_FROM_H
