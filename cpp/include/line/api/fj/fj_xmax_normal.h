/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_NORMAL_H
#define LINE_API_FJ_XMAX_NORMAL_H

/**
 * Expected maximum and variance of K i.i.d. normal samples.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_normal.m, cross-checked against
 * FJ_xmax.fj_xmax_normal in jar/src/main/java/jline/api/fj/FJ_xmax.java
 * (identical).
 *
 *   arnold:    G(K) = sqrt(2 ln K)
 *   johnson:   G(K) = sqrt(2 ln K) - (ln ln K - ln(4 pi) + 2 gamma)/(2 sqrt(2 ln K))
 *   corrected: johnson minus the Petzold bias 0.1727 K^-0.2750
 *   Var        ~ 1.64492 sigma^2 / (2 ln K)
 *
 * static_assert(num_traits<T>::has_transcendental) -- log, sqrt and a real
 * power throughout. Note ln ln K is -inf at K = 2 in MATLAB (ln 2 < 1 makes
 * ln ln K finite and negative; it is K = 1 that diverges), and the function
 * requires K >= 2 for that reason.
 */

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K      number of samples, K >= 2
 * @param mu     mean of the normal
 * @param sigma  standard deviation, >= 0
 * @param method which correction to apply
 * @return       [Xmax, Vmax]
 */
template <class T>
FJXmaxNormalResult<T> fj_xmax_normal(unsigned K, const T& mu, const T& sigma,
                                     FJNormalMethod method = FJNormalMethod::Johnson) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_xmax_normal requires transcendental arithmetic");
    if (K < 2) throw InputError("fj_xmax_normal: the normal approximation requires K >= 2");
    if (sigma < num_traits<T>::from_int(0))
        throw InputError("fj_xmax_normal: sigma must be non-negative");

    const T two = num_traits<T>::from_int(2);
    const T gamma_em = num_traits<T>::from_double(0.5772156649015329);
    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    const T lnK = detail::num_log(Kt);
    const T sqrt_2lnK = detail::num_sqrt(T(two * lnK));

    T GK = sqrt_2lnK;
    if (method != FJNormalMethod::Arnold) {
        const T corr = (detail::num_log(lnK) - detail::num_log(T(num_traits<T>::from_int(4) * detail::num_pi<T>())) +
                        two * gamma_em) /
                       (two * sqrt_2lnK);
        GK = sqrt_2lnK - corr;
        if (method == FJNormalMethod::Corrected)
            GK -= num_traits<T>::from_double(0.1727) *
                  detail::num_pow(Kt, T(num_traits<T>::from_double(-0.2750)));
    }
    const T Xmax = mu + sigma * GK;
    const T Vmax = num_traits<T>::from_double(1.64492) * sigma * sigma / (two * lnK);
    return {Xmax, Vmax};
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_NORMAL_H
