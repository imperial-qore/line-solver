/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_RMAX_EVD_H
#define LINE_API_FJ_RMAX_EVD_H

/**
 * Extreme-value approximation to the maximum of K branch response times, from
 * their mean and standard deviation.
 *
 * Templated port of matlab/src/api/fj/fj_rmax_evd.m, cross-checked against
 * FJ_rmax.fj_rmax_evd in jar/src/main/java/jline/api/fj/FJ_rmax.java
 * (identical).
 *
 *   Rmax = R + (sqrt(6) ln K / pi) sigma_R
 *
 * with the correction term divided by 1.27 in the calibrated variant of
 * Thomasian et al. (2007).
 *
 * static_assert(num_traits<T>::has_transcendental) -- sqrt and log. Note that
 * at K = 1 the correction vanishes and Rmax = R exactly, which is the only
 * value of K at which the approximation is exact.
 */

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K          number of branches, K >= 1
 * @param R          mean branch response time, > 0
 * @param sigma_R    standard deviation of the branch response time, >= 0
 * @param calibrated apply the 1/1.27 calibration of Thomasian et al. (2007)
 */
template <class T>
T fj_rmax_evd(unsigned K, const T& R, const T& sigma_R, bool calibrated = false) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_rmax_evd requires transcendental arithmetic");
    detail::require_positive_K(K, "fj_rmax_evd");
    if (R <= num_traits<T>::from_int(0)) throw InputError("fj_rmax_evd: the mean response time R must be positive");
    if (sigma_R < num_traits<T>::from_int(0)) throw InputError("fj_rmax_evd: sigma_R must be non-negative");

    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    T corr = detail::num_sqrt(T(num_traits<T>::from_int(6))) * detail::num_log(Kt) / detail::num_pi<T>();
    if (calibrated) corr /= num_traits<T>::from_double(1.27);
    return R + corr * sigma_R;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_RMAX_EVD_H
