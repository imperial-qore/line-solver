/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSPR_GIM1_H
#define LINE_API_AOI_LCFSPR_GIM1_H

/**
 * Mean Age of Information, its transform and the peak age of a GI/M/1
 * preemptive LCFS queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfspr_gim1.m, cross-checked
 * against jar/src/main/java/jline/api/aoi/Aoi_lcfspr_gim1.java (identical).
 *
 *   E[A]      = E[Y] + 1/mu
 *   q         = P(S < Y) = 1 - Y*(mu)
 *   E[S|succ] = (1/mu + Y*'(mu) - Y*(mu)/mu) / q
 *   E[Apeak]  = E[S|succ] + 1/(lambda q)
 *   A*(s)     = Y*(s) mu/(s + mu)
 *
 * from Inoue et al. (2019, Section IV).
 *
 * static_assert(num_traits<T>::has_transcendental) -- MATLAB's
 * finite-difference derivative of Y* at mu. The mean is rational in E_Y and mu
 * and is exact; the transform is exact wherever Y* is.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace aoi {

/**
 * @param Y_lst LST of the interarrival time
 * @param mu    service rate, > 0
 * @param E_Y   mean interarrival time, > 0
 * @return      [meanAoI, peakAoI, A*(s)]
 */
template <class T>
AoiLstResult<T> aoi_lcfspr_gim1(const Lst<T>& Y_lst, const T& mu, const T& E_Y) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lcfspr_gim1 requires transcendental arithmetic");
    detail::require_positive(mu, "aoi_lcfspr_gim1", "the service rate mu");
    detail::require_positive(E_Y, "aoi_lcfspr_gim1", "the mean interarrival time E_Y");
    if (!Y_lst) throw InputError("aoi_lcfspr_gim1: the interarrival LST must be callable");
    const T one = num_traits<T>::from_int(1);
    const T lambda = one / E_Y;
    const T rho = lambda / mu;
    detail::require_stable(rho, "aoi_lcfspr_gim1");

    const T meanAoI = E_Y + one / mu;
    const T dYstar = detail::lst_derivative<T>(Y_lst, mu);
    const T q = one - Y_lst(mu);
    if (q == num_traits<T>::from_int(0))
        throw NumericError("aoi_lcfspr_gim1: Y*(mu) = 1, no update ever completes");
    const T ES_succ = (one / mu + dYstar - Y_lst(mu) / mu) / q;
    const T peakAoI = ES_succ + one / (lambda * q);

    Lst<T> lstAoI = [Y_lst, mu](const T& s) { return T(Y_lst(s) * (mu / (s + mu))); };
    return {meanAoI, peakAoI, lstAoI, true};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSPR_GIM1_H
