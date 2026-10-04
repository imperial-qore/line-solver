/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSD_GIM1_H
#define LINE_API_AOI_LCFSD_GIM1_H

/**
 * @file
 * @ingroup api_aoi
 * Mean Age of Information, its transform and the peak age of a GI/M/1
 * non-preemptive LCFS queue with discarding (LCFS-D).
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfsd_gim1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfsd_gim1.java (identical).
 *
 * Exact, Inoue et al. (2019, Section 3.3, NP-LCFS (C)), G* the interarrival
 * LST and rho = 1/(E_Y mu):
 *
 *   E[A]     = 1/mu + E[G^2]/(2E[G])
 *              + rho (-G*'(mu) + mu G*(mu) G*''(mu)/(1 + mu G*'(mu)))      (eq. 69)
 *   E[Apeak] = P0 (E[G] + (1+G*(mu))/mu) + Pw (E[G]/(1-G*(mu)) + 1/mu)
 *   P0       = q/(q + G*(mu)), Pw = 1 - P0, q = 1 + mu G*'(mu)/(1-G*(mu))   (eq. 91)
 *   A*(s)    = eq. (65)
 *
 * Discarding keeps at most one update waiting, so rho >= 1 is accepted.
 *
 * static_assert(num_traits<T>::has_transcendental) -- finite-difference first
 * and second derivatives of G* at mu, with the steps MATLAB uses.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace aoi {

/**
 * @brief Mean Age of Information, its transform and the peak age of a GI/M/1
 *        non-preemptive LCFS queue with discarding (LCFS-D).
 *
 * @param Y_lst LST of the interarrival time
 * @param mu    service rate, > 0
 * @param E_Y   mean interarrival time, > 0
 * @param E_Y2  second moment of the interarrival time, >= E_Y^2
 * @return      [meanAoI, peakAoI, A*(s)]
 */
template <class T>
AoiLstResult<T> aoi_lcfsd_gim1(const Lst<T>& Y_lst, const T& mu, const T& E_Y, const T& E_Y2) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lcfsd_gim1 requires transcendental arithmetic");
    detail::require_positive(mu, "aoi_lcfsd_gim1", "the service rate mu");
    detail::require_positive(E_Y, "aoi_lcfsd_gim1", "the mean interarrival time E_Y");
    if (!Y_lst) throw InputError("aoi_lcfsd_gim1: the interarrival LST must be callable");
    if (E_Y2 < E_Y * E_Y) throw InputError("aoi_lcfsd_gim1: E_Y2 must be at least E_Y^2");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T lambda = one / E_Y;
    const T rho = lambda / mu;

    const T gM = Y_lst(mu);
    const T mdY = T(-detail::lst_derivative<T>(Y_lst, mu));
    const T d2Y = detail::lst_second_derivative<T>(Y_lst, mu);

    const T meanAoI =
        one / mu + E_Y2 / (two * E_Y) + rho * (mdY + mu * gM * d2Y / (one - mu * mdY));

    // two-state chain "no wait"/"wait" of the informative updates (eq. 91)
    const T q = one - mu * mdY / (one - gM);
    const T P0 = q / (q + gM);
    const T Pw = gM / (q + gM);
    const T peakAoI = P0 * (E_Y + (one + gM) / mu) + Pw * (E_Y / (one - gM) + one / mu);

    Lst<T> lstAoI = [Y_lst, mu, E_Y, rho, gM, mdY, one](const T& s) {
        if (num_traits<T>::to_double(s < T(0) ? T(-s) : s) < 1e-12)
            return one;  // removable singularity of the residual-interarrival LST
        const T ys = s + mu;
        const T mdYs = T(-detail::lst_derivative<T>(Y_lst, ys));
        const T gres = (one - Y_lst(s)) / (s * E_Y);
        return T((gres + rho * mu / ys * (Y_lst(ys) - gM * (one - mu * mdYs) / (one - mu * mdY))) *
                 mu / ys);
    };
    return {meanAoI, peakAoI, lstAoI, true};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSD_GIM1_H
