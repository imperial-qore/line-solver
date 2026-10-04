/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSS_GIM1_H
#define LINE_API_AOI_LCFSS_GIM1_H

/**
 * @file
 * @ingroup api_aoi
 * Mean Age of Information, its transform and the peak age of a GI/M/1
 * non-preemptive LCFS queue with set-aside (LCFS-S), i.e. LCFS without
 * discarding: older waiting updates are still served, after fresher ones.
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfss_gim1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfss_gim1.java (identical).
 *
 * Exact, Inoue et al. (2019, Section 3.3, NP-LCFS (D)), G* the interarrival
 * LST, rho = 1/(E_Y mu) and gamma the root of G*(mu - mu x) = x:
 *
 *   E[A]     = 1/mu + E[G^2]/(2E[G]) + rho (-G*'(mu - mu gamma))            (eq. 71)
 *   E[Apeak] = P0 (E[G] + (1+G*(mu))/mu)
 *              + Pw (1/mu + (E[G] + G*'(mu) + (gamma - G*(mu))/(gamma mu))/(1 - G*(mu)))
 *   P0       = (1-gamma)/(1 - gamma G*(mu)), Pw = 1 - P0             (eqs. 100-103)
 *   A*(s)    = eq. (67)
 *
 * static_assert(num_traits<T>::has_transcendental) -- gamma is a bracketed
 * root of a transcendental equation and G*' is a finite difference. A
 * bracketing failure is an error, never a gamma = rho fallback.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace aoi {

/**
 * @brief Mean Age of Information, its transform and the peak age of a GI/M/1
 *        non-preemptive LCFS queue with set-aside (LCFS-S).
 *
 * @param Y_lst LST of the interarrival time
 * @param mu    service rate, > 0
 * @param E_Y   mean interarrival time, > 0
 * @param E_Y2  second moment of the interarrival time, >= E_Y^2
 * @return      [meanAoI, peakAoI, A*(s)]
 */
template <class T>
AoiLstResult<T> aoi_lcfss_gim1(const Lst<T>& Y_lst, const T& mu, const T& E_Y, const T& E_Y2) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lcfss_gim1 requires transcendental arithmetic");
    detail::require_positive(mu, "aoi_lcfss_gim1", "the service rate mu");
    detail::require_positive(E_Y, "aoi_lcfss_gim1", "the mean interarrival time E_Y");
    if (!Y_lst) throw InputError("aoi_lcfss_gim1: the interarrival LST must be callable");
    if (E_Y2 < E_Y * E_Y) throw InputError("aoi_lcfss_gim1: E_Y2 must be at least E_Y^2");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T lambda = one / E_Y;
    const T rho = lambda / mu;
    detail::require_stable(rho, "aoi_lcfss_gim1");

    const T gam = detail::gim1_sigma<T>(Y_lst, mu, "aoi_lcfss_gim1");
    const T gM = Y_lst(mu);
    const T mdY = T(-detail::lst_derivative<T>(Y_lst, mu));
    const T mdYg = T(-detail::lst_derivative<T>(Y_lst, T(mu - mu * gam)));

    const T meanAoI = one / mu + E_Y2 / (two * E_Y) + rho * mdYg;

    const T P0 = (one - gam) / (one - gam * gM);
    const T Pw = gam * (one - gM) / (one - gam * gM);
    const T E0 = E_Y + (one + gM) / mu;
    const T Ew = one / mu + (E_Y - mdY + (gam - gM) / (gam * mu)) / (one - gM);
    const T peakAoI = P0 * E0 + Pw * Ew;

    Lst<T> lstAoI = [Y_lst, mu, E_Y, rho, gam, one](const T& s) {
        if (num_traits<T>::to_double(s < T(0) ? T(-s) : s) < 1e-12)
            return one;  // removable singularity of the residual-interarrival LST
        const T gres = (one - Y_lst(s)) / (s * E_Y);
        return T((gres + rho * (Y_lst(T(s + mu - mu * gam)) - gam) * mu / (s + mu)) * mu / (s + mu));
    };
    return {meanAoI, peakAoI, lstAoI, true};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSS_GIM1_H
