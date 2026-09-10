/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSS_GIM1_H
#define LINE_API_AOI_LCFSS_GIM1_H

/**
 * Mean and peak Age of Information of a GI/M/1 non-preemptive LCFS queue with
 * set-aside (LCFS-S).
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfss_gim1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfss_gim1.java (identical).
 *
 *   sigma solves Y*(mu - mu sigma) = sigma in (0,1)
 *   E[D]     = 1/(mu (1 - sigma))
 *   E[A]     = E[Y] + 1/mu + sigma E[D]
 *   E[Apeak] = E[Y] + E[D]
 *
 * from Inoue et al. (2019, Section V) adapted to GI/M/1.
 *
 * static_assert(num_traits<T>::has_transcendental) -- sigma is a bracketed
 * root of a transcendental equation; the rest is rational in sigma.
 *
 * As in aoi_lcfsd_gim1, MATLAB's sigma = rho fallback on a bracketing failure
 * is replaced by an error. MATLAB returns an empty LST handle, reported here
 * as has_lst = false.
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
 * @return      [meanAoI, peakAoI], no LST
 */
template <class T>
AoiLstResult<T> aoi_lcfss_gim1(const Lst<T>& Y_lst, const T& mu, const T& E_Y) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lcfss_gim1 requires transcendental arithmetic");
    detail::require_positive(mu, "aoi_lcfss_gim1", "the service rate mu");
    detail::require_positive(E_Y, "aoi_lcfss_gim1", "the mean interarrival time E_Y");
    if (!Y_lst) throw InputError("aoi_lcfss_gim1: the interarrival LST must be callable");
    const T one = num_traits<T>::from_int(1);
    const T lambda = one / E_Y;
    const T rho = lambda / mu;
    detail::require_stable(rho, "aoi_lcfss_gim1");

    const T E_S = one / mu;
    const T sigma = detail::gim1_sigma<T>(Y_lst, mu, "aoi_lcfss_gim1");
    const T E_D = one / (mu * (one - sigma));
    const T meanAoI = E_Y + E_S + sigma * E_D;
    const T peakAoI = E_Y + E_D;
    return {meanAoI, peakAoI, Lst<T>(), false};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSS_GIM1_H
