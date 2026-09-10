/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSD_GIM1_H
#define LINE_API_AOI_LCFSD_GIM1_H

/**
 * Mean and peak Age of Information of a GI/M/1 non-preemptive LCFS queue with
 * discarding (LCFS-D, equivalently GI/M/1/2*).
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfsd_gim1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfsd_gim1.java (identical).
 *
 *   sigma solves Y*(mu - mu sigma) = sigma in (0,1)
 *   E[S]     = 1/mu,  E[T_eff] = E[S](1 + sigma)
 *   E[A]     = E[Y] + E[S](1 + sigma) + sigma E[S]/(1 + sigma)
 *   E[Apeak] = E[Y] + E[T_eff]
 *
 * from Inoue et al. (2019, Section VI) adapted to GI/M/1: an arriving update
 * finds the server busy with probability sigma, and by memorylessness the
 * remaining service is again exponential with rate mu.
 *
 * static_assert(num_traits<T>::has_transcendental) -- sigma is a bracketed
 * root of a transcendental equation. Everything downstream of sigma is
 * rational.
 *
 * MATLAB falls back to sigma = rho when fzero fails to bracket; this port
 * instead reports the failure, since a silent substitution of the utilization
 * for the busy probability changes the answer without saying so. MATLAB
 * returns an empty LST handle, reported here as has_lst = false.
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
AoiLstResult<T> aoi_lcfsd_gim1(const Lst<T>& Y_lst, const T& mu, const T& E_Y) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lcfsd_gim1 requires transcendental arithmetic");
    detail::require_positive(mu, "aoi_lcfsd_gim1", "the service rate mu");
    detail::require_positive(E_Y, "aoi_lcfsd_gim1", "the mean interarrival time E_Y");
    if (!Y_lst) throw InputError("aoi_lcfsd_gim1: the interarrival LST must be callable");
    const T one = num_traits<T>::from_int(1);
    const T lambda = one / E_Y;
    const T rho = lambda / mu;
    detail::require_stable(rho, "aoi_lcfsd_gim1");

    const T E_S = one / mu;
    const T sigma = detail::gim1_sigma<T>(Y_lst, mu, "aoi_lcfsd_gim1");
    const T E_T_eff = E_S + sigma * E_S;
    const T meanAoI = E_Y + E_S * (one + sigma) + sigma * E_S / (one + sigma);
    const T peakAoI = E_Y + E_T_eff;
    return {meanAoI, peakAoI, Lst<T>(), false};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSD_GIM1_H
