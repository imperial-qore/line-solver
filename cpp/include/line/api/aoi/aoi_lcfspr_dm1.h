/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSPR_DM1_H
#define LINE_API_AOI_LCFSPR_DM1_H

/**
 * Mean, variance and peak Age of Information of a D/M/1 preemptive LCFS queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfspr_dm1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfspr_dm1.java (identical).
 *
 *   E[A]     = tau + 1/mu
 *   q        = P(S < tau) = 1 - exp(-mu tau)
 *   E[S|succ]= (1/mu - tau exp(-mu tau) - exp(-mu tau)/mu) / q
 *   E[Apeak] = E[S|succ] + tau/q
 *   Var[A]   = 1/mu^2                    (the interarrival time is constant)
 *
 * from Inoue et al. (2019, Section IV).
 *
 * static_assert(num_traits<T>::has_transcendental) -- the exponentials in the
 * success probability and the conditional service time.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"

namespace line {
namespace aoi {

/**
 * @param tau deterministic interarrival time, > 0
 * @param mu  service rate, > 0
 * @return    [meanAoI, varAoI, peakAoI]
 */
template <class T>
AoiResult<T> aoi_lcfspr_dm1(const T& tau, const T& mu) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lcfspr_dm1 requires transcendental arithmetic");
    detail::require_positive(tau, "aoi_lcfspr_dm1", "the interarrival time tau");
    detail::require_positive(mu, "aoi_lcfspr_dm1", "the service rate mu");
    const T one = num_traits<T>::from_int(1);
    const T lambda = one / tau;
    const T rho = lambda / mu;
    detail::require_stable(rho, "aoi_lcfspr_dm1");

    const T meanAoI = tau + one / mu;
    const T e = detail::num_exp(T(-mu * tau));
    const T q = one - e;
    const T ES_succ = (one / mu - tau * e - e / mu) / q;
    const T peakAoI = ES_succ + tau / q;
    const T varAoI = one / (mu * mu);
    return {meanAoI, varAoI, peakAoI};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSPR_DM1_H
