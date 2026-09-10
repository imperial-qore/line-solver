/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_FCFS_DM1_H
#define LINE_API_AOI_FCFS_DM1_H

/**
 * Mean, variance and peak Age of Information of a D/M/1 FCFS queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_fcfs_dm1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_fcfs_dm1.java (identical; both solve for
 * sigma by a bracketed root search, MATLAB with fzero on [0.001, 0.999] and
 * the JAR with bisection, which this port follows).
 *
 *   sigma solves  sigma = exp(-mu tau (1 - sigma))   in (0,1)
 *   E[D]     = 1/(mu (1 - sigma))
 *   E[A]     = tau/2 + E[D]
 *   E[Apeak] = tau + E[D]
 *   Var[A]   = (E[D^2] - E[D]^2) + (sigma/(mu(1-sigma)))^2,  E[D^2] = 2 E[D]^2
 *
 * from Inoue et al. (2019). With deterministic interarrivals the correlation
 * term in the mean vanishes, which is why E[A] is simply E[Y^2]/(2E[Y]) + E[D].
 *
 * static_assert(num_traits<T>::has_transcendental) -- sigma is the root of a
 * transcendental equation, so nothing downstream of it is in the field.
 *
 * The variance, as in the MATLAB source, is not the exact second moment of the
 * AoI: it adds the system-delay variance to a squared sigma term rather than
 * differentiating the AoI transform.
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
AoiResult<T> aoi_fcfs_dm1(const T& tau, const T& mu) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_fcfs_dm1 requires transcendental arithmetic");
    detail::require_positive(tau, "aoi_fcfs_dm1", "the interarrival time tau");
    detail::require_positive(mu, "aoi_fcfs_dm1", "the service rate mu");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T lambda = one / tau;
    const T rho = lambda / mu;
    detail::require_stable(rho, "aoi_fcfs_dm1");

    const T sigma = detail::bisect<T>(
        [&](const T& s) { return T(detail::num_exp(T(-mu * tau * (one - s))) - s); },
        num_traits<T>::from_double(0.001), num_traits<T>::from_double(0.999), "aoi_fcfs_dm1");

    const T E_D = one / (mu * (one - sigma));
    const T meanAoI = tau / two + E_D;
    const T peakAoI = tau + E_D;

    const T E_D2 = two * E_D * E_D;
    const T Var_D = E_D2 - E_D * E_D;
    const T extra = sigma / (mu * (one - sigma));
    T varAoI = Var_D + extra * extra;
    if (varAoI < num_traits<T>::from_int(0)) varAoI = num_traits<T>::from_int(0);
    return {meanAoI, varAoI, peakAoI};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_FCFS_DM1_H
