/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSPR_MM1_H
#define LINE_API_AOI_LCFSPR_MM1_H

/**
 * Mean, variance and peak Age of Information of an M/M/1 preemptive LCFS queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfspr_mm1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfspr_mm1.java (identical).
 *
 *   E[A]     = (1/mu)(1 + 1/rho) = 1/mu + 1/lambda
 *   E[Apeak] = 1/(lambda+mu) + 1/lambda + 1/mu
 *   E[A^2]   = 2 (1/lambda^2 + 1/(lambda mu) + 1/mu^2)
 *
 * from Inoue et al. (2019, Section IV). Rational throughout, hence exact in
 * the field. The identity worth checking is E[A]_LCFSPR < E[A]_FCFS at every
 * rho in (0,1): both sides are exact rationals, so the comparison is a
 * theorem rather than a numerical observation.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"

namespace line {
namespace aoi {

/**
 * @param lambda arrival rate, > 0
 * @param mu     service rate, > 0
 * @return       [meanAoI, varAoI, peakAoI]
 */
template <class T>
AoiResult<T> aoi_lcfspr_mm1(const T& lambda, const T& mu) {
    detail::require_positive(lambda, "aoi_lcfspr_mm1", "the arrival rate lambda");
    detail::require_positive(mu, "aoi_lcfspr_mm1", "the service rate mu");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    detail::require_stable(rho, "aoi_lcfspr_mm1");

    const T meanAoI = (one / mu) * (one + one / rho);
    const T peakAoI = one / (lambda + mu) + one / lambda + one / mu;
    const T E_A2 = two * (one / (lambda * lambda) + one / (lambda * mu) + one / (mu * mu));
    T varAoI = E_A2 - meanAoI * meanAoI;
    if (varAoI < num_traits<T>::from_int(0)) varAoI = num_traits<T>::from_int(0);
    return {meanAoI, varAoI, peakAoI};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSPR_MM1_H
