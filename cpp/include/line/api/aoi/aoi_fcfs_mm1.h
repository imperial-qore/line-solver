/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_FCFS_MM1_H
#define LINE_API_AOI_FCFS_MM1_H

/**
 * Mean, variance and peak Age of Information of an M/M/1 FCFS queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_fcfs_mm1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_fcfs_mm1.java (identical).
 *
 *   E[A]     = (1/mu)(1 + 1/rho + rho^2/(1-rho))
 *   E[Apeak] = (1/mu)(1 + 1/rho + rho/(1-rho))
 *   E[A^2]   = (2/mu^2)(1 - rho - rho^3 + 4 rho^4 - 2 rho^5)/(rho^2 (1-rho)^2)
 *   Var[A]   = E[A^2] - E[A]^2
 *
 * from Inoue, Masuyama, Takine and Tanaka (IEEE Trans. IT 65(12), 2019).
 * Every quantity is a rational function of rho, so the triple is exact in the
 * field. This is the reference AoI closed form: E[A] has an interior minimum
 * in rho near 0.53, and the exact instantiation locates it as the root of a
 * polynomial rather than by a rounded search.
 *
 * MATLAB clamps a negative variance to zero as a numerical safety net; the
 * port keeps the clamp so the two agree, but note that in exact arithmetic the
 * clamp can never fire, because E[A^2] - E[A]^2 is then evaluated without
 * cancellation error.
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
AoiResult<T> aoi_fcfs_mm1(const T& lambda, const T& mu) {
    detail::require_positive(lambda, "aoi_fcfs_mm1", "the arrival rate lambda");
    detail::require_positive(mu, "aoi_fcfs_mm1", "the service rate mu");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    detail::require_stable(rho, "aoi_fcfs_mm1");

    const T meanAoI = (one / mu) * (one + one / rho + rho * rho / (one - rho));
    const T peakAoI = (one / mu) * (one + one / rho + rho / (one - rho));

    const T r2 = rho * rho, r3 = r2 * rho, r4 = r3 * rho, r5 = r4 * rho;
    const T E_A2 = (two / (mu * mu)) * (one - rho - r3 + num_traits<T>::from_int(4) * r4 - two * r5) /
                   (r2 * (one - rho) * (one - rho));
    T varAoI = E_A2 - meanAoI * meanAoI;
    if (varAoI < num_traits<T>::from_int(0)) varAoI = num_traits<T>::from_int(0);
    return {meanAoI, varAoI, peakAoI};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_FCFS_MM1_H
