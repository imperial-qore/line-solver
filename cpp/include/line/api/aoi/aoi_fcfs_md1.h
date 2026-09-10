/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_FCFS_MD1_H
#define LINE_API_AOI_FCFS_MD1_H

/**
 * Mean, variance and peak Age of Information of an M/D/1 FCFS queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_fcfs_md1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_fcfs_md1.java (identical).
 *
 *   rho = lambda d,  E[W] = lambda d^2 / (2(1-rho)),  E[T] = E[W] + d
 *   E[A]     = d (1/2 + 1/(2(1-rho)) + ((1-rho)/rho) exp(rho))
 *   E[Apeak] = E[T] + 1/lambda
 *   Var[A]   = 1/lambda^2 + 2 E[W] d/(1-rho) + d^2 rho/(1-rho)^2
 *
 * from Inoue et al. (2019). The exp(rho) in the mean is what carries the
 * negative correlation between the interarrival time and the waiting time.
 *
 * static_assert(num_traits<T>::has_transcendental) -- the single exp(rho).
 * Everything else in the file is rational; a caller who wants the exact
 * peak AoI and variance can get them from the M/GI/1 route with an Erlang
 * approximating the constant.
 *
 * The variance is documented in the MATLAB source as an approximation, not the
 * exact second moment: it is assembled from E[Y^2] - E[Y]^2 plus two waiting-
 * time terms, not from the AoI transform. Treat it as such.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"

namespace line {
namespace aoi {

/**
 * @param lambda arrival rate, > 0
 * @param d      deterministic service time, > 0
 * @return       [meanAoI, varAoI, peakAoI]
 */
template <class T>
AoiResult<T> aoi_fcfs_md1(const T& lambda, const T& d) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_fcfs_md1 requires transcendental arithmetic");
    detail::require_positive(lambda, "aoi_fcfs_md1", "the arrival rate lambda");
    detail::require_positive(d, "aoi_fcfs_md1", "the service time d");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T rho = lambda * d;
    detail::require_stable(rho, "aoi_fcfs_md1");

    const T E_H = d, E_H2 = d * d;
    const T E_W = lambda * E_H2 / (two * (one - rho));
    const T E_T = E_W + E_H;
    const T E_Y = one / lambda;
    const T E_Y2 = two / (lambda * lambda);

    const T meanAoI = d * (num_traits<T>::from_rational(1, 2) + one / (two * (one - rho)) +
                           ((one - rho) / rho) * detail::num_exp(rho));
    const T peakAoI = E_T + E_Y;

    T varAoI = E_Y2 - E_Y * E_Y + two * E_W * E_H / (one - rho) +
               E_H2 * rho / ((one - rho) * (one - rho));
    if (varAoI < num_traits<T>::from_int(0)) varAoI = num_traits<T>::from_int(0);
    return {meanAoI, varAoI, peakAoI};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_FCFS_MD1_H
