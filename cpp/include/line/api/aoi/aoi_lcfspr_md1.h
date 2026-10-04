/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSPR_MD1_H
#define LINE_API_AOI_LCFSPR_MD1_H

/**
 * @file
 * @ingroup api_aoi
 * Mean, variance and peak Age of Information of an M/D/1 preemptive LCFS queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfspr_md1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfspr_md1.java (identical).
 *
 *   E[A]     = exp(lambda d)/lambda
 *   E[Apeak] = d + exp(lambda d)/lambda
 *   Var[A]   = (exp(2 lambda d) - 2 lambda d exp(lambda d))/lambda^2
 *
 * from Inoue et al. (2019, Section IV). Under preemption an update is
 * delivered only if no arrival occurs during its service, which happens with
 * probability exp(-lambda d); the reciprocal of that probability is the
 * exp(lambda d) in the mean and peak age. The moments follow from the M/GI/1
 * preemptive transform A*(s) = f(s)/(s + f(s)), f(s) = lambda H*(s+lambda).
 *
 * static_assert(num_traits<T>::has_transcendental) -- the exp in all three.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"

namespace line {
namespace aoi {

/**
 * @brief Mean, variance and peak Age of Information of an M/D/1 preemptive
 *        LCFS queue.
 *
 * @param lambda arrival rate, > 0
 * @param d      deterministic service time, > 0
 * @return       [meanAoI, varAoI, peakAoI]
 */
template <class T>
AoiResult<T> aoi_lcfspr_md1(const T& lambda, const T& d) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lcfspr_md1 requires transcendental arithmetic");
    detail::require_positive(lambda, "aoi_lcfspr_md1", "the arrival rate lambda");
    detail::require_positive(d, "aoi_lcfspr_md1", "the service time d");
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda * d;
    detail::require_stable(rho, "aoi_lcfspr_md1");

    const T two = num_traits<T>::from_int(2);
    const T e = detail::num_exp(T(lambda * d));
    const T meanAoI = e / lambda;
    const T peakAoI = d + e / lambda;
    const T varAoI = (e * e - two * lambda * d * e) / (lambda * lambda);
    return {meanAoI, varAoI, peakAoI};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSPR_MD1_H
