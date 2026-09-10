/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSPR_MD1_H
#define LINE_API_AOI_LCFSPR_MD1_H

/**
 * Mean, variance and peak Age of Information of an M/D/1 preemptive LCFS queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfspr_md1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfspr_md1.java (identical).
 *
 *   E[A]     = 1/lambda + d
 *   E[Apeak] = d + exp(lambda d)/lambda
 *   Var[A]   = 1/lambda^2                       (the service time is constant)
 *
 * from Inoue et al. (2019, Section IV). Under preemption an update is
 * delivered only if no arrival occurs during its service, which happens with
 * probability exp(-lambda d); the reciprocal of that probability is the
 * exp(lambda d) in the peak age.
 *
 * static_assert(num_traits<T>::has_transcendental) -- the exp in the peak age.
 * The mean and the variance are rational and could be evaluated exactly, but
 * the triple is returned as a unit, as in MATLAB.
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
AoiResult<T> aoi_lcfspr_md1(const T& lambda, const T& d) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lcfspr_md1 requires transcendental arithmetic");
    detail::require_positive(lambda, "aoi_lcfspr_md1", "the arrival rate lambda");
    detail::require_positive(d, "aoi_lcfspr_md1", "the service time d");
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda * d;
    detail::require_stable(rho, "aoi_lcfspr_md1");

    const T meanAoI = one / lambda + d;
    const T peakAoI = d + detail::num_exp(T(lambda * d)) / lambda;
    const T varAoI = one / (lambda * lambda);
    return {meanAoI, varAoI, peakAoI};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSPR_MD1_H
