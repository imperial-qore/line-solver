/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSS_MGI1_H
#define LINE_API_AOI_LCFSS_MGI1_H

/**
 * Mean and peak Age of Information of an M/GI/1 non-preemptive LCFS queue with
 * set-aside (LCFS-S).
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfss_mgi1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfss_mgi1.java (identical).
 *
 *   E[B]  = E[H]/(1-rho),   E[B^2] = E[H^2]/(1-rho)^3       (M/G/1 busy period)
 *   E[A]     = 1/lambda + E[H] + lambda E[H^2] / (2 (1-rho)^2)
 *   E[Apeak] = 1/lambda + E[H] + lambda E[B^2] / (2 E[B])
 *
 * from Inoue et al. (2019, Section V).
 *
 * The service-time LST is NOT used: the formula depends on the service
 * distribution only through its first two moments, so the port takes E_H and
 * E_H2 and nothing else. (The MATLAB and Java signatures both accept an LST
 * argument and both ignore it -- see the report note.) Everything is rational
 * in lambda, E_H and E_H2, so both means are exact in the field.
 *
 * MATLAB returns an empty LST handle for this discipline, which the port
 * reports as has_lst = false.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace aoi {

/**
 * @param lambda arrival rate, > 0
 * @param E_H    mean service time, > 0
 * @param E_H2   second raw moment of the service time, >= E_H^2
 * @return       [meanAoI, peakAoI], no LST
 */
template <class T>
AoiLstResult<T> aoi_lcfss_mgi1(const T& lambda, const T& E_H, const T& E_H2) {
    detail::require_positive(lambda, "aoi_lcfss_mgi1", "the arrival rate lambda");
    detail::require_positive(E_H, "aoi_lcfss_mgi1", "the mean service time E_H");
    if (E_H2 < E_H * E_H) throw InputError("aoi_lcfss_mgi1: E_H2 must be at least E_H^2");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T rho = lambda * E_H;
    detail::require_stable(rho, "aoi_lcfss_mgi1");

    const T E_Y = one / lambda;
    const T E_B = E_H / (one - rho);
    const T E_B2 = E_H2 / num_pow_int(T(one - rho), 3);

    const T meanAoI = E_Y + E_H + lambda * E_H2 / (two * (one - rho) * (one - rho));
    const T peakAoI = E_Y + E_H + lambda * E_B2 / (two * E_B);
    return {meanAoI, peakAoI, Lst<T>(), false};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSS_MGI1_H
