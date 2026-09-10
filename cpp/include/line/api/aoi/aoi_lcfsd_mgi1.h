/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSD_MGI1_H
#define LINE_API_AOI_LCFSD_MGI1_H

/**
 * Mean and peak Age of Information of an M/GI/1 non-preemptive LCFS queue with
 * discarding (LCFS-D, equivalently M/GI/1/2*).
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfsd_mgi1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfsd_mgi1.java (identical).
 *
 *   E[H_res] = E[H^2] / (2 E[H])              (residual service time)
 *   E[T_eff] = E[H] + rho E[H_res]
 *   E[A]     = 1/lambda + E[H] + rho E[H^2]/(2 E[H]) + rho E[H]/(1+rho)
 *   E[Apeak] = 1/lambda + E[T_eff]
 *
 * from Inoue et al. (2019, Section VI).
 *
 * As with LCFS-S, the service-time LST is not used: only the first two moments
 * enter, so the port takes them directly. Rational throughout, hence exact in
 * the field. MATLAB returns an empty LST handle, reported here as
 * has_lst = false.
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
AoiLstResult<T> aoi_lcfsd_mgi1(const T& lambda, const T& E_H, const T& E_H2) {
    detail::require_positive(lambda, "aoi_lcfsd_mgi1", "the arrival rate lambda");
    detail::require_positive(E_H, "aoi_lcfsd_mgi1", "the mean service time E_H");
    if (E_H2 < E_H * E_H) throw InputError("aoi_lcfsd_mgi1: E_H2 must be at least E_H^2");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T rho = lambda * E_H;
    detail::require_stable(rho, "aoi_lcfsd_mgi1");

    const T E_Y = one / lambda;
    const T E_H_res = E_H2 / (two * E_H);
    const T E_T_eff = E_H + rho * E_H_res;

    const T meanAoI = E_Y + E_H + rho * E_H2 / (two * E_H) + rho * E_H / (one + rho);
    const T peakAoI = E_Y + E_T_eff;
    return {meanAoI, peakAoI, Lst<T>(), false};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSD_MGI1_H
