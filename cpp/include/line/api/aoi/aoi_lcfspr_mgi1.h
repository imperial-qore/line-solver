/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSPR_MGI1_H
#define LINE_API_AOI_LCFSPR_MGI1_H

/**
 * @file
 * @ingroup api_aoi
 * Mean Age of Information, its transform and the peak age of an M/GI/1
 * preemptive LCFS queue.
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfspr_mgi1.m, cross-checked
 * against jar/src/main/java/jline/api/aoi/Aoi_lcfspr_mgi1.java (identical).
 *
 *   E[A]     = 1/(lambda H*(lambda))
 *   q        = P(S < Y) = H*(lambda)
 *   E[Apeak] = -H*'(lambda)/H*(lambda) + 1/(lambda H*(lambda))
 *   A*(s)    = lambda H*(s+lambda) / (s + lambda H*(s+lambda))
 *
 * Looking back from any instant, each earlier update is delivered iff its
 * service ends before the next arrival, so the age is a geometric sum of
 * Exp(lambda) gaps ending at the first delivered update. The old 1/lambda + E[H]
 * and lambda/(s+lambda) H*(s) hold only for exponential service.
 *
 * static_assert(num_traits<T>::has_transcendental) -- the peak age uses
 * MATLAB's finite-difference derivative of H* at lambda, with step
 * 1e-6 max(1,lambda). The mean and the transform are exact wherever H* is.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace aoi {

/**
 * @brief Mean Age of Information, its transform and the peak age of an M/GI/1
 *        preemptive LCFS queue.
 *
 * @param lambda arrival rate, > 0
 * @param H_lst  LST of the service time
 * @param E_H    mean service time, > 0
 * @return       [meanAoI, peakAoI, A*(s)]
 */
template <class T>
AoiLstResult<T> aoi_lcfspr_mgi1(const T& lambda, const Lst<T>& H_lst, const T& E_H) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lcfspr_mgi1 requires transcendental arithmetic");
    detail::require_positive(lambda, "aoi_lcfspr_mgi1", "the arrival rate lambda");
    detail::require_positive(E_H, "aoi_lcfspr_mgi1", "the mean service time E_H");
    if (!H_lst) throw InputError("aoi_lcfspr_mgi1: the service-time LST must be callable");
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda * E_H;
    detail::require_stable(rho, "aoi_lcfspr_mgi1");

    const T dHstar = detail::lst_derivative<T>(H_lst, lambda);
    const T Hlam = H_lst(lambda);
    if (Hlam == num_traits<T>::from_int(0))
        throw NumericError("aoi_lcfspr_mgi1: H*(lambda) = 0, no update ever completes");
    const T meanAoI = one / (lambda * Hlam);
    const T peakAoI = -dHstar / Hlam + one / (lambda * Hlam);

    Lst<T> lstAoI = [lambda, H_lst](const T& s) {
        const T f = lambda * H_lst(T(s + lambda));
        return T(f / (s + f));  // A*(0) = 1 identically
    };
    return {meanAoI, peakAoI, lstAoI, true};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSPR_MGI1_H
