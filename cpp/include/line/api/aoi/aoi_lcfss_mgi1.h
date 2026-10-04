/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LCFSS_MGI1_H
#define LINE_API_AOI_LCFSS_MGI1_H

/**
 * @file
 * @ingroup api_aoi
 * Mean Age of Information, its transform and the peak age of an M/GI/1
 * non-preemptive LCFS queue with set-aside (LCFS-S), i.e. LCFS without
 * discarding: older waiting updates are still served, after fresher ones.
 *
 * Templated port of matlab/src/api/aoi/aoi_lcfss_mgi1.m, cross-checked against
 * jar/src/main/java/jline/api/aoi/Aoi_lcfss_mgi1.java (identical).
 *
 * Exact, Inoue et al. (2019, Section 3.3, NP-LCFS (D)), H* the service LST:
 *
 *   E[A]     = lambda E[H^2]/2 + ((1-rho)^2/(rho H*(lambda)) + 2) E[H]     (eq. 70)
 *   E[Apeak] = E[H] + E[W] + 1/(lambda (2 - rho - H*(lambda)))
 *   E[W]     = ((1 - H*(lambda))/lambda + H*'(lambda))/(2 - rho - H*(lambda))
 *   A*(s)    = eq. (66)
 *
 * static_assert(num_traits<T>::has_transcendental) -- MATLAB's
 * finite-difference derivative of H* at lambda.
 */

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace aoi {

/**
 * @brief Mean Age of Information, its transform and the peak age of an M/GI/1
 *        non-preemptive LCFS queue with set-aside (LCFS-S).
 *
 * @param lambda arrival rate, > 0
 * @param H_lst  LST of the service time
 * @param E_H    mean service time, > 0
 * @param E_H2   second raw moment of the service time, >= E_H^2
 * @return       [meanAoI, peakAoI, A*(s)]
 */
template <class T>
AoiLstResult<T> aoi_lcfss_mgi1(const T& lambda, const Lst<T>& H_lst, const T& E_H, const T& E_H2) {
    static_assert(num_traits<T>::has_transcendental,
                  "aoi_lcfss_mgi1 requires transcendental arithmetic");
    detail::require_positive(lambda, "aoi_lcfss_mgi1", "the arrival rate lambda");
    detail::require_positive(E_H, "aoi_lcfss_mgi1", "the mean service time E_H");
    if (!H_lst) throw InputError("aoi_lcfss_mgi1: the service-time LST must be callable");
    if (E_H2 < E_H * E_H) throw InputError("aoi_lcfss_mgi1: E_H2 must be at least E_H^2");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T rho = lambda * E_H;
    detail::require_stable(rho, "aoi_lcfss_mgi1");

    const T hL = H_lst(lambda);
    const T mdH = T(-detail::lst_derivative<T>(H_lst, lambda));

    const T meanAoI = lambda * E_H2 / two + ((one - rho) * (one - rho) / (rho * hL) + two) * E_H;
    // informative arrivals: lambda_dagger = lambda (2 - rho - H*(lambda)); E[W] from eq. (98)
    const T den = two - rho - hL;
    const T E_W = ((one - hL) / lambda - mdH) / den;
    const T peakAoI = E_H + E_W + one / (lambda * den);

    Lst<T> lstAoI = [lambda, H_lst, E_H, rho, one](const T& s) {
        if (num_traits<T>::to_double(s < T(0) ? T(-s) : s) < 1e-12)
            return one;  // removable singularity of the residual-service LST
        const T Hs = H_lst(s);
        const T HsL = H_lst(T(s + lambda));
        const T hres_s = (one - Hs) / (s * E_H);
        return T(lambda / (s + lambda) * Hs *
                 (rho * hres_s + (one - rho) * (s + lambda) * (one - Hs + HsL) / (s + lambda * HsL)));
    };
    return {meanAoI, peakAoI, lstAoI, true};
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LCFSS_MGI1_H
