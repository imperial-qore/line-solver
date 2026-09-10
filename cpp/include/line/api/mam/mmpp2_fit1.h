/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMPP2_FIT1_H
#define LINE_API_MAM_MMPP2_FIT1_H

/**
 * MAP(2) matching mean, SCV, skewness and the index of dispersion for counts
 * (matlab/lib/kpctoolbox/mmpp/mmpp2_fit1.m).
 *
 * The four characteristics are converted to the (e1, e2, e3, g2) coordinates
 * of map2_fit:
 *   E2 = (1 + scv) E1^2,
 *   g2 = -(scv - idc)/(idc - 1),
 *   E3 = -(2 E1^3 - 3 E1 E2 - skew (E2 - E1^2)^(3/2)),
 * with skew = -1 acting as the sentinel that lets map2_fit choose E3 itself.
 *
 * Gated on transcendental arithmetic: the third moment needs (E2 - E1^2)^(3/2)
 * and map2_fit is itself gated.
 */

#include "line/api/mam/map2_fit.h"
#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** MAP(2) with the given mean, SCV, skewness and IDC. */
template <class T>
Map2FitResult<T> mmpp2_fit1(const T& mean, const T& scv, const T& skew, const T& idc) {
    static_assert(num_traits<T>::has_transcendental,
                  "mmpp2_fit1 requires transcendental arithmetic");
    using fitdetail::num_sqrt;
    using fitdetail::pw;
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T E1 = mean;
    const T E2 = (one + scv) * E1 * E1;
    if (idc == one) throw InputError("mmpp2_fit1: unit IDC admits no MMPP(2)");
    const T g2 = -(scv - idc) / (idc - one);
    T E3;
    if (skew == -one) {
        E3 = -one;  // map2_fit sentinel: pick the third moment automatically
    } else {
        const T v = E2 - E1 * E1;
        if (v < num_traits<T>::from_int(0)) throw InputError("mmpp2_fit1: negative variance");
        E3 = -(two * pw(E1, 3) - three * E1 * E2 - skew * v * num_sqrt(v));
    }
    return map2_fit(E1, E2, E3, g2);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMPP2_FIT1_H
