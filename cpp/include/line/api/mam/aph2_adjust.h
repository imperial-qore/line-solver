/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_APH2_ADJUST_H
#define LINE_API_MAM_APH2_ADJUST_H

/**
 * Nearest APH(2)-feasible values of the second and third moments
 * (matlab/lib/m3a/m3a/aph2/aph2_adjust.m, method 'simple').
 *
 * Applies the closed-form bounds of Telek and Heindl (2002): the SCV is
 * lifted to 1/2 when it falls below it, and the third moment is clamped into
 * [lb, ub] for SCV <= 1 and above lb (times 1 + tol) for SCV > 1.
 *
 * Only the 'simple' method is ported. The other four methods of the reference
 * ('opt_param', 'opt_param_gads', 'opt_char', 'opt_char_gads') minimize a
 * distance with fmincon or GlobalSearch and are therefore out of scope; they
 * are documented as skipped rather than approximated.
 *
 * Gated on transcendental arithmetic: the lower bound for SCV <= 1 contains
 * sqrt(2) and (1 - scv)^(3/2).
 */

#include "line/api/mam/map_fit_detail.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** Result of aph2_adjust. */
template <class T>
struct Aph2AdjustResult {
    T M2a;
    T M3a;
};

/**
 * Feasible (M2, M3) closest to the input, holding M1 fixed. tol is the
 * relative slack applied above the SCV > 1 lower bound (MATLAB uses 1e-4).
 */
template <class T>
Aph2AdjustResult<T> aph2_adjust(const T& M1, const T& M2, const T& M3, const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "aph2_adjust requires transcendental arithmetic");
    using fitdetail::num_sqrt;
    using fitdetail::pw;

    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T six = num_traits<T>::from_int(6);
    const T half = num_traits<T>::from_rational(1, 2);

    if (M1 <= num_traits<T>::from_int(0)) throw InputError("aph2_adjust: non-positive first moment");

    const T M1sq = M1 * M1;
    const T scv = (M2 - M1sq) / M1sq;

    Aph2AdjustResult<T> r;
    T scva;
    if (scv < half) {
        r.M2a = three / two * M1sq;
        scva = (r.M2a - M1sq) / M1sq;
    } else {
        r.M2a = M2;
        scva = scv;
    }

    if (scva <= one) {
        const T d = one - scva;
        const T lb = three * pw(M1, 3) * (three * scva - one + num_sqrt(two) * d * num_sqrt(d));
        const T ub = six * pw(M1, 3) * scva;
        if (M3 < lb)
            r.M3a = lb;
        else if (M3 > ub)
            r.M3a = ub;
        else
            r.M3a = M3;
    } else {
        const T lb = three / two * pw(M1, 3) * (one + scva) * (one + scva);
        r.M3a = (M3 <= lb) ? T(lb * (one + tol)) : M3;
    }
    return r;
}

/** aph2_adjust with the MATLAB default slack tol = 1e-4. */
template <class T>
Aph2AdjustResult<T> aph2_adjust(const T& M1, const T& M2, const T& M3) {
    return aph2_adjust(M1, M2, M3, T(num_traits<T>::from_double(1e-4)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_APH2_ADJUST_H
