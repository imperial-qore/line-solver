/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_APH2_FITALL_H
#define LINE_API_MAM_APH2_FITALL_H

/**
 * All APH(2) representations matching three moments
 * (matlab/lib/m3a/m3a/aph2/aph2_fitall.m).
 *
 * The two phase means solve a quadratic whose discriminant is
 *   tmp0 = M3^2/9 + (8 M1^3/3 - 2 M1 M2) M3 - 3 M1^2 M2^2 + 2 M2^3,
 * giving one solution when tmp0 vanishes (identical diagonal entries of D0)
 * and two otherwise. Each is retained only when both phase means are positive
 * and the branching probability lies in [0, 1] up to a degeneracy tolerance.
 *
 * As in the in-tree MATLAB (the "added by GC" branches), an infeasible
 * discriminant or an empty feasible set falls back to a single aph_fit(M1, M2,
 * M3, 2), so the result is never empty.
 *
 * Gated on transcendental arithmetic: the discriminant square root and the
 * SCV <= 1 lower bound of the third moment.
 */

#include <vector>

#include "line/api/mam/aph2_assemble.h"
#include "line/api/mam/aph_fit.h"
#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/**
 * All feasible APH(2) fits of (M1, M2, M3). degentol is the tolerance used
 * both for the "M3 sits on its lower bound" degeneracy and for accepting a
 * branching probability marginally outside [0, 1] (MATLAB uses 1e-8).
 */
template <class T>
std::vector<Map<T>> aph2_fitall(const T& M1, const T& M2, const T& M3, const T& degentol) {
    static_assert(num_traits<T>::has_transcendental,
                  "aph2_fitall requires transcendental arithmetic");
    using fitdetail::num_sqrt;
    using fitdetail::pw;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);

    if (M1 <= zero) throw InputError("aph2_fitall: non-positive first moment");

    const T SCV = (M2 - M1 * M1) / (M1 * M1);
    bool degenerate = false;
    if (SCV <= one) {
        const T d = one - SCV;
        const T M3lb = three * pw(M1, 3) * (three * SCV - one + num_sqrt(two) * d * num_sqrt(d));
        if (num_abs(T(M3 - M3lb)) < degentol) degenerate = true;
    }

    T tmp0 = zero;
    if (!degenerate) {
        tmp0 = M3 * M3 / num_traits<T>::from_int(9) +
               ((num_traits<T>::from_int(8) * pw(M1, 3)) / three - two * M2 * M1) * M3 -
               three * M1 * M1 * M2 * M2 + two * pw(M2, 3);
        if (tmp0 < zero) {
            std::vector<Map<T>> out;
            out.push_back(aph_fit(M1, M2, M3, 2u).aph);
            return out;
        }
    }

    const T tmp1 = three * num_sqrt(tmp0);
    const T tmp2 = M3 - three * M1 * M2;
    const T tmp3 = num_traits<T>::from_int(6) * M2 - num_traits<T>::from_int(12) * M1 * M1;
    // M2 == 2 M1^2 is SCV == 1: the EXPONENTIAL, and the commonest input there
    // is, not a malformed moment set. The reference does not guard it -- at
    // SCV <= 1 with M3 on its lower bound it takes the tmp0 == 0 path and
    // evaluates tmp2/tmp3, which for an exponential is 0/0 and yields NaN that
    // the caller later discards as unfeasible. Throwing was worse than the NaN,
    // since it took down a Poisson split/merge that has an exact answer; but
    // returning the NaN would be worse still. Both roots of this quadratic are
    // unusable when tmp3 vanishes (0/0 when tmp2 does too, otherwise infinite),
    // so hand the moments to the general n-phase fitter, exactly as the
    // infeasible branch above already does.
    if (tmp3 == zero) {
        std::vector<Map<T>> out;
        out.push_back(aph_fit(M1, M2, M3, 2u).aph);
        return out;
    }

    const std::size_t n = (tmp0 == zero) ? 1u : 2u;
    std::vector<T> h1v(n, zero), h2v(n, zero);
    if (n == 1) {
        h2v[0] = tmp2 / tmp3;
        h1v[0] = h2v[0];
    } else {
        h2v[0] = (tmp2 + tmp1) / tmp3;
        h2v[1] = (tmp2 - tmp1) / tmp3;
        h1v[1] = h2v[0];
        h1v[0] = h2v[1];
    }

    std::vector<Map<T>> out;
    for (std::size_t j = 0; j < n; ++j) {
        const T h1 = h1v[j];
        const T h2 = h2v[j];
        if (h2 == zero) continue;
        T r1 = (M1 - h1) / h2;
        if (h1 > zero && h2 > zero && r1 >= -degentol && r1 <= one + degentol) {
            if (r1 > one) r1 = one;
            if (r1 < zero) r1 = zero;
            out.push_back(aph2_assemble(h1, h2, r1));
        }
    }
    if (out.empty()) out.push_back(aph_fit(M1, M2, M3, 2u).aph);
    return out;
}

/** aph2_fitall with the MATLAB default degentol = 1e-8. */
template <class T>
std::vector<Map<T>> aph2_fitall(const T& M1, const T& M2, const T& M3) {
    return aph2_fitall(M1, M2, M3, T(num_traits<T>::from_double(1e-8)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_APH2_FITALL_H
