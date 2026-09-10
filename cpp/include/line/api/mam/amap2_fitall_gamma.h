/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_AMAP2_FITALL_GAMMA_H
#define LINE_API_MAM_AMAP2_FITALL_GAMMA_H

/**
 * All AMAP(2) representations matching three moments and the autocorrelation
 * decay rate (matlab/lib/m3a/m3a/amap2/amap2_fitall_gamma.m).
 *
 * The two phase means h1, h2 come from the same quadratic as aph2_fitall. For
 * GAMMA >= 0 the first canonical form leaves a further quadratic in the
 * branching probability r2, with discriminant
 *   z = M1^2 G^2 + (2 M1 h1 + 2 M1 h2 - 4 h1 h2 - 2 M1^2) G
 *       + M1^2 - 2 M1 h1 - 2 M1 h2 + h1^2 + 2 h1 h2 + h2^2,
 * so up to four AMAP(2)s can match the same characteristics; for GAMMA < 0 the
 * second canonical form determines r1 and r2 uniquely per (h1, h2). Solutions
 * whose probabilities fall outside [0, 1] beyond r12tol are discarded, the
 * others are clamped into [0, 1].
 *
 * Unlike aph2_fitall this returns an empty vector when the characteristics are
 * infeasible: no approximate fitting is performed, exactly as in the reference.
 *
 * Gated on transcendental arithmetic: the two moment discriminants and the
 * SCV <= 1 lower bound of the third moment are square roots.
 */

#include <vector>

#include "line/api/mam/amap2_assemble.h"
#include "line/api/mam/map_fit_detail.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/**
 * Every AMAP(2) matching (M1, M2, M3, GAMMA). degentol screens the degenerate
 * discriminants (MATLAB 1e-8); r12tol is the slack allowed on the branching
 * probabilities before a solution is rejected (MATLAB 1e-6).
 */
template <class T>
std::vector<Map<T>> amap2_fitall_gamma(const T& M1, const T& M2, const T& M3, const T& GAMMA,
                                       const T& degentol, const T& r12tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "amap2_fitall_gamma requires transcendental arithmetic");
    using fitdetail::num_sqrt;
    using fitdetail::pw;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T three = num_traits<T>::from_int(3);
    const T four = num_traits<T>::from_int(4);

    std::vector<Map<T>> out;
    if (M1 <= zero) throw InputError("amap2_fitall_gamma: non-positive first moment");

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
        if (tmp0 < zero) return out;
    }

    const T tmp1 = three * num_sqrt(tmp0);
    const T tmp2 = M3 - three * M1 * M2;
    const T tmp3 = num_traits<T>::from_int(6) * M2 - num_traits<T>::from_int(12) * M1 * M1;
    if (tmp3 == zero) throw NumericError("amap2_fitall_gamma: degenerate moment set (M2 = 2 M1^2)");

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
    for (std::size_t j = 0; j < n; ++j)
        if (h2v[j] <= zero) return out;

    const T lo = -r12tol;
    const T hi = one + r12tol;

    for (std::size_t j = 0; j < n; ++j) {
        const T h1 = h1v[j];
        const T h2 = h2v[j];
        if (GAMMA >= zero) {
            const T z = M1 * M1 * GAMMA * GAMMA +
                        (two * M1 * h1 + two * M1 * h2 - four * h1 * h2 - two * M1 * M1) * GAMMA +
                        M1 * M1 - two * M1 * h1 - two * M1 * h2 + h1 * h1 + two * h1 * h2 + h2 * h2;
            std::vector<T> r2v;
            if (num_abs(z) < degentol) {
                if (h1 == zero) continue;
                r2v.push_back(T((h1 - M1 + h2 + GAMMA * M1) / (two * h1)));
            } else if (z > zero) {
                if (h1 == zero) continue;
                const T s = num_sqrt(z);
                r2v.push_back(T((h1 - M1 + h2 - s + GAMMA * M1) / (two * h1)));
                r2v.push_back(T((h1 - M1 + h2 + s + GAMMA * M1) / (two * h1)));
            }
            for (std::size_t i = 0; i < r2v.size(); ++i) {
                T r2 = r2v[i];
                const T den = h2 - M1 * r2;
                if (den == zero) continue;
                T r1 = (M1 - h1 - M1 * r2 + h1 * r2) / den;
                if (!(r1 >= lo && r1 <= hi && r2 >= lo && r2 <= hi)) continue;
                if (r1 > one) r1 = one;
                if (r1 < zero) r1 = zero;
                if (r2 > one) r2 = one;
                if (r2 < zero) r2 = zero;
                out.push_back(amap2_assemble(h1, h2, r1, r2, 1));
            }
        } else {
            if (h1 == zero) continue;
            T r2 = (h1 - M1 + h2 + GAMMA * M1) / h1;
            if (r2 == one) continue;
            T r1 = (r2 + (h1 + h2 - h1 * r2) / M1 - two) / (r2 - one);
            if (!(r1 >= lo && r1 <= hi && r2 >= lo && r2 <= hi)) continue;
            if (r1 > one) r1 = one;
            if (r1 < zero) r1 = zero;
            if (r2 > one) r2 = one;
            if (r2 < zero) r2 = zero;
            out.push_back(amap2_assemble(h1, h2, r1, r2, 2));
        }
    }
    return out;
}

/** amap2_fitall_gamma with the MATLAB defaults degentol = 1e-8, r12tol = 1e-6. */
template <class T>
std::vector<Map<T>> amap2_fitall_gamma(const T& M1, const T& M2, const T& M3, const T& GAMMA) {
    return amap2_fitall_gamma(M1, M2, M3, GAMMA, T(num_traits<T>::from_double(1e-8)),
                              T(num_traits<T>::from_double(1e-6)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_AMAP2_FITALL_GAMMA_H
