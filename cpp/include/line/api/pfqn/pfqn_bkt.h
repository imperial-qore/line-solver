/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_BKT_H
#define LINE_API_PFQN_PFQN_BKT_H

/**
 * Knessl-Tier expansion with the Stirling-remainder correction (BKT).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_bkt.m. pfqn_kt extracts N from the
 * generating function of G by steepest descent on
 *
 *   F(xi) = sum_r Z_r xi_r - sum_k log(1 - U_k) - sum_r N_r log xi_r,  U = L xi.
 *
 * On the demand-free integral the exact coefficient is [xi^N] exp(Z xi) = Z^N/N!,
 * whereas the expansion returns N log Z - (N log N - N + log(2 pi N)/2), which is
 * Stirling's approximation of log(N!) in place of log(N!). So KT lies ABOVE the
 * exact value by the remainder
 *
 *   s(N) = log(N!) - (N log N - N + log(2 pi N)/2)
 *        = lgamma(N+1) - (N + 1/2) log N + N - log(2 pi)/2
 *
 * per Laplaced class direction, and BKT subtracts sum_r s(N_r). s(1) is the
 * constant pfqn_ble adds per STATION direction, and s(N) = 1/(12 N) + O(N^-2);
 * the remainder is evaluated exactly, since truncating it at 1/(12 N) loses an
 * order of magnitude (median |error| on the 1562 models of Cas17 sec5.3.1: 0.083
 * nats for KT, 1.9e-4 for the truncation, 1.4e-5 for the exact remainder). With
 * a think time BKT is the SAME estimator as pfqn_ble, to the accuracy of the
 * two saddle-point solvers; without one they differ by the constant
 * kappa - r(N+M). See _kb/03-api-layer.md.
 *
 * WHICH CLASSES. Only the classes pfqn_kt actually Laplaces are corrected: a class
 * with no jobs is dropped by its recursion, and a self-looping class (one nonzero
 * demand and no think time, folded only when more than one class remains) has its
 * coefficient extracted exactly, so neither carries a remainder. The predicate
 * below is pfqn_kt's own, re-derived here rather than trusting R.
 *
 * ARITHMETIC. Inherited from pfqn_kt: gated on num_traits<T>::has_transcendental.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_kt.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_bkt, mirroring [Gn, lGn] (X and Q are pfqn_kt's seeds). */
template <class T>
using BktResult = KtResult<T>;

/** s(N) = log(N!) - (N log N - N + log(2 pi N)/2), exactly, for N >= 1. */
template <class T>
T pfqn_stirling_remainder(const T& n) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_stirling_remainder requires transcendental arithmetic");
    using std::log;
    const T half = num_traits<T>::from_double(0.5);
    const T twopi = num_traits<T>::from_double(6.283185307179586476925286766559);
    return T(detail::num_factln<T>(n) - (n + half) * log(n) + n - log(twopi) / num_traits<T>::from_int(2));
}

/**
 * @param L (M x R) demands, @param N (R) population, @param Z (R) think times
 *          (pass an empty vector or all zeros for the Z = 0 branch)
 */
template <class T>
BktResult<T> pfqn_bkt(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_bkt requires transcendental arithmetic (steepest-descent expansion of log G)");
    using std::exp;
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = L.rows(), R = L.cols();
    std::vector<T> Zc = Z;
    if (Zc.empty()) Zc.assign(R, zero);

    BktResult<T> res = pfqn_kt(L, N, Zc);

    std::size_t nkeep = 0;
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] > zero) ++nkeep;
    T corr = zero;
    for (std::size_t r = 0; r < R; ++r) {
        if (!(N[r] > zero)) continue;  // dropped by pfqn_kt's recursion
        if (nkeep > 1) {               // pfqn_kt folds self-loops only with Rorig > 1
            std::size_t nnz = 0;
            for (std::size_t i = 0; i < M; ++i)
                if (L(i, r) != zero) ++nnz;
            if (nnz == 1 && Zc[r] == zero) continue;  // extracted exactly, no remainder
        }
        corr += pfqn_stirling_remainder<T>(N[r]);
    }
    res.lG -= corr;
    res.G = exp(res.lG);
    return res;
}

template <class T>
BktResult<T> pfqn_bkt(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_bkt(L, N, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_BKT_H
