/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_PARETO_H
#define LINE_API_FJ_XMAX_PARETO_H

/**
 * Expected maximum and characteristic maximum of K i.i.d. shifted-Pareto
 * samples with survival S(x) = (k/(k+x))^beta.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_pareto.m, cross-checked against
 * FJ_xmax.fj_xmax_pareto and fj_xmax_pareto_char_max in
 * jar/src/main/java/jline/api/fj/FJ_xmax.java (identical formulas; the JAR
 * replaces MATLAB's adaptive `integral` with a 10001-point composite Simpson
 * rule on the same truncated range, which this port also does).
 *
 *   Xmax = int_0^inf [1 - F(x)^K] dx,  truncated at k K^{2/beta} * 10
 *   m_K  = k (K^{1/beta} - 1)
 *   M_K  = m_K + K k^beta (k + m_K)^{1-beta} / (beta - 1)
 *
 * static_assert(num_traits<T>::has_transcendental) -- real powers throughout
 * plus the quadrature. Worth flagging: the truncation point is a heuristic and
 * the Pareto tail is heavy, so Xmax is systematically underestimated; the
 * shortfall grows as beta approaches 2 and the function rejects beta <= 2
 * outright because the mean of the maximum is then the only finite moment
 * left. M_K, by contrast, is a closed form and is exact.
 */

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K    number of samples, K >= 1
 * @param beta Pareto shape, must exceed 2 for finite moments
 * @param k    Pareto scale, > 0 (MATLAB defaults it to beta - 1, the value
 *             that makes the branch mean equal to 1)
 * @return     [Xmax, MK]
 */
template <class T>
FJXmaxParetoResult<T> fj_xmax_pareto(unsigned K, const T& beta, const T& k) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_xmax_pareto requires transcendental arithmetic");
    detail::require_positive_K(K, "fj_xmax_pareto");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (beta <= num_traits<T>::from_int(2))
        throw InputError("fj_xmax_pareto: the shape beta must exceed 2 for finite moments");
    if (k <= zero) throw InputError("fj_xmax_pareto: the scale k must be positive");

    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    const T upper = k * detail::num_pow(Kt, T(num_traits<T>::from_int(2) / beta)) * num_traits<T>::from_int(10);
    const T Xmax = detail::simpson<T>(
        [&](const T& x) {
            const T xx = x > zero ? x : zero;
            const T F = one - detail::num_pow(T(k / (k + xx)), beta);
            return T(one - num_pow_int(F, K));
        },
        zero, upper);

    const T mK = k * (detail::num_pow(Kt, T(one / beta)) - one);
    const T tail = detail::num_pow(k, beta) * detail::num_pow(T(k + mK), T(one - beta)) / (beta - one);
    return {Xmax, T(mK + Kt * tail)};
}

/** MATLAB's default scale k = beta - 1, which normalizes the branch mean to 1. */
template <class T>
FJXmaxParetoResult<T> fj_xmax_pareto(unsigned K, const T& beta) {
    return fj_xmax_pareto(K, beta, T(beta - num_traits<T>::from_int(1)));
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_PARETO_H
