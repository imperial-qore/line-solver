/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_CHAR_MAX_H
#define LINE_API_FJ_CHAR_MAX_H

/**
 * Gravey's characteristic maximum M_K, an upper bound on the expected maximum
 * of K i.i.d. random variables.
 *
 * Templated port of matlab/src/api/fj/fj_char_max.m, cross-checked against
 * FJ_char_max.fj_char_max_exp and fj_char_max_erlang in
 * jar/src/main/java/jline/api/fj/FJ_char_max.java (identical; the JAR splits
 * the MATLAB dist_type switch into two entry points, which this port mirrors
 * as two overloads).
 *
 * m_K is the greatest lower bound with P(X > m_K) <= 1/K, and
 * M_K = m_K + K int_{m_K}^inf P(X > x) dx.
 *
 *   exponential: m_K = ln(K)/mu,  M_K = H_K/mu   (the bound is tight here)
 *   Erlang-k:    m_K solves exp(-mu m) sum_{i<k} (mu m)^i/i! = 1/K,
 *                M_K = (k/mu)[1 + K exp(-mu m_K) (mu m_K)^k / k!]
 *
 * static_assert(num_traits<T>::has_transcendental) -- a log in the
 * exponential case and a bracketed root of a transcendental equation in the
 * Erlang one.
 *
 * The MATLAB 'general' mode, which takes an arbitrary survival function
 * handle and integrates its tail to infinity by adaptive quadrature after a
 * doubling search for the truncation point, is NOT ported: the truncation
 * search has no termination guarantee for a heavy tail and the JAR omits it
 * too. Callers with a specific tail should use fj_order_stat.
 */

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/api/fj/fj_xmax_erlang.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * Exponential branch.
 *
 * @param K  number of samples, K >= 1
 * @param mu rate, > 0
 * @return   [MK, mK]
 */
template <class T>
FJCharMaxResult<T> fj_char_max(unsigned K, const T& mu) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_char_max requires transcendental arithmetic");
    detail::require_positive_K(K, "fj_char_max");
    if (mu <= num_traits<T>::from_int(0)) throw InputError("fj_char_max: the rate mu must be positive");
    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    return {T(fj_harmonic<T>(K) / mu), T(detail::num_log(Kt) / mu)};
}

/**
 * Erlang-k branch.
 *
 * @param K  number of samples, K >= 1
 * @param k  Erlang stages, k >= 1
 * @param mu per-stage rate, > 0
 * @return   [MK, mK]
 */
template <class T>
FJCharMaxResult<T> fj_char_max(unsigned K, unsigned k, const T& mu) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_char_max requires transcendental arithmetic");
    detail::require_positive_K(K, "fj_char_max");
    if (k < 1) throw InputError("fj_char_max: the Erlang stage count k must be positive");
    if (mu <= num_traits<T>::from_int(0)) throw InputError("fj_char_max: the rate mu must be positive");

    const T one = num_traits<T>::from_int(1);
    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    const T kt = num_traits<T>::from_int(static_cast<long>(k));
    const T target = one / Kt;

    T hi = kt / mu + detail::num_log(Kt) / mu;
    if (hi <= num_traits<T>::from_int(0)) hi = kt / mu;
    while (detail::erlang_survival(hi, k, mu) > target) hi *= num_traits<T>::from_int(2);
    const T mK = detail::bisect<T>(
        [&](const T& x) { return T(detail::erlang_survival(x, k, mu) - target); },
        num_traits<T>::from_int(0), hi, "fj_char_max");

    const T MK = (kt / mu) * (one + Kt * detail::num_exp(T(-mu * mK)) * num_pow_int(T(mu * mK), k) /
                                        num_factorial<T>(k));
    return {MK, mK};
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_CHAR_MAX_H
