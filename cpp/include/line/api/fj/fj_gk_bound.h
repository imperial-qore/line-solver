/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_GK_BOUND_H
#define LINE_API_FJ_GK_BOUND_H

/**
 * G(K) factors for the standardized-maximum approximation X_K^max ~ mu + sigma G(K).
 *
 * Templated port of matlab/src/api/fj/fj_gk_bound.m. The JAR carries the same
 * four values in jline.api.fj.GKBoundResult, computed inside
 * FJ_xmax.fj_xmax_approx (identical formulas).
 *
 *   exponential: G(K) = H_K - 1
 *   uniform:     G(K) = sqrt(3) (K-1)/(K+1)
 *   evd:         G(K) = sqrt(6) ln(K) / pi
 *   upper bound: G(K) = (K-1)/sqrt(2K-1)                     (David 1970)
 *
 * static_assert(num_traits<T>::has_transcendental) -- three of the four
 * involve sqrt or log, so the struct as a whole is only defined for the
 * inexact number types. The exponential entry alone is rational and is
 * reachable exactly through fj_harmonic.
 */

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K number of branches, K >= 1
 * @return  all four G(K) factors, matching MATLAB's 'all' mode
 */
template <class T>
FJGKBoundResult<T> fj_gk_bound(unsigned K) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_gk_bound requires transcendental arithmetic");
    detail::require_positive_K(K, "fj_gk_bound");
    const T one = num_traits<T>::from_int(1);
    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    FJGKBoundResult<T> r;
    r.K = K;
    r.exponential = fj_harmonic<T>(K) - one;
    r.uniform = detail::num_sqrt(T(num_traits<T>::from_int(3))) * (Kt - one) / (Kt + one);
    r.evd = detail::num_sqrt(T(num_traits<T>::from_int(6))) * detail::num_log(Kt) / detail::num_pi<T>();
    r.upper_bound = (Kt - one) / detail::num_sqrt(T(num_traits<T>::from_int(2) * Kt - one));
    return r;
}

/** Single-family accessor, matching MATLAB's 'exp'/'uniform'/'evd'/'bound' modes. */
template <class T>
T fj_gk_bound(unsigned K, FJDistType type) {
    const FJGKBoundResult<T> all = fj_gk_bound<T>(K);
    switch (type) {
        case FJDistType::Exp: return all.exponential;
        case FJDistType::Uniform: return all.uniform;
        case FJDistType::Evd: return all.evd;
        case FJDistType::Bound: return all.upper_bound;
    }
    throw InputError("fj_gk_bound: unknown distribution type");
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_GK_BOUND_H
