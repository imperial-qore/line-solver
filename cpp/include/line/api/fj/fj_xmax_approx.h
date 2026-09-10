/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_APPROX_H
#define LINE_API_FJ_XMAX_APPROX_H

/**
 * Two-moment approximation to the expected maximum of K i.i.d. samples,
 * X_K^max ~ mu_X + sigma_X G(K).
 *
 * Templated port of matlab/src/api/fj/fj_xmax_approx.m, cross-checked against
 * FJ_xmax.fj_xmax_approx in jar/src/main/java/jline/api/fj/FJ_xmax.java
 * (identical).
 *
 * MIXED ARITHMETIC. The exponential family uses G(K) = H_K - 1, which is
 * rational and exact in any field; the uniform, EVD and bound families need
 * sqrt or log and are only available when T carries transcendental functions.
 * Asking for one of those at exact arithmetic throws UnsupportedError.
 */

#include "line/api/fj/fj_gk_bound.h"
#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K       number of branches, K >= 1
 * @param mu_X    mean of the branch distribution
 * @param sigma_X standard deviation of the branch distribution, >= 0
 * @param type    which G(K) family to use
 * @return        [Xmax, GK]
 */
template <class T>
FJXmaxApproxResult<T> fj_xmax_approx(unsigned K, const T& mu_X, const T& sigma_X,
                                     FJDistType type = FJDistType::Exp) {
    detail::require_positive_K(K, "fj_xmax_approx");
    if (sigma_X < num_traits<T>::from_int(0))
        throw InputError("fj_xmax_approx: sigma_X must be non-negative");

    T GK = num_traits<T>::from_int(0);
    if (type == FJDistType::Exp) {
        GK = fj_harmonic<T>(K) - num_traits<T>::from_int(1);
    } else {
        if constexpr (num_traits<T>::has_transcendental) {
            GK = fj_gk_bound<T>(K, type);
        } else {
            throw UnsupportedError(
                "fj_xmax_approx: only the exponential G(K) is rational; the uniform, EVD and "
                "bound families need sqrt or log and therefore transcendental arithmetic");
        }
    }
    return {T(mu_X + sigma_X * GK), GK};
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_APPROX_H
