/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_EMMA_H
#define LINE_API_FJ_XMAX_EMMA_H

/**
 * EMMA (Extreme-value Maximum Moment Approximation) to the expected maximum
 * of K i.i.d. samples.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_emma.m, cross-checked against
 * FJ_xmax.fj_xmax_emma in jar/src/main/java/jline/api/fj/FJ_xmax.java
 * (identical; both use the same rounded constant phi = 0.570376, which is
 * exp(-exp(-gamma)) to six places).
 *
 *   exponential: E[Y_K] = -(1/mu) ln(1 - phi^{1/K})
 *   general:     E[Y_K] = F^{-1}(phi^{1/K})
 *
 * static_assert(num_traits<T>::has_transcendental) -- a real root phi^{1/K}
 * and a log. The general form takes the quantile function as a callable, so
 * the caller supplies whatever inverse CDF it has.
 */

#include <functional>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * Exponential branch.
 *
 * @param K  number of samples, K >= 1
 * @param mu rate of the exponential, > 0
 */
template <class T>
T fj_xmax_emma(unsigned K, const T& mu) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_xmax_emma requires transcendental arithmetic");
    detail::require_positive_K(K, "fj_xmax_emma");
    if (mu <= num_traits<T>::from_int(0)) throw InputError("fj_xmax_emma: the rate mu must be positive");
    const T one = num_traits<T>::from_int(1);
    const T phi = num_traits<T>::from_double(0.570376);
    const T root = detail::num_pow(phi, T(one / num_traits<T>::from_int(static_cast<long>(K))));
    return -(one / mu) * detail::num_log(T(one - root));
}

/**
 * General branch: the caller supplies the quantile function F^{-1}.
 *
 * @param K    number of samples, K >= 1
 * @param Finv inverse CDF of the branch distribution
 */
template <class T>
T fj_xmax_emma(unsigned K, const std::function<T(const T&)>& Finv) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_xmax_emma requires transcendental arithmetic");
    detail::require_positive_K(K, "fj_xmax_emma");
    if (!Finv) throw InputError("fj_xmax_emma: the inverse CDF must be callable");
    const T one = num_traits<T>::from_int(1);
    const T phi = num_traits<T>::from_double(0.570376);
    return Finv(T(detail::num_pow(phi, T(one / num_traits<T>::from_int(static_cast<long>(K))))));
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_EMMA_H
