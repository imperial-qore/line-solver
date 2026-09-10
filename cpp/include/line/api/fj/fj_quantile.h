/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_QUANTILE_H
#define LINE_API_FJ_QUANTILE_H

/**
 * Quantile of the maximum of K i.i.d. samples.
 *
 * Templated port of matlab/src/api/fj/fj_quantile.m, cross-checked against
 * jar/src/main/java/jline/api/fj/FJ_quantile.java (identical).
 *
 *   Gumbel approximation: x(K,q) = ln K - ln ln(1/q)
 *   exact, given F^-1:    x(K,q) = F^{-1}(q^{1/K})
 *
 * static_assert(num_traits<T>::has_transcendental) -- a log in the first form
 * and a real root in the second. The Gumbel form is documented as inaccurate
 * at small K and is not a bound in either direction, so it should not be used
 * to certify a service-level target.
 */

#include <functional>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * Gumbel approximation.
 *
 * @param K number of samples, K >= 1
 * @param q quantile level, 0 < q < 1
 */
template <class T>
T fj_quantile(unsigned K, const T& q) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_quantile requires transcendental arithmetic");
    detail::require_positive_K(K, "fj_quantile");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (q <= zero || q >= one) throw InputError("fj_quantile: q must satisfy 0 < q < 1");
    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    return detail::num_log(Kt) - detail::num_log(T(detail::num_log(T(one / q))));
}

/**
 * Exact quantile through the supplied inverse CDF.
 *
 * @param K    number of samples, K >= 1
 * @param q    quantile level, 0 < q < 1
 * @param Finv inverse CDF of the branch distribution
 */
template <class T>
T fj_quantile(unsigned K, const T& q, const std::function<T(const T&)>& Finv) {
    static_assert(num_traits<T>::has_transcendental,
                  "fj_quantile requires transcendental arithmetic");
    detail::require_positive_K(K, "fj_quantile");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (q <= zero || q >= one) throw InputError("fj_quantile: q must satisfy 0 < q < 1");
    if (!Finv) throw InputError("fj_quantile: the inverse CDF must be callable");
    return Finv(T(detail::num_pow(q, T(one / num_traits<T>::from_int(static_cast<long>(K))))));
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_QUANTILE_H
