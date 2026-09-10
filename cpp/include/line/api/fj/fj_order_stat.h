/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_ORDER_STAT_H
#define LINE_API_FJ_ORDER_STAT_H

/**
 * CDF and expected value of the k-th order statistic of K i.i.d. samples.
 *
 * Templated port of matlab/src/api/fj/fj_order_stat.m. The JAR carries the
 * same CDF in jline.api.fj.FJ_order_stat (identical).
 *
 *   F_{Y_k}(y) = sum_{j=k..K} C(K,j) F(y)^j (1 - F(y))^{K-j}
 *   F_{Y_K}(y) = F(y)^K                                       (the maximum)
 *
 * MIXED ARITHMETIC. The CDF is a polynomial in the value of the base CDF, so
 * it is exact in any field once F(y) is known -- and it satisfies the exact
 * identity sum_{k=1..K} F_{Y_k} = K F, plus F_{Y_1} = 1 - (1-F)^K, both of
 * which are bit-exact only in rational arithmetic.
 *
 * The expected value is a quadrature. For k = 1 and k = K it integrates the
 * survival function directly; for an interior k the MATLAB file differentiates
 * the supplied CDF by a central difference with a fixed step 1e-8, which is a
 * tolerance-driven approximation. Both are therefore produced only when T
 * carries transcendental functions, and mean_available reports which branch
 * ran. The truncation point follows MATLAB: double from 1 until F(u) >= 1-1e-6,
 * at most 100 times.
 */

#include <functional>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param y   evaluation point of the CDF
 * @param k   order of the statistic, 1 = minimum, K = maximum
 * @param K   number of samples
 * @param F_X CDF of the branch distribution
 * @return    [F_Yk, E_Yk, mean_available]
 */
template <class T>
FJOrderStatResult<T> fj_order_stat(const T& y, unsigned k, unsigned K,
                                   const std::function<T(const T&)>& F_X) {
    if (k < 1 || k > K) throw InputError("fj_order_stat: k must satisfy 1 <= k <= K");
    if (!F_X) throw InputError("fj_order_stat: the CDF must be callable");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    const T F_y = F_X(y);
    T F_Yk = zero;
    if (k == K) {
        F_Yk = num_pow_int(F_y, K);
    } else {
        for (unsigned j = k; j <= K; ++j)
            F_Yk += detail::fj_binom<T>(K, j) * num_pow_int(F_y, j) * num_pow_int(T(one - F_y), K - j);
    }

    if constexpr (num_traits<T>::has_transcendental) {
        // Truncation point: double from 1 until the CDF is essentially 1.
        const T target = one - num_traits<T>::from_double(1e-6);
        T upper = one;
        for (unsigned it = 0; it < 100 && F_X(upper) < target; ++it) upper *= num_traits<T>::from_int(2);

        T E_Yk = zero;
        if (k == K) {
            E_Yk = detail::simpson<T>([&](const T& t) { return T(one - num_pow_int(F_X(t), K)); }, zero, upper);
        } else if (k == 1) {
            E_Yk = detail::simpson<T>([&](const T& t) { return num_pow_int(T(one - F_X(t)), K); }, zero, upper);
        } else {
            const T eps = num_traits<T>::from_double(1e-8);
            const T coeff = num_traits<T>::from_int(static_cast<long>(K)) * detail::fj_binom<T>(K - 1, k - 1);
            E_Yk = detail::simpson<T>(
                [&](const T& t) {
                    const T f = (F_X(T(t + eps)) - F_X(T(t - eps))) / (num_traits<T>::from_int(2) * eps);
                    const T Ft = F_X(t);
                    return T(t * coeff * f * num_pow_int(Ft, k - 1) * num_pow_int(T(one - Ft), K - k));
                },
                zero, upper);
        }
        return {F_Yk, E_Yk, true};
    } else {
        return {F_Yk, zero, false};
    }
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_ORDER_STAT_H
