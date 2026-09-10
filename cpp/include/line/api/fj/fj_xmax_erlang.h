/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_ERLANG_H
#define LINE_API_FJ_XMAX_ERLANG_H

/**
 * Expected maximum of K i.i.d. Erlang-k service times.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_erlang.m, cross-checked against
 * FJ_xmax.fj_xmax_erlang in jar/src/main/java/jline/api/fj/FJ_xmax.java
 * (identical formulas; the JAR replaces MATLAB's adaptive `integral` with a
 * 10001-point composite Simpson rule, which this port also does).
 *
 * MIXED ARITHMETIC. At k = 2 the MATLAB file uses the closed form
 *
 *   X_K^max = (1/mu) sum_{n=1..K} C(K,n) (-1)^{n-1} sum_{m=1..n} C(n,m) m! / (2 n^{m+1})
 *
 * which is rational and therefore exact in any field, cancellation included:
 * the alternating outer sum is exactly the kind that double cannot carry past
 * about K = 25. For any other k the mean is a quadrature of 1 - F(t)^K against
 * the Erlang CDF and needs exp, so it is only available when T carries
 * transcendental functions; asking for it at exact arithmetic throws
 * UnsupportedError rather than silently substituting something else.
 */

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

namespace detail {

/** Erlang-k CDF, 1 - exp(-mu t) sum_{j<k} (mu t)^j / j!. */
template <class T>
T erlang_cdf(const T& t, unsigned k, const T& mu) {
    static_assert(num_traits<T>::has_transcendental,
                  "erlang_cdf requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    if (t <= zero) return zero;
    const T x = mu * t;
    T S = zero;
    for (unsigned j = 0; j < k; ++j) S += num_pow_int(x, j) / num_factorial<T>(j);
    return num_traits<T>::from_int(1) - num_exp(T(-x)) * S;
}

/** Erlang-k survival function, exp(-mu x) sum_{j<k} (mu x)^j / j!. */
template <class T>
T erlang_survival(const T& x, unsigned k, const T& mu) {
    static_assert(num_traits<T>::has_transcendental,
                  "erlang_survival requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0);
    if (x <= zero) return num_traits<T>::from_int(1);
    const T y = mu * x;
    T S = zero;
    for (unsigned j = 0; j < k; ++j) S += num_pow_int(y, j) / num_factorial<T>(j);
    return num_exp(T(-y)) * S;
}

/** Quadrature branch of fj_xmax_erlang, for k != 2. */
template <class T>
T fj_xmax_erlang_quad(unsigned K, unsigned k, const T& mu) {
    const T kk = num_traits<T>::from_int(static_cast<long>(k));
    const T upper = kk / mu * num_traits<T>::from_int(10) + num_traits<T>::from_int(10) * num_sqrt(kk) / mu;
    const T one = num_traits<T>::from_int(1);
    return simpson<T>([&](const T& t) { return T(one - num_pow_int(erlang_cdf(t, k, mu), K)); },
                      num_traits<T>::from_int(0), upper);
}

}  // namespace detail

/**
 * @param K  number of branches, K >= 1
 * @param k  Erlang stages, k >= 1
 * @param mu per-stage rate, > 0 (branch mean is k/mu)
 * @return   expected maximum of K Erlang-k samples
 */
template <class T>
T fj_xmax_erlang(unsigned K, unsigned k, const T& mu) {
    detail::require_positive_K(K, "fj_xmax_erlang");
    if (k < 1) throw InputError("fj_xmax_erlang: the stage count k must be a positive integer");
    if (mu <= num_traits<T>::from_int(0)) throw InputError("fj_xmax_erlang: the rate mu must be positive");

    if (k == 2) {
        T outer = num_traits<T>::from_int(0);
        const T two = num_traits<T>::from_int(2);
        for (unsigned n = 1; n <= K; ++n) {
            const T nn = num_traits<T>::from_int(static_cast<long>(n));
            T inner = num_traits<T>::from_int(0);
            for (unsigned m = 1; m <= n; ++m)
                inner += detail::fj_binom<T>(n, m) * num_factorial<T>(m) / (two * num_pow_int(nn, m + 1));
            const T term = detail::fj_binom<T>(K, n) * inner;
            if ((n - 1) % 2 == 0) outer += term;
            else outer -= term;
        }
        return outer / mu;
    }

    if constexpr (num_traits<T>::has_transcendental) {
        return detail::fj_xmax_erlang_quad<T>(K, k, mu);
    } else {
        throw UnsupportedError(
            "fj_xmax_erlang: only k = 2 has a closed form; any other stage count needs a "
            "quadrature of the Erlang CDF and therefore transcendental arithmetic");
    }
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_ERLANG_H
