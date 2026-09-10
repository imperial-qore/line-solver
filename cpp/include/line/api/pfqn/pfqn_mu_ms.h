/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MU_MS_H
#define LINE_API_PFQN_MU_MS_H

/**
 * Aggregate load-dependent rate of m identical c-server FCFS stations.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mu_ms.m.
 *
 * A single c-server station holding k jobs has balance function
 * beta_1(k) = 1 / prod_{j=1}^{k} min(j, c). The flow-equivalent aggregate of m
 * such stations in parallel has the convolution
 *
 *   beta_m(n) = sum_{k=0}^{n} beta_1(k) beta_{m-1}(n - k),   beta_m(0) = 1,
 *
 * and the load-dependent rate of the aggregate is the ratio of consecutive
 * balance-function values,
 *
 *   mu(n) = beta_m(n - 1) / beta_m(n),   n = 1, ..., N.
 *
 * The reference writes the inner term as 1/(prod(a) * prod(b)) with
 * b = 1/beta_{m-1}(n-k), i.e. as beta_{m-1}(n-k)/prod_{j<=k} min(j,c), which is
 * the convolution above; the MATLAB idiom relies on prod([]) = 1 to cover
 * k = 0. It also fills the table with the population outermost and the station
 * count innermost, so beta_{m-1}(n) is available at the same n; that ordering
 * is preserved here even though the port could iterate either way.
 *
 * Arithmetic: EXACT-CAPABLE. Only additions, multiplications and divisions in
 * the field of the inputs, all of them on values built from the integers, so
 * at T = Rational the returned rates are exact rationals.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/**
 * @param N maximum population
 * @param m number of identical stations
 * @param c number of servers per station
 * @return (N) rates mu(1), ..., mu(N)
 */
template <class T>
std::vector<T> pfqn_mu_ms(int N, int m, int c) {
    if (N < 0) throw InputError("pfqn_mu_ms: negative population");
    if (m < 1) throw InputError("pfqn_mu_ms: the station count must be at least one");
    if (c < 1) throw InputError("pfqn_mu_ms: the server count must be at least one");

    const T one = num_traits<T>::from_int(1);
    const T zero = num_traits<T>::from_int(0);

    // beta_1(k) = 1 / prod_{j=1}^{k} min(j, c), tabulated once.
    std::vector<T> beta1(static_cast<std::size_t>(N) + 1, one);
    for (int k = 1; k <= N; ++k) {
        const int mn = k < c ? k : c;
        beta1[static_cast<std::size_t>(k)] =
            beta1[static_cast<std::size_t>(k - 1)] / num_traits<T>::from_int(mn);
    }

    // g[i][n] = beta_{i+1}(n); population outermost, as in the reference.
    std::vector<std::vector<T>> g(static_cast<std::size_t>(m),
                                  std::vector<T>(static_cast<std::size_t>(N) + 1, zero));
    for (int n = 0; n <= N; ++n) {
        for (int i = 1; i <= m; ++i) {
            if (n == 0) {
                g[static_cast<std::size_t>(i - 1)][0] = one;
            } else if (i == 1) {
                g[0][static_cast<std::size_t>(n)] = beta1[static_cast<std::size_t>(n)];
            } else {
                T s = zero;
                for (int k = 0; k <= n; ++k)
                    s += beta1[static_cast<std::size_t>(k)] *
                         g[static_cast<std::size_t>(i - 2)][static_cast<std::size_t>(n - k)];
                g[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(n)] = s;
            }
        }
    }

    std::vector<T> mu(static_cast<std::size_t>(N), one);
    for (int n = 1; n <= N; ++n) {
        const T& den = g[static_cast<std::size_t>(m - 1)][static_cast<std::size_t>(n)];
        if (den == zero) throw NumericError("pfqn_mu_ms: zero balance function");
        mu[static_cast<std::size_t>(n - 1)] =
            g[static_cast<std::size_t>(m - 1)][static_cast<std::size_t>(n - 1)] / den;
    }
    return mu;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MU_MS_H
