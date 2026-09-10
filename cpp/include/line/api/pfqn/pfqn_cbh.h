/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_CBH_H
#define LINE_API_PFQN_PFQN_CBH_H

/**
 * Convolutional Bound Hierarchy (Dowdy, Eager, Gordon and Saxton 1984) on the
 * throughput of a single-class closed product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_cbh.m. Column c = M - level of
 * Buzen's g array is filled from a Balanced Job Bound estimate of the first c
 * stations (e_0 = 1, e_i = e_{i-1}/B_i), the remaining `level` stations are
 * convolved exactly, and the infinite-server delay is convolved exactly as
 * well. The BJB upper fill yields the upper bound, the lower fill the lower
 * bound, and both meet the exact solution at level = M.
 *
 * ARITHMETIC. Every step is an addition, a multiplication or a division of
 * field elements -- the balanced fill uses the arithmetic mean and the maximum
 * of the demands, not a root or a logarithm -- so the bound is EXACT in
 * rational arithmetic and is deliberately left ungated. The only non-rational
 * ingredient, Z^j/j!, is a rational for rational Z.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_cbh, mirroring [Xlo, Xhi]. */
template <class T>
struct CbhBounds {
    T Xlo;
    T Xhi;
};

namespace detail {

/** One side of the hierarchy: BJB-filled column c, then exact convolution. */
template <class T>
T cbh_hier(const std::vector<T>& L, int N, const T& Z, std::size_t c, bool upper) {
    const std::size_t M = L.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    T Rc = zero, Lbc = L[0];
    for (std::size_t i = 0; i < c; ++i) {
        Rc += L[i];
        if (L[i] > Lbc) Lbc = L[i];
    }
    const T Lac = T(Rc / num_traits<T>::from_int(static_cast<long>(c)));

    std::vector<T> g(static_cast<std::size_t>(N) + 1, zero);
    if (c == 1) {
        // The single-server column is already exact: g(n) = L_1^n.
        for (int n = 0; n <= N; ++n) g[static_cast<std::size_t>(n)] = num_pow_int(L[0], static_cast<unsigned>(n));
    } else {
        g[0] = one;
        for (int i = 1; i <= N; ++i) {
            const T iT = num_traits<T>::from_int(i);
            const T den = upper ? T(Rc + T(iT - one) * Lac) : T(Rc + T(iT - one) * Lbc);
            const T B = T(iT / den);  // BJB fill; upper fill -> upper bound
            if (B == zero) throw NumericError("pfqn_cbh: zero balanced-job-bound fill");
            g[static_cast<std::size_t>(i)] = T(g[static_cast<std::size_t>(i) - 1] / B);
        }
    }
    for (std::size_t m = c; m < M; ++m)
        for (int n = 1; n <= N; ++n)
            g[static_cast<std::size_t>(n)] += L[m] * g[static_cast<std::size_t>(n) - 1];

    if (Z > zero) {
        std::vector<T> gd(static_cast<std::size_t>(N) + 1, zero);
        // Poisson weight Z^j/j!. In double the naive ratio overflows for j >~ 171
        // and j runs to the POPULATION; an exact T cannot overflow and keeps it.
        for (int j = 0; j <= N; ++j) {
            if constexpr (num_traits<T>::has_transcendental) {
                using std::exp;
                using std::log;
                const T jT = num_traits<T>::from_int(j);
                gd[static_cast<std::size_t>(j)] = T(exp(jT * log(Z) - detail::num_factln<T>(jT)));
            } else {
                gd[static_cast<std::size_t>(j)] =
                    T(num_pow_int(Z, static_cast<unsigned>(j)) / num_factorial<T>(static_cast<unsigned>(j)));
            }
        }
        std::vector<T> gfull(static_cast<std::size_t>(N) + 1, zero);
        for (int n = 0; n <= N; ++n) {
            T acc = zero;
            for (int j = 0; j <= n; ++j)
                acc += g[static_cast<std::size_t>(j)] * gd[static_cast<std::size_t>(n - j)];
            gfull[static_cast<std::size_t>(n)] = acc;
        }
        g = gfull;
    }
    if (g[static_cast<std::size_t>(N)] == zero) throw NumericError("pfqn_cbh: zero normalizing constant");
    return T(g[static_cast<std::size_t>(N) - 1] / g[static_cast<std::size_t>(N)]);
}

}  // namespace detail

/**
 * @param L     (M) per-station demands
 * @param N     population
 * @param Z     think time
 * @param level number of exactly convolved stations, clamped to [1, M]
 */
template <class T>
CbhBounds<T> pfqn_cbh(const std::vector<T>& L, int N, const T& Z, int level) {
    const std::size_t M = L.size();
    if (M == 0) throw InputError("pfqn_cbh: empty demand vector");
    if (N < 1) throw InputError("pfqn_cbh: population must be at least one");
    int lv = std::max(1, std::min(level, static_cast<int>(M)));
    const std::size_t c = static_cast<std::size_t>(std::max(1, static_cast<int>(M) - lv));
    CbhBounds<T> r;
    r.Xlo = detail::cbh_hier(L, N, Z, c, false);
    r.Xhi = detail::cbh_hier(L, N, Z, c, true);
    return r;
}

template <class T>
CbhBounds<T> pfqn_cbh(const std::vector<T>& L, int N, const T& Z) {
    return pfqn_cbh(L, N, Z, 2);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_CBH_H
