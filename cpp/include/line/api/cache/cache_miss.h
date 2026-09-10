/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_MISS_H
#define LINE_API_CACHE_MISS_H

/**
 * Exact cache miss rates from the recursive normalizing constant.
 *
 * Templated port of matlab/src/api/cache/cache_miss.m, cross-checked against
 * jar/src/main/java/jline/api/cache/Cache_miss.java.
 *
 * Two quantities are computed, both as ratios of cache_erec constants:
 *   M       = E(gamma, m + e_1) / E(gamma, m), the global miss rate; and
 *   pi0(k)  = E(gamma without item k, m) / E(gamma, m), the probability that
 *             item k is absent from the cache,
 * from which the per-user rate MU(v) = sum_k lambda(v,k) pi0(k) and the
 * per-item rate MI(k) = (sum_v lambda(v,k)) pi0(k) follow.
 *
 * Pure field operations throughout, so the exact instantiation returns rates
 * with no rounding.
 *
 * DIVERGENCE, MATLAB vs JAR: MATLAB conditions on the absence of item k by
 * deleting ROW k of gamma, gamma(setdiff(1:n,k),:) -- gamma is item-by-list,
 * so a row is an item. The JAR's Cache_miss deletes COLUMN k instead (it
 * builds gammaWithoutK by iterating over gamma.getNumCols() and skipping
 * j == k), which removes a cache LIST, not an item, and additionally reads the
 * item count as lambda.getNumCols() while indexing gamma by it. The JAR's
 * pi0/MU/MI are therefore wrong whenever h != n, and meaningless in general.
 * This port follows MATLAB.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_erec.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Return value of cache_miss, mirroring [M,MU,MI,pi0]. */
template <class T>
struct CacheMissResult {
    T M;                  ///< global miss rate
    std::vector<T> MU;    ///< (u) per-user miss rate; empty when no lambda given
    std::vector<T> MI;    ///< (n) per-item miss rate; empty when no lambda given
    std::vector<T> pi0;   ///< (n) per-item miss probability; empty when no lambda given
};

/**
 * @param gamma  (n x h) access factors
 * @param m      (h) list capacities
 * @param lambda (u x n) per-user per-item request rates, MATLAB's
 *               lambda(:,:,1); pass an empty matrix for the miss rate alone
 */
template <class T>
CacheMissResult<T> cache_miss(const Matrix<T>& gamma, const std::vector<int>& m,
                              const Matrix<T>& lambda) {
    if (gamma.cols() != m.size())
        throw InputError("cache_miss: gamma and m disagree on the number of lists");
    if (m.empty()) throw InputError("cache_miss: empty capacity vector");

    std::vector<int> ma = m;
    ma[0] += 1;

    const T Em = cache_erec(gamma, m);
    if (Em == num_traits<T>::from_int(0))
        throw NumericError("cache_miss: the normalizing constant is zero");

    CacheMissResult<T> r;
    r.M = cache_erec(gamma, ma) / Em;
    if (lambda.empty()) return r;

    const std::size_t u = lambda.rows();
    const std::size_t n = gamma.rows();
    if (lambda.cols() != n)
        throw InputError("cache_miss: lambda and gamma disagree on the number of items");

    r.pi0.assign(n, num_traits<T>::from_int(0));
    for (std::size_t k = 0; k < n; ++k)
        r.pi0[k] = cache_erec(detail::gamma_without_row(gamma, k), m) / Em;

    r.MU.assign(u, num_traits<T>::from_int(0));
    for (std::size_t v = 0; v < u; ++v)
        for (std::size_t k = 0; k < n; ++k) r.MU[v] += lambda(v, k) * r.pi0[k];

    r.MI.assign(n, num_traits<T>::from_int(0));
    for (std::size_t k = 0; k < n; ++k) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t v = 0; v < u; ++v) s += lambda(v, k);
        r.MI[k] = s * r.pi0[k];
    }
    return r;
}

/** Overload without request rates: only the global miss rate is defined. */
template <class T>
CacheMissResult<T> cache_miss(const Matrix<T>& gamma, const std::vector<int>& m) {
    return cache_miss(gamma, m, Matrix<T>());
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_MISS_H
