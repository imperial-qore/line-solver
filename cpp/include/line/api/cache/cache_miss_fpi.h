/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_MISS_FPI_H
#define LINE_API_CACHE_MISS_FPI_H

/**
 * Cache miss rates from the fixed-point multipliers.
 *
 * Templated port of matlab/src/api/cache/cache_miss_fpi.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_miss_fpi.java.
 *
 * With xi from cache_xi_fp and S(i) = sum_l gamma(i,l) xi(l), the probability
 * that item i is absent from the cache is pi0(i) = 1/(1+S(i)), so
 *
 *   MI(i) = (sum_v lambda(v,i)) pi0(i),   MU(v) = sum_i lambda(v,i) pi0(i),
 *   M     = sum_i MI(i).
 *
 * ARITHMETIC: transcendental, inherited from cache_xi_fp.
 *
 * Note the miss probability used here, 1/(1+S), is the correct per-item form;
 * it is NOT the 1 - sum_l pij(i,l) returned by cache_xi_fp, which is floored at
 * 1e-14. The two agree up to that floor.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_miss.h"
#include "line/api/cache/cache_xi_fp.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/**
 * @param gamma  (n x h) access factors
 * @param m      (h) list capacities
 * @param lambda (u x n) per-user per-item request rates
 */
template <class T>
CacheMissResult<T> cache_miss_fpi(const Matrix<T>& gamma, const std::vector<int>& m,
                                  const Matrix<T>& lambda) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_miss_fpi requires transcendental arithmetic");
    const std::size_t n = gamma.rows();
    const std::size_t h = gamma.cols();
    if (lambda.cols() != n)
        throw InputError("cache_miss_fpi: lambda and gamma disagree on the number of items");
    const std::size_t u = lambda.rows();

    const CacheXiFpResult<T> f = cache_xi_fp(gamma, m);
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    CacheMissResult<T> r;
    r.pi0.assign(n, zero);
    r.MI.assign(n, zero);
    r.MU.assign(u, zero);
    r.M = zero;

    for (std::size_t i = 0; i < n; ++i) {
        T s = zero;
        for (std::size_t l = 0; l < h; ++l) s += gamma(i, l) * f.xi[l];
        r.pi0[i] = one / (one + s);
        T lam = zero;
        for (std::size_t v = 0; v < u; ++v) lam += lambda(v, i);
        r.MI[i] = lam * r.pi0[i];
        r.M += r.MI[i];
        for (std::size_t v = 0; v < u; ++v) r.MU[v] += lambda(v, i) * r.pi0[i];
    }
    return r;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_MISS_FPI_H
