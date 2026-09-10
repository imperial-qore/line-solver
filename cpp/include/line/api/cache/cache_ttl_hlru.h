/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_TTL_HLRU_H
#define LINE_API_CACHE_TTL_HLRU_H

/**
 * TTL (characteristic-time) approximation of an h-LRU / LRU(m) cache.
 *
 * Templated port of matlab/src/api/cache/cache_ttl_hlru.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_ttl_hlru.java.
 *
 * The policy is h LRU lists of capacities m(1..h): a miss inserts the item at
 * the head of list 1 and a hit in list l exchanges it with the tail of list
 * l+1. Under the characteristic-time approximation the level of an item with
 * request rate lam is a birth-death chain with up-probability 1 - e(l) and
 * down-probability e(l), e(l) = exp(-lam T(l)), so
 *
 *   pi(l) proportional to prod_{s<=l} (1 - e(s))/e(s),
 *
 * with the times T solved from the capacity constraints by cache_t_hlru. For
 * h = 1 this is exactly the Che approximation of LRU, so an M/LRU/1 cache can
 * be checked against the closed form 1 - exp(-lam T) directly.
 *
 * ARITHMETIC: transcendental, as cache_t_hlru.
 *
 * The MATLAB reference takes lambda as the (u x n x h+1) array built by
 * solver_mva_cache_analyzer and sums slice min(2,h+1) over the user classes;
 * that slice carries the same per-item rate as every other one. The port takes
 * that slice directly as a (u x n) matrix, which is the same reduction without
 * the three-dimensional container.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_t_hlru.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/**
 * @param lambda (u x n) per-user per-item request rates
 * @param m      (h) list capacities
 * @return (n x (h+1)); column 0 is "not cached", column 1+l is "in list l"
 */
template <class T>
Matrix<T> cache_ttl_hlru(const Matrix<T>& lambda, const std::vector<int>& m) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_ttl_hlru requires transcendental arithmetic");
    if (lambda.empty()) throw InputError("cache_ttl_hlru: empty request-rate matrix");
    const std::size_t n = lambda.cols();
    std::vector<T> lam(n, num_traits<T>::from_int(0));
    for (std::size_t v = 0; v < lambda.rows(); ++v)
        for (std::size_t k = 0; k < n; ++k) lam[k] += lambda(v, k);

    const std::vector<T> Tv = detail::hlru_solve_times(lam, m);
    return detail::hlru_levelprobs(lam, Tv);
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_TTL_HLRU_H
