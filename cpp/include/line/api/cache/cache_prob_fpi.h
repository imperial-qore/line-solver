/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_PROB_FPI_H
#define LINE_API_CACHE_PROB_FPI_H

/**
 * Cache hit and miss probabilities from the fixed-point multipliers.
 *
 * Templated port of matlab/src/api/cache/cache_prob_fpi.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_prob_fpi.java.
 *
 * With S(i) = sum_l gamma(i,l) xi(l) from cache_xi_fp,
 *
 *   prob(i,1)   = 1 / (1 + S(i))          (miss)
 *   prob(i,1+l) = S(i) / (1 + S(i))       (hit, as written by the references)
 *
 * ARITHMETIC: transcendental, inherited from cache_xi_fp's tolerance-stopped
 * iteration.
 *
 * REFERENCE DEFECT (both codebases): the hit entry is the AGGREGATE hit
 * probability S/(1+S), written identically into all h list columns, instead of
 * the per-list gamma(i,l) xi(l)/(1+S(i)). MATLAB's
 * `prob(i,2:(1+h)) = gamma(i,:)*xi(:) ./ (1+gamma(i,:)*xi(:))` is a scalar
 * broadcast over the row, and the JAR reproduces it with an explicit loop that
 * stores the same `mul.get(0)` in every column. The consequence is that the
 * row sums to 1 + (h-1) S/(1+S), not to 1, whenever h > 1; only h == 1 is
 * correct. This port is faithful to the references -- correcting it here would
 * silently diverge from MATLAB and the JAR -- and the conservation test in the
 * suite is asserted only for h == 1 for exactly this reason.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_xi_fp.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/**
 * @param gamma (n x h) access factors
 * @param m     (h) list capacities
 * @return (n x (h+1)); column 0 miss, columns 1..h hit (see the defect note)
 */
template <class T>
Matrix<T> cache_prob_fpi(const Matrix<T>& gamma, const std::vector<int>& m) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_prob_fpi requires transcendental arithmetic");
    const std::size_t n = gamma.rows();
    const std::size_t h = gamma.cols();
    const CacheXiFpResult<T> f = cache_xi_fp(gamma, m);

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    Matrix<T> prob(n, h + 1, zero);
    for (std::size_t i = 0; i < n; ++i) {
        T s = zero;
        for (std::size_t l = 0; l < h; ++l) s += gamma(i, l) * f.xi[l];
        const T den = one + s;
        prob(i, 0) = one / den;
        for (std::size_t l = 0; l < h; ++l) prob(i, 1 + l) = s / den;
    }
    return prob;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_PROB_FPI_H
