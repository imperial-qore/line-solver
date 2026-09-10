/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_PROB_SPM_H
#define LINE_API_CACHE_PROB_SPM_H

/**
 * Saddle-point approximation of the per-item cache hit probabilities.
 *
 * Templated port of matlab/src/api/cache/cache_prob_spm.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_prob_rayint.java (the JAR's
 * only implementation of this shape; see the note below).
 *
 * Same ratio-of-constants identity as cache_prob_erec,
 *
 *   prob(i,1+j) = m(j) gamma(i,j) exp(lE_i - lE),
 *
 * with lE the log normalizing constant from cache_spm and lE_i the same
 * quantity for the model with item i deleted and list j one slot smaller. The
 * miss probability closes the row. Because both constants are approximations
 * the rows only sum to one up to the saddle-point error, which is why the
 * reference wraps the miss entry in abs().
 *
 * ARITHMETIC: transcendental, inherited from cache_spm.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_spm.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/**
 * @param gamma (n x h) access factors
 * @param m     (h) list capacities
 * @return (n x (h+1)); column 0 miss, column 1+j hit in list j
 */
template <class T>
Matrix<T> cache_prob_spm(const Matrix<T>& gamma, const std::vector<int>& m) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_prob_spm requires transcendental arithmetic");
    using std::exp;
    const std::size_t n = gamma.rows();
    const std::size_t h = gamma.cols();
    if (h != m.size()) throw InputError("cache_prob_spm: gamma and m disagree on the list count");

    const T lE = cache_spm(gamma, m).lZ;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    Matrix<T> prob(n, h + 1, zero);
    for (std::size_t i = 0; i < n; ++i) {
        T hits = zero;
        for (std::size_t j = 0; j < h; ++j) {
            std::vector<int> mj = m;
            mj[j] -= 1;
            const T lEi = cache_spm(detail::gamma_without_row(gamma, i), mj).lZ;
            const T p =
                num_traits<T>::from_int(static_cast<long>(m[j])) * gamma(i, j) * exp(T(lEi - lE));
            prob(i, j + 1) = p;
            hits += p;
        }
        prob(i, 0) = num_abs(T(one - hits));
    }
    return prob;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_PROB_SPM_H
