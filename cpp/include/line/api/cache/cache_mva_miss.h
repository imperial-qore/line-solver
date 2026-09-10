/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_MVA_MISS_H
#define LINE_API_CACHE_MVA_MISS_H

/**
 * Per-item and global cache miss probabilities by mean value analysis.
 *
 * Templated port of matlab/src/api/cache/cache_mva_miss.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_mva_miss.java.
 *
 * Recursion on the capacity vector: at capacity m the weight of item k in list
 * j is w(k,j) = prod_{i<=j} R(i,k) p(k)^j |Mk(k; m - e_j)|, i.e. the item's
 * miss probability one slot down, scaled by the probability of reaching list j
 * through the routing chain R and by j independent requests. The per-list
 * normalization x(j) = 1/sum_k |w(k,j)| turns the weights into occupancy
 * shares and the miss probability of item k is what is left after removing the
 * m(j) slots of every list:
 *
 *   Mk(k) = |1 - sum_j x(j) m(j) w(k,j)|,   M = sum_k p(k) Mk(k).
 *
 * The base case sum(m) == 0 (or any negative capacity, which the recursion
 * reaches from a list of capacity zero) is Mk == 1: nothing is cached.
 *
 * Only products, sums, integer powers and divisions, so the algorithm stays in
 * the field and instantiates at exact arithmetic. The abs() calls of the
 * reference are kept: they are no-ops on a correct instance (every weight is
 * non-negative) but they change the answer on an over-committed one, so
 * dropping them would silently diverge from MATLAB and the JAR.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Return value of cache_mva_miss, mirroring [M,Mk]. */
template <class T>
struct CacheMvaMissResult {
    T M;                ///< global miss rate, sum_k p(k) Mk(k)
    std::vector<T> Mk;  ///< (n) per-item miss probability
};

/**
 * @param p (n) item popularities
 * @param m (h) list capacities
 * @param R (h x n) per-list routing probabilities
 */
template <class T>
CacheMvaMissResult<T> cache_mva_miss(const std::vector<T>& p, const std::vector<int>& m,
                                     const Matrix<T>& R) {
    const std::size_t n = p.size();
    const std::size_t h = m.size();
    if (R.rows() < h || R.cols() != n)
        throw InputError("cache_mva_miss: R must be at least (h x n)");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    long mt = 0;
    int mmin = 0;
    bool first = true;
    for (int v : m) {
        mt += v;
        if (first || v < mmin) mmin = v;
        first = false;
    }

    CacheMvaMissResult<T> r;
    if (mt == 0 || mmin < 0) {
        r.Mk.assign(n, one);
        r.M = zero;
        for (std::size_t k = 0; k < n; ++k) r.M += p[k];
        return r;
    }

    Matrix<T> w(n, h, zero);
    for (std::size_t j = 0; j < h; ++j) {
        std::vector<int> mj = m;
        mj[j] -= 1;
        const CacheMvaMissResult<T> rec = cache_mva_miss(p, mj, R);
        for (std::size_t k = 0; k < n; ++k) {
            T prod = one;
            for (std::size_t i = 0; i <= j; ++i) prod *= R(i, k);
            w(k, j) = prod * num_pow_int(p[k], static_cast<unsigned>(j + 1)) * num_abs(rec.Mk[k]);
        }
    }

    std::vector<T> x(h, zero);
    for (std::size_t j = 0; j < h; ++j) {
        T s = zero;
        for (std::size_t k = 0; k < n; ++k) s += num_abs(w(k, j));
        if (s == zero)
            throw NumericError("cache_mva_miss: list has zero total weight, x is undefined");
        x[j] = one / s;
    }

    r.Mk.assign(n, zero);
    r.M = zero;
    for (std::size_t k = 0; k < n; ++k) {
        T v = one;
        for (std::size_t j = 0; j < h; ++j)
            v -= x[j] * num_traits<T>::from_int(static_cast<long>(m[j])) * w(k, j);
        r.Mk[k] = num_abs(v);
        r.M += p[k] * r.Mk[k];
    }
    return r;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_MVA_MISS_H
