/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_MISS_ASY_H
#define LINE_API_CACHE_MISS_ASY_H

/**
 * Asymptotic (large-cache) miss ratio by a rank-threshold fixed point.
 *
 * Templated port of jar/src/main/java/jline/api/cache/Cache_miss_asy.java.
 * MATLAB has no counterpart, so the JAR is the reference.
 *
 * The deterministic limit of a multi-list cache: as the item count grows, list
 * l holds exactly the m(l) items of largest effective popularity, so an item's
 * membership becomes a THRESHOLD test rather than a probability. Writing
 * pi(k) for the miss probability of item k, the effective popularity of item j
 * in list l is gamma(l,j) (1 - pi(j)) and the fixed point is
 *
 *   pi(k) = sum_l gamma(l,k) 1{item k is outside the top m(l)} / sum_l gamma(l,k),
 *
 * iterated to a sup-norm tolerance from the uniform start pi = 1/n. The
 * returned scalar is the request-weighted miss ratio
 * sum_{l,k} gamma(l,k) pi(k) / sum_{l,k} gamma(l,k).
 *
 * INDEX CONVENTION, AND IT IS THE REVERSE OF EVERY OTHER CACHE FUNCTION HERE.
 * The reference reads n = gamma.getNumCols() and h = gamma.getNumRows(), so
 * its gamma is (h x n), LIST-major, whereas cache_spm, cache_erec, cache_miss
 * and the rest all take gamma as (n x h), ITEM-major. That is reproduced,
 * because silently transposing would make the two conventions disagree about
 * which of a non-square gamma's dimensions is the item count and return a
 * plausible wrong ratio rather than an error. Callers holding an item-major
 * gamma must transpose before calling.
 *
 * The threshold is strict (`>`), so ties at the cutoff resolve as "in cache".
 * With the ranking taken over the OTHER n-1 items and then compared against
 * item k, an item exactly at the boundary is admitted; this reproduces the
 * reference and matters only on gamma matrices with repeated entries.
 *
 * A degenerate capacity (zero total, or any negative entry) returns 1, i.e.
 * every request misses, which is the reference's early exit.
 *
 * Arithmetic: EXACT-CAPABLE in its operations, but the answer is defined by a
 * sup-norm tolerance and an iteration cap, so it is inexact by construction
 * and gated on transcendental arithmetic like the other fixed-point cache
 * routines (cache_xi_iter, cache_miss_fpi).
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/**
 * @param gamma     (h x n) LIST-major access factors; see the note above
 * @param m         (h) list capacities
 * @param maxIter   cap on fixed-point sweeps
 * @param tolerance sup-norm stopping tolerance on pi
 * @return the request-weighted asymptotic miss ratio
 */
template <class T>
T cache_miss_asy(const Matrix<T>& gamma, const std::vector<int>& m, int maxIter,
                 const T& tolerance) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_miss_asy is a tolerance-stopped fixed point and needs transcendental "
                  "arithmetic");
    const std::size_t h = gamma.rows();
    const std::size_t n = gamma.cols();
    if (m.size() != h)
        throw InputError("cache_miss_asy: gamma is list-major and disagrees with m on the list "
                         "count");
    if (n == 0) throw InputError("cache_miss_asy: no items");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    long mtot = 0;
    for (std::size_t l = 0; l < h; ++l) {
        if (m[l] < 0) return one;
        mtot += m[l];
    }
    if (mtot == 0) return one;

    // GlobalConstants.Zero, the reference's test for a vanishing rate sum
    const T tiny = num_traits<T>::from_double(1e-14);

    std::vector<T> pi(n, T(one / num_traits<T>::from_int(static_cast<long>(n))));
    std::vector<T> prev(n, zero);
    std::vector<T> next(n, zero);
    std::vector<T> pop;  // effective popularities of the other items, reused
    pop.reserve(n);

    for (int iter = 0; iter < maxIter; ++iter) {
        prev = pi;
        for (std::size_t k = 0; k < n; ++k) {
            T numer = zero, denom = zero;
            for (std::size_t l = 0; l < h; ++l) {
                const int cap = m[l];
                if (cap <= 0) continue;
                pop.clear();
                for (std::size_t j = 0; j < n; ++j)
                    if (j != k) pop.push_back(T(gamma(l, j) * T(one - prev[j])));
                // descending, so the take-th entry is the weakest still cached
                std::sort(pop.begin(), pop.end(), [](const T& a, const T& b) { return b < a; });
                const std::size_t take =
                    std::min(static_cast<std::size_t>(cap), pop.size());

                T notInCache = one;
                if (take < static_cast<std::size_t>(cap)) {
                    // fewer competitors than slots: item k is always cached
                    notInCache = zero;
                } else if (take > 0) {
                    const T weakest = pop[take - 1];
                    if (T(gamma(l, k) * T(one - prev[k])) > weakest) notInCache = zero;
                }
                numer += gamma(l, k) * notInCache;
                denom += gamma(l, k);
            }
            next[k] = denom > tiny ? T(numer / denom) : one;
        }
        pi = next;

        T diff = zero;
        for (std::size_t k = 0; k < n; ++k) {
            const T d = num_abs(T(pi[k] - prev[k]));
            if (d > diff) diff = d;
        }
        if (diff < tolerance) break;
    }

    T missRate = zero, totalRate = zero;
    for (std::size_t l = 0; l < h; ++l) {
        for (std::size_t k = 0; k < n; ++k) {
            missRate += gamma(l, k) * pi[k];
            totalRate += gamma(l, k);
        }
    }
    return totalRate > tiny ? T(missRate / totalRate) : one;
}

/** Reference defaults: 1000 sweeps at a 1e-8 sup-norm tolerance. */
template <class T>
T cache_miss_asy(const Matrix<T>& gamma, const std::vector<int>& m) {
    return cache_miss_asy(gamma, m, 1000, num_traits<T>::from_double(1e-8));
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_MISS_ASY_H
