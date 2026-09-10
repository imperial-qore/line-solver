/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_COST_H
#define LINE_API_CACHE_COST_H

/**
 * Mean per-list storage cost of a cache with item sizes, and the screen for
 * promotion paths that storage cost caps make unreachable.
 *
 * Templated port of matlab/src/api/cache/cache_cost.m and
 * cache_cost_pathcheck.m, cross-checked against
 * jar/src/main/java/jline/api/cache/Cache_cost.java and
 * Cache_cost_pathcheck.java.
 *
 * K_j = sum_i sigma_i pi_ij is the expected storage cost of the items resident
 * in list j at steady state (Casale-Gast, IEEE/ACM Trans. Networking 29(2),
 * 2021, Sec. IX). Every operation is a field operation, so the exact
 * instantiation returns K as a rational.
 *
 * The path check exists because E(m,k) sums over every size-feasible state
 * while under RR-C(m) an item only reaches list j by being promoted one list
 * at a time from the miss list: a cap on an intermediate list can leave
 * size-feasible states unreachable, and E(m,k) then normalizes over states the
 * cache never visits. An EMPTY report is a necessary, not sufficient,
 * condition for the two sets to agree -- a list of capacity above one may
 * still be unreachable when its cap admits no combination containing the item.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_prob_erec.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** One (item, list, blocking list) triple, all 0-based. */
struct CacheBlockedPair {
    std::size_t item;
    std::size_t list;
    std::size_t blocking_list;
};

/**
 * @param gamma (n x h) access factors
 * @param m     (h) list capacities
 * @param sigma (n) per-item storage costs
 * @param k     (h) per-list storage cost caps; empty for none
 * @param pij   (n x (h+1)) occupancy, column 0 the miss probability; empty to recompute
 * @return (h) mean storage cost held by each list
 */
template <class T>
std::vector<T> cache_cost(const Matrix<T>& gamma, const std::vector<int>& m,
                          const std::vector<int>& sigma, const std::vector<int>& k,
                          const Matrix<T>& pij) {
    const std::size_t n = gamma.rows();
    const std::size_t h = gamma.cols();
    if (sigma.size() != n)
        throw InputError("cache_cost: the item size vector must have one entry per item");
    const Matrix<T> occ = pij.empty() ? cache_prob_erec(gamma, m, sigma, k) : pij;
    std::vector<T> K(h, num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < h; ++j) {
        T c = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < n; ++i)
            c += num_traits<T>::from_int(static_cast<long>(sigma[i])) * occ(i, j + 1);
        K[j] = c;
    }
    return K;
}

template <class T>
std::vector<T> cache_cost(const Matrix<T>& gamma, const std::vector<int>& m,
                          const std::vector<int>& sigma, const std::vector<int>& k) {
    return cache_cost(gamma, m, sigma, k, Matrix<T>());
}

/**
 * @param gamma  (n x h) access factors
 * @param sigma  (n) per-item storage costs
 * @param k      (h) per-list storage cost caps
 * @param parent (h) parent list of each list, 0-based, -1 for lists rooted in the miss list
 * @return the blocked (item, list, blocking list) triples, empty when none
 */
template <class T>
std::vector<CacheBlockedPair> cache_cost_pathcheck(const Matrix<T>& gamma,
                                                   const std::vector<int>& sigma,
                                                   const std::vector<int>& k,
                                                   const std::vector<int>& parent) {
    const std::size_t n = gamma.rows();
    const std::size_t h = gamma.cols();
    std::vector<CacheBlockedPair> viol;
    if (sigma.size() != n || k.size() != h || parent.size() != h) return viol;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < h; ++j) {
            if (gamma(i, j) == num_traits<T>::from_int(0) || sigma[i] > k[j])
                continue;  // item i never resides in list j anyway
            int l = parent[j];
            while (l >= 0) {
                if (sigma[i] > k[static_cast<std::size_t>(l)]) {
                    CacheBlockedPair p;
                    p.item = i;
                    p.list = j;
                    p.blocking_list = static_cast<std::size_t>(l);
                    viol.push_back(p);
                    break;
                }
                l = parent[static_cast<std::size_t>(l)];
            }
        }
    }
    return viol;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_COST_H
