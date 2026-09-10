/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_TTL_LRUM_MAP_H
#define LINE_API_CACHE_TTL_LRUM_MAP_H

/**
 * TTL approximation of an LRU(m) cache whose items are requested by Markovian
 * arrival processes.
 *
 * Templated port of matlab/src/api/cache/cache_ttl_lrum_map.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_ttl_lrum_map.java. Solves the
 * characteristic times with cache_t_lrum_map and then reports, per item, the
 * request-weighted hit probabilities of each list and the time-stationary
 * level occupancies (Gast and Van Houdt, Performance Evaluation 2017).
 *
 * This is the model to use when the items have genuinely distinct or
 * correlated request processes; when items are i.i.d. marks of a common stream
 * the request sequence is IRM and the Poisson-based approximations
 * (cache_ttl_lrua, cache_ttl_hlru) already apply.
 *
 * ARITHMETIC: transcendental, through cache_lrum_map_levelstats.
 *
 * The miss column follows MATLAB exactly, pij(k,0) = max(0, 1 - sum_l
 * hitfrac(k,l)): the hit fractions come from a fixed point and can overshoot
 * one by a rounding-scale amount, and the clamp is part of the reference
 * definition rather than a workaround added here.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_lrum_map_levelstats.h"
#include "line/api/cache/cache_t_lrum_map.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Request-weighted and time-stationary level probabilities. */
template <class T>
struct CacheTtlLrumMapResult {
    Matrix<T> pij;      ///< (n x h+1) column 0 = miss probability, column 1+l = hit in list l
    Matrix<T> pijtime;  ///< (n x h+1) time-stationary level occupancy probabilities
    std::vector<T> t;   ///< (h) characteristic times used
};

/**
 * @param items (n) per-item request MAPs
 * @param m     (h) list capacities
 * @param tol   relative tolerance passed to cache_t_lrum_map
 */
template <class T>
CacheTtlLrumMapResult<T> cache_ttl_lrum_map(const std::vector<mam::Map<T>>& items,
                                            const std::vector<T>& m, const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_ttl_lrum_map requires transcendental arithmetic");
    const std::size_t n = items.size();
    const std::size_t h = m.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    CacheTtlLrumMapResult<T> out;
    out.t = cache_t_lrum_map(items, m, tol);
    out.pij = Matrix<T>(n, h + 1, zero);
    out.pijtime = Matrix<T>(n, h + 1, zero);

    for (std::size_t k = 0; k < n; ++k) {
        const CacheLrumMapLevelStats<T> st =
            cache_lrum_map_levelstats(items[k].D0, items[k].D1, out.t);
        for (std::size_t l = 0; l <= h; ++l) out.pijtime(k, l) = st.prob[l];
        T hits = zero;
        for (std::size_t l = 0; l < h; ++l) {
            out.pij(k, 1 + l) = st.hitfrac[l];
            hits += st.hitfrac[l];
        }
        const T miss = one - hits;
        out.pij(k, 0) = miss > zero ? miss : zero;
    }
    return out;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_TTL_LRUM_MAP_H
