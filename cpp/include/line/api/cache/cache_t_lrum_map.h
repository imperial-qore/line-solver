/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_T_LRUM_MAP_H
#define LINE_API_CACHE_T_LRUM_MAP_H

/**
 * Characteristic times of the LRU(m)-MAP TTL approximation.
 *
 * Templated port of matlab/src/api/cache/cache_t_lrum_map.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_t_lrum_map.java. The times
 * T_1..T_h are fixed by equating the expected occupancy of each list to its
 * capacity,
 *
 *   sum_k occ_l(item k; T) = m_l,   l = 1..h,
 *
 * with occ from cache_lrum_map_levelstats (Gast and Van Houdt, Performance
 * Evaluation 2017, Section 3.1.2).
 *
 * HOW THE SYSTEM IS SOLVED, and why it is not a multivariate optimizer:
 * occ_l is strictly increasing in T_l (a longer timer keeps the item in list l
 * longer), so each equation is a well-posed SCALAR root problem in its own
 * unknown once the other times are held fixed. The system is therefore solved
 * by Gauss-Seidel sweeps of bracketed scalar solves -- the same structure
 * cache_t_hlru.h already uses for h-LRU -- with the bracket found by doubling
 * and closed by bisection. The result is deterministic and needs no derivative,
 * no line search and no trust region.
 *
 * That is a real improvement on both references. MATLAB calls fsolve on
 * log(T) from a fixed zero start, so it needs the Optimization Toolbox and
 * stops on fsolve's default tolerances (the reference instance in the tests
 * leaves a capacity residual of about 1e-10). The JAR does not solve the
 * system at all: it hands the residual NORM to COBYLA, a derivative-free
 * constrained optimizer, which turns h independent monotone equations into one
 * nonconvex minimization -- and COBYLA is stopped at rhoend = 1e-6, so its
 * times carry a much larger error than either the MATLAB or this port.
 *
 * ARITHMETIC: transcendental, through cache_lrum_map_levelstats.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_lrum_map_levelstats.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/rootfind.h"

namespace line {
namespace cache {

namespace detail {

/** Total occupancy of every list at the characteristic times Tv. */
template <class T>
std::vector<T> lrum_map_occupancy(const std::vector<mam::Map<T>>& items, const std::vector<T>& Tv) {
    std::vector<T> occ(Tv.size(), num_traits<T>::from_int(0));
    for (std::size_t k = 0; k < items.size(); ++k) {
        const CacheLrumMapLevelStats<T> st =
            cache_lrum_map_levelstats(items[k].D0, items[k].D1, Tv);
        for (std::size_t l = 0; l < Tv.size(); ++l) occ[l] += st.occ[l];
    }
    return occ;
}

}  // namespace detail

/**
 * @param items  (n) per-item request MAPs
 * @param m      (h) list capacities, 0 < sum(m) < n
 * @param tol    relative tolerance on the times, e.g. 1e-12
 * @param maxswp cap on Gauss-Seidel sweeps
 * @return (h) characteristic times
 */
template <class T>
std::vector<T> cache_t_lrum_map(const std::vector<mam::Map<T>>& items, const std::vector<T>& m,
                                const T& tol, unsigned maxswp = 200) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_t_lrum_map requires transcendental arithmetic");
    const std::size_t n = items.size();
    const std::size_t h = m.size();
    if (n == 0) throw InputError("cache_t_lrum_map: no items");
    if (h == 0) throw InputError("cache_t_lrum_map: no lists");
    const T zero = num_traits<T>::from_int(0);
    T total = zero;
    for (std::size_t l = 0; l < h; ++l) {
        if (!(m[l] > zero)) throw InputError("cache_t_lrum_map: list capacities must be positive");
        total += m[l];
    }
    if (!(total < num_traits<T>::from_int(static_cast<long>(n))))
        throw InputError("cache_t_lrum_map: the cache is not smaller than the item catalogue");

    std::vector<T> Tv(h, num_traits<T>::from_int(1));
    for (unsigned sweep = 0; sweep < maxswp; ++sweep) {
        const std::vector<T> Told = Tv;
        for (std::size_t l = 0; l < h; ++l) {
            std::vector<T> work = Tv;
            // Residual of list l as a function of its own characteristic time,
            // increasing, negative at zero.
            auto resid = [&](const T& x) {
                work[l] = x;
                return T(detail::lrum_map_occupancy(items, work)[l] - m[l]);
            };
            T lo = num_traits<T>::from_double(1e-12);
            T hi = Tv[l] > lo ? Tv[l] : num_traits<T>::from_int(1);
            if (!(resid(lo) < zero))
                throw NumericError("cache_t_lrum_map: list occupancy exceeds its capacity even at "
                                   "a vanishing characteristic time");
            bracket_expand<T>(resid, lo, hi, 200);
            const RootResult<T> r = root_bisect<T>(resid, lo, hi, T(tol * hi), 400);
            Tv[l] = r.root;
        }
        T rel = zero;
        for (std::size_t l = 0; l < h; ++l) {
            const T den = Told[l] > tol ? Told[l] : tol;
            const T d = num_abs(T(Tv[l] - Told[l])) / den;
            if (d > rel) rel = d;
        }
        if (rel < tol) break;
    }
    return Tv;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_T_LRUM_MAP_H
