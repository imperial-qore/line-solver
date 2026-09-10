/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_MVA_H
#define LINE_API_CACHE_MVA_H

/**
 * Exact mean value analysis of a multi-list cache.
 *
 * Templated port of matlab/src/api/cache/cache_mva.m, cross-checked against
 * jar/src/main/java/jline/api/cache/Cache_mva.java.
 *
 * Sweeps the box of cache occupancies {0..m(1)} x ... x {0..m(h)} in the order
 * MATLAB's State.cartesian produces (first list varying fastest, so the
 * predecessor state m - e_l always precedes m). At each state the MVA
 * arrival-theorem step is
 *
 *   x(l)      = m(l) / sum_k gamma(k,l) (1 - pi(k; m - e_l))
 *   pij(k,l)  = gamma(k,l) (1 - pi(k; m - e_l)) x(l)
 *   pi(k)     = sum_l pij(k,l),
 *
 * i.e. the list-l throughput normalizes the "item k is not already cached"
 * probabilities seen at the previous population. Only divisions and products,
 * so the exact instantiation carries the whole sweep in rationals; the
 * conservation law sum_l pij(k,l) + pi0(k) = 1 then holds identically.
 *
 * REFERENCE DEFECT (both codebases): the normalizing constant E of the return
 * list is initialized to 1 and never updated -- MATLAB sets `E=1;` before the
 * sweep and returns it, and the JAR mirrors that with `int E = 1;`. The value
 * is therefore not the cache normalizing constant (use cache_erec for that).
 * It is reproduced here for interface compatibility and flagged in the result
 * struct, not silently recomputed, since callers may rely on the constant.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Return value of cache_mva, mirroring [pi,pi0,pij,x,u,E]. */
template <class T>
struct CacheMvaResult {
    std::vector<T> pi;   ///< (n) probability that item k is cached anywhere
    std::vector<T> pi0;  ///< (n) miss probability, 1 - pi
    Matrix<T> pij;       ///< (n x h) probability that item k sits in list l
    std::vector<T> x;    ///< (h) per-list throughput
    Matrix<T> u;         ///< (n x h) utilization, x(l) gamma(k,l)
    T E;                 ///< always 1: see the reference defect noted above
};

/**
 * @param gamma (n x h) access factors
 * @param m     (h) list capacities
 */
template <class T>
CacheMvaResult<T> cache_mva(const Matrix<T>& gamma, const std::vector<int>& m) {
    const std::size_t n = gamma.rows();
    const std::size_t h = gamma.cols();
    if (h != m.size()) throw InputError("cache_mva: gamma and m disagree on the number of lists");
    for (int v : m)
        if (v < 0) throw InputError("cache_mva: negative capacity");

    // Mixed-radix enumeration of the occupancy box, list 0 varying fastest.
    std::vector<std::size_t> stride(h, 1);
    std::size_t total = 1;
    for (std::size_t l = 0; l < h; ++l) {
        stride[l] = total;
        total *= static_cast<std::size_t>(m[l]) + 1;
    }

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    // pi(s,k) over the whole box; pij only needs the current and final states,
    // but is kept per-state to mirror the reference exactly.
    Matrix<T> pi(total, n, zero);
    std::vector<Matrix<T>> pij(total, Matrix<T>(n, h, zero));
    std::vector<T> x(h, zero);

    std::vector<int> mcur(h, 0);
    for (std::size_t s = 0; s < total; ++s) {
        for (std::size_t l = 0; l < h; ++l) {
            if (mcur[l] == 0) continue;  // m - e_l is outside the box
            const std::size_t s_l = s - stride[l];
            T den = zero;
            for (std::size_t k = 0; k < n; ++k) den += gamma(k, l) * (one - pi(s_l, k));
            if (den == zero)
                throw NumericError(
                    "cache_mva: list has zero aggregate access factor, the throughput is undefined");
            x[l] = num_traits<T>::from_int(static_cast<long>(mcur[l])) / den;
            for (std::size_t k = 0; k < n; ++k) {
                const T p = gamma(k, l) * (one - pi(s_l, k)) * x[l];
                pij[s](k, l) = p;
                pi(s, k) += p;
            }
        }
        // advance the mixed-radix counter
        for (std::size_t l = 0; l < h; ++l) {
            if (++mcur[l] <= m[l]) break;
            mcur[l] = 0;
        }
    }

    const std::size_t sfin = total - 1;  // the state equal to m
    CacheMvaResult<T> r;
    r.pi.assign(n, zero);
    r.pi0.assign(n, zero);
    for (std::size_t k = 0; k < n; ++k) {
        r.pi[k] = pi(sfin, k);
        r.pi0[k] = one - r.pi[k];
    }
    r.pij = pij[sfin];
    r.x = x;
    r.u = Matrix<T>(n, h, zero);
    for (std::size_t l = 0; l < h; ++l)
        for (std::size_t k = 0; k < n; ++k) r.u(k, l) = x[l] * gamma(k, l);
    r.E = one;
    return r;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_MVA_H
