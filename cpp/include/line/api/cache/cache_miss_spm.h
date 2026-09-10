/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_MISS_SPM_H
#define LINE_API_CACHE_MISS_SPM_H

/**
 * Saddle-point approximation of the cache miss rates.
 *
 * Templated port of matlab/src/api/cache/cache_miss_spm.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_miss_rayint.java.
 *
 * Identical structure to cache_miss, with every exact constant replaced by the
 * cache_spm approximation of its logarithm:
 *
 *   M      = exp(lE(gamma, m + e_1) - lE(gamma, m)),
 *   pi0(k) = exp(lE(gamma without item k, m) - lE(gamma, m)),
 *   MU(v)  = sum_k lambda(v,k) pi0(k),   MI(k) = (sum_v lambda(v,k)) pi0(k).
 *
 * Items with an all-zero access-factor row are never cached and are skipped,
 * contributing pi0 = 0 and MI = 0 exactly as in the reference.
 *
 * ARITHMETIC: transcendental, inherited from cache_spm.
 *
 * REFERENCE DEFECT (MATLAB): the "recompute xi" fallback taken when a pi0(k)
 * falls outside [0,1] re-invokes cache_spm without the warm start, but the
 * warm start never had any effect in the first place -- cache_spm forwards it
 * to cache_xi_iter's third argument, which that function declares as `tmax`
 * and never reads. The two branches are therefore the same computation, so
 * this port evaluates lE1(k) once. This changes no value; it only removes a
 * duplicated call.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_miss.h"
#include "line/api/cache/cache_spm.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Return value of cache_miss_spm, mirroring [M,MU,MI,pi0,lE]. */
template <class T>
struct CacheMissSpmResult {
    T M;                 ///< global miss rate
    std::vector<T> MU;   ///< (u) per-user miss rate
    std::vector<T> MI;   ///< (n) per-item miss rate
    std::vector<T> pi0;  ///< (n) per-item miss probability
    T lE;                ///< log normalizing constant at capacity m
};

/**
 * @param gamma  (n x h) access factors
 * @param m      (h) list capacities
 * @param lambda (u x n) per-user per-item request rates
 */
template <class T>
CacheMissSpmResult<T> cache_miss_spm(const Matrix<T>& gamma, const std::vector<int>& m,
                                     const Matrix<T>& lambda) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_miss_spm requires transcendental arithmetic");
    using std::exp;
    if (m.empty()) throw InputError("cache_miss_spm: empty capacity vector");
    const std::size_t n = gamma.rows();
    const std::size_t h = gamma.cols();
    if (h != m.size()) throw InputError("cache_miss_spm: gamma and m disagree on the list count");

    std::vector<int> ma = m;
    ma[0] += 1;

    const T zero = num_traits<T>::from_int(0);
    CacheMissSpmResult<T> r;
    r.lE = cache_spm(gamma, m).lZ;
    const T lEa = cache_spm(gamma, ma).lZ;
    r.M = exp(T(lEa - r.lE));

    if (lambda.empty()) return r;
    if (lambda.cols() != n)
        throw InputError("cache_miss_spm: lambda and gamma disagree on the number of items");
    const std::size_t u = lambda.rows();

    r.pi0.assign(n, zero);
    r.MI.assign(n, zero);
    r.MU.assign(u, zero);

    for (std::size_t k = 0; k < n; ++k) {
        T rowsum = zero;
        for (std::size_t l = 0; l < h; ++l) rowsum += gamma(k, l);
        if (!(rowsum > zero)) continue;  // never cached: pi0 = 0, MI = 0

        const T lE1 = cache_spm(detail::gamma_without_row(gamma, k), m).lZ;
        r.pi0[k] = exp(T(lE1 - r.lE));
        T lam = zero;
        for (std::size_t v = 0; v < u; ++v) {
            r.MU[v] += lambda(v, k) * r.pi0[k];
            lam += lambda(v, k);
        }
        r.MI[k] = lam * r.pi0[k];
    }
    return r;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_MISS_SPM_H
