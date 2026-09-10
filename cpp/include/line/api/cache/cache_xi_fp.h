/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_XI_FP_H
#define LINE_API_CACHE_XI_FP_H

/**
 * Lagrange multipliers of a multi-list cache by fixed-point iteration.
 *
 * Templated port of matlab/src/api/cache/cache_xi_fp.m, cross-checked against
 * jar/src/main/java/jline/api/cache/Cache_xi_fp.java.
 *
 * The asymptotic (large-cache) form of the list occupancy distribution is
 *
 *   pij(k,l) = gamma(k,l) xi(l) / (1 + sum_s gamma(k,s) xi(s)),
 *   pi0(k)   = 1 - sum_l pij(k,l),
 *
 * and the multipliers xi are fixed by the capacity constraints, which in this
 * decoupled form read xi(l) = m(l) / sum_k pi0(k) gamma(k,l). Alternating the
 * two updates from pi0 = 1/(h+1) is the iteration below; it stops when the
 * relative change of pi0 falls under 1e-14 in the 1-norm.
 *
 * ARITHMETIC: every update is a field operation, but the iteration is stopped
 * by a tolerance and the fixed point is not reached in finitely many exact
 * steps, so the algorithm is inexact by nature. It is gated on transcendental
 * arithmetic rather than offered at exact arithmetic, where a "1e-14" cutoff
 * would be a rounding-free computation of a rounded answer.
 *
 * REFERENCE DEFECT (both codebases): the optional initial-xi argument is
 * inert. MATLAB writes `xi=zeros(1,h);` before testing `nargin<3`, discarding
 * the caller's vector, and in any case the first statement of the loop
 * recomputes xi from pi0, so no initial xi can influence the result. The JAR
 * has the same structure. The argument is therefore not offered here.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Return value of cache_xi_fp, mirroring [xi,pi0,pij,it]. */
template <class T>
struct CacheXiFpResult {
    std::vector<T> xi;   ///< (h) Lagrange multipliers
    std::vector<T> pi0;  ///< (n) per-item miss probability
    Matrix<T> pij;       ///< (n x h) per-item per-list hit probability
    int it;              ///< iterations performed
};

/**
 * @param gamma (n x h) access factors
 * @param m     (h) list capacities
 */
template <class T>
CacheXiFpResult<T> cache_xi_fp(const Matrix<T>& gamma, const std::vector<int>& m) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_xi_fp requires transcendental arithmetic");
    const std::size_t n = gamma.rows();
    const std::size_t h = gamma.cols();
    if (h != m.size()) throw InputError("cache_xi_fp: gamma and m disagree on the number of lists");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T tol = num_traits<T>::from_double(1e-14);

    std::vector<T> pi0(n, one / num_traits<T>::from_int(static_cast<long>(h + 1)));
    std::vector<T> xi(h, zero);
    Matrix<T> pij(n, h, zero);

    int it = 1;
    for (; it <= 10000; ++it) {
        const std::vector<T> pi0_1 = pi0;

        // xi(l) = m(l) / sum_k pi0(k) gamma(k,l)
        for (std::size_t l = 0; l < h; ++l) {
            T d = zero;
            for (std::size_t k = 0; k < n; ++k) d += pi0_1[k] * gamma(k, l);
            if (d == zero)
                throw NumericError("cache_xi_fp: a list has zero aggregate access factor");
            xi[l] = num_traits<T>::from_int(static_cast<long>(m[l])) / d;
        }

        for (std::size_t k = 0; k < n; ++k) {
            T s = zero;
            for (std::size_t l = 0; l < h; ++l) s += gamma(k, l) * xi[l];
            const T den = num_abs(T(one + s));
            T hits = zero;
            for (std::size_t l = 0; l < h; ++l) {
                pij(k, l) = num_abs(T(gamma(k, l) * xi[l])) / den;
                hits += pij(k, l);
            }
            const T v = one - hits;
            pi0[k] = v > tol ? v : tol;
        }

        T delta = zero;
        for (std::size_t k = 0; k < n; ++k) delta += num_abs(T(one - pi0[k] / pi0_1[k]));
        if (delta < tol) break;
    }
    if (it > 10000) it = 10000;

    for (std::size_t l = 0; l < h; ++l)
        if (xi[l] < zero) xi[l] = tol;

    CacheXiFpResult<T> r;
    r.xi = xi;
    r.pi0 = pi0;
    r.pij = pij;
    r.it = it;
    return r;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_XI_FP_H
