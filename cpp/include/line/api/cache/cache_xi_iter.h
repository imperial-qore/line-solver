/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_XI_ITER_H
#define LINE_API_CACHE_XI_ITER_H

/**
 * Lagrange multipliers of a multi-list cache by the Gast-Van Houdt iteration.
 *
 * Templated port of matlab/src/api/cache/cache_xi_iter.m, cross-checked
 * against jar/src/main/java/jline/api/cache/Cache_xi_iter.java (and its
 * verbatim duplicate Cache_xi_bvh.java).
 *
 * Writing pp(0,k) = 1 and pp(l,k) = gamma(k,l-1), the stationary occupancy of
 * list l under the asymptotic independence approximation is
 *
 *   F_l(z_l) = sum_k z_l pp(l,k) / (n z_l pp(l,k) + a_l(k)),
 *   a_l(k)   = n sum_{s != l} z_s pp(s,k),
 *
 * which is increasing in z_l, so the capacity constraint F_l(z_l) = m(l)/n has
 * a unique root. The outer loop is a Gauss-Seidel sweep over the lists, each
 * inner solve a bracketed bisection (the bracket [0,1] when the unit point
 * already overshoots, otherwise doubled from 1 until it does), refined by a
 * fixed 50 halvings. Convergence of the sweep is declared when the multipliers
 * move by less than 1e-12 relative.
 *
 * ARITHMETIC: only field operations, but the answer is defined by two nested
 * tolerances -- a fixed 50-step bisection and a 1e-12 sweep test -- so it is
 * inexact by construction and is gated on transcendental arithmetic. Note in
 * particular that the 50 halvings cap the achievable accuracy at about 1e-15
 * of the bracket regardless of the precision of T.
 *
 * REFERENCE DEFECT (MATLAB): the third argument is declared as `tmax` and is
 * assigned Inf when absent, but is never read in the body, so it has no
 * effect. cache_spm passes its own `xi0` warm start into that slot, meaning
 * the warm start silently does nothing. The argument is not offered here.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

namespace detail {

/** sum_k z pp(l,k) / (n z pp(l,k) + a(k)), the occupancy of list l at z. */
template <class T>
T xi_iter_occupancy(const std::vector<T>& ppl, const std::vector<T>& a, const T& z, long n) {
    const T zero = num_traits<T>::from_int(0);
    const T nT = num_traits<T>::from_int(n);
    T s = zero;
    for (std::size_t k = 0; k < ppl.size(); ++k) {
        const T den = nT * z * ppl[k] + a[k];
        if (den == zero) throw NumericError("cache_xi_iter: singular occupancy denominator");
        s += z * ppl[k] / den;
    }
    return s;
}

}  // namespace detail

/**
 * @param gamma (n x h) access factors
 * @param m     (h) list capacities
 * @return (h) multipliers xi
 */
template <class T>
std::vector<T> cache_xi_iter(const Matrix<T>& gamma, const std::vector<int>& m) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_xi_iter requires transcendental arithmetic");
    const std::size_t n = gamma.rows();
    const std::size_t h = m.size();
    if (gamma.cols() != h)
        throw InputError("cache_xi_iter: gamma and m disagree on the number of lists");
    if (n == 0) throw InputError("cache_xi_iter: no items");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T nT = num_traits<T>::from_int(static_cast<long>(n));

    std::vector<T> f(h, zero);
    for (std::size_t l = 0; l < h; ++l) f[l] = num_traits<T>::from_int(static_cast<long>(m[l])) / nT;

    // pp: (h+1) x n, row 0 all ones, row l+1 the l-th column of gamma.
    Matrix<T> pp(h + 1, n, one);
    for (std::size_t l = 0; l < h; ++l)
        for (std::size_t k = 0; k < n; ++k) pp(l + 1, k) = gamma(k, l);

    std::vector<T> z(h + 1, one), zold(h + 1, zero);
    const T reltol = num_traits<T>::from_double(1e-12);

    for (int sweep = 0;; ++sweep) {
        T dmax = zero, omax = zero;
        for (std::size_t l = 0; l <= h; ++l) {
            const T d = num_abs(T(z[l] - zold[l]));
            if (d > dmax) dmax = d;
            const T o = num_abs(zold[l]);
            if (o > omax) omax = o;
        }
        if (!(dmax > reltol * omax)) break;
        if (sweep > 10000)
            throw NumericError("cache_xi_iter: the Gauss-Seidel sweep did not converge");
        zold = z;

        // temp(k) = n sum_s z(s) pp(s,k)
        std::vector<T> temp(n, zero);
        for (std::size_t k = 0; k < n; ++k) {
            T s = zero;
            for (std::size_t l = 0; l <= h; ++l) s += z[l] * pp(l, k);
            temp[k] = nT * s;
        }

        for (std::size_t l = 0; l < h; ++l) {
            std::vector<T> ppl(n), a(n);
            for (std::size_t k = 0; k < n; ++k) {
                ppl[k] = pp(l + 1, k);
                a[k] = temp[k] - nT * z[l + 1] * ppl[k];
            }

            const T Fi = detail::xi_iter_occupancy(ppl, a, one, static_cast<long>(n));
            T zmin, zmax;
            if (Fi > f[l]) {
                zmin = zero;
                zmax = one;
            } else {
                zmin = one;
                zmax = two;
                while (detail::xi_iter_occupancy(ppl, a, zmax, static_cast<long>(n)) < f[l]) {
                    zmin = zmax;
                    zmax = zmax * two;
                }
            }
            for (int b = 0; b < 50; ++b) {
                const T mid = (zmin + zmax) / two;
                z[l + 1] = mid;
                if (detail::xi_iter_occupancy(ppl, a, mid, static_cast<long>(n)) < f[l])
                    zmin = mid;
                else
                    zmax = mid;
            }
        }
    }

    return std::vector<T>(z.begin() + 1, z.end());
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_XI_ITER_H
