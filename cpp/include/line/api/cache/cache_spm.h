/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_SPM_H
#define LINE_API_CACHE_SPM_H

/**
 * Saddle-point approximation of the cache normalizing constant.
 *
 * Templated port of matlab/src/api/cache/cache_spm.m, cross-checked against
 * jar/src/main/java/jline/api/cache/Cache_spm.java.
 *
 * The constant E(gamma,m) of cache_erec is the coefficient extraction
 *
 *   E = [prod_l z_l^{m_l}] prod_k (1 + sum_l gamma(k,l) z_l) * prod_l m_l!,
 *
 * evaluated here by a multidimensional saddle point at the multipliers xi from
 * cache_xi_iter. With S(k) = sum_l gamma(k,l) xi(l) and
 *
 *   phi = sum_k log(1 + S(k)) - sum_l m(l) log(xi(l)),
 *   C(j,l) = delta(j,l) sum_k gamma(k,j)/(1+S(k))
 *            - xi(j) sum_k gamma(k,j) gamma(k,l)/(1+S(k))^2,
 *
 * the Gaussian integral around the saddle gives
 *
 *   log E ~ phi + sum_l log(m_l!) - (h/2) log(2 pi) - (1/2) sum_l log xi(l)
 *           - (1/2) log det C.
 *
 * TWO BOUNDARIES, both fixed in all four codebases on 2026-08-29; the notes
 * below say what the behaviour used to be, because saved results predating the
 * fix carry it.
 *
 * n == sum(m): every item is cached, so the capacity equations force every
 * multiplier to infinity and there is no interior saddle. Z comes from the
 * exact cache_erec and xi is reported as +infinity, its limit, WITHOUT running
 * the iteration. Before the fix MATLAB computed Z exactly but then called
 * cache_xi_iter for the third return value and hung (`cache_spm([.5 .25;.4
 * .2;.3 .15],[2 1])` did not return in two minutes); the JAR had no fallback at
 * all and ran the same non-terminating loop for every value; this port raised
 * NumericError from cache_xi_iter's sweep cap, losing the exact Z with it.
 *
 * m_l == 0: list l has xi(l)=0, which is a boundary of the Laplace integral
 * rather than a direction of it, so list l is dropped before the saddle solve
 * and h shrinks with it. Dropping is exact: setting z_l=0 in the generating
 * function removes list l from E(m), and prod_l m_l! is unchanged because
 * 0!=1. Before the fix no codebase dropped it -- only all-zero ROWS of gamma
 * were filtered, never zero-capacity COLUMNS -- so the bisection floored xi(l)
 * at ~2^-50 and the -(1/2) sum_l log xi(l) prefactor gained ~+17 per empty
 * list, silently: n=12 items at gamma=(0.8,0.6,0.4) and m=(0,3,3) returned
 * lZ=24.814 against an exact 9.127, a factor of 6.5e8. cache_prob_spm reaches
 * this on any list with m_l==1, since it evaluates E at oner(m,l).
 * All m_l zero is the empty cache: Z=1, lZ=0, xi=0.
 *
 * ARITHMETIC: log, exp and sqrt throughout, plus the tolerance-stopped
 * cache_xi_iter, so transcendental arithmetic is required.
 *
 * MATLAB's `lZ=real(lZ)` discards the imaginary part that appears when det C
 * or some xi(l) is negative; that is implemented here as taking the modulus
 * inside the logarithm, which is the same real branch. Z is returned on the
 * same branch (MATLAB would return a complex Z there).
 *
 * The item count is gamma's ROW count in every codebase. MATLAB used to read
 * `n = length(gamma)`, which is max(n,h), and used it both for the n == mt
 * degenerate test and as the loop bound over gamma's rows, reading past the end
 * whenever h > n; that was corrected to size(gamma,1) alongside the two fixes
 * above. cache_erec still carries the same length() misuse.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_xi_iter.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace cache {

/** Return value of cache_spm, mirroring [Z,lZ,xi]. */
template <class T>
struct CacheSpmResult {
    T Z;                ///< normalizing constant
    T lZ;               ///< its logarithm
    std::vector<T> xi;  ///< (h) saddle-point multipliers
};

namespace detail {

/** Determinant by Gaussian elimination with partial pivoting; 0 if singular. */
template <class T>
T cache_det(Matrix<T> A) {
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("cache_spm: determinant of a non-square matrix");
    const T zero = num_traits<T>::from_int(0);
    T det = num_traits<T>::from_int(1);
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t p = k;
        T amax = num_abs(A(k, k));
        for (std::size_t i = k + 1; i < n; ++i) {
            const T a = num_abs(A(i, k));
            if (a > amax) {
                amax = a;
                p = i;
            }
        }
        if (amax == zero) return zero;
        if (p != k) {
            for (std::size_t j = 0; j < n; ++j) std::swap(A(k, j), A(p, j));
            det = -det;
        }
        det *= A(k, k);
        for (std::size_t i = k + 1; i < n; ++i) {
            const T fct = A(i, k) / A(k, k);
            for (std::size_t j = k; j < n; ++j) A(i, j) -= fct * A(k, j);
        }
    }
    return det;
}

/** gamma restricted to the rows with a strictly positive row sum. */
template <class T>
Matrix<T> gamma_nonzero_rows(const Matrix<T>& gamma) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < gamma.rows(); ++i) {
        T s = zero;
        for (std::size_t j = 0; j < gamma.cols(); ++j) s += gamma(i, j);
        if (s > zero) keep.push_back(i);
    }
    Matrix<T> g(keep.size(), gamma.cols());
    for (std::size_t a = 0; a < keep.size(); ++a)
        for (std::size_t j = 0; j < gamma.cols(); ++j) g(a, j) = gamma(keep[a], j);
    return g;
}

}  // namespace detail

/**
 * @param gamma_in (n x h) access factors; all-zero rows are dropped first
 * @param m     (h) list capacities
 */
template <class T>
CacheSpmResult<T> cache_spm(const Matrix<T>& gamma_in, const std::vector<int>& m) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_spm requires transcendental arithmetic");
    using std::acos;
    using std::log;
    using std::exp;
    using std::sqrt;

    const Matrix<T> gamma = detail::gamma_nonzero_rows(gamma_in);
    const std::size_t h = m.size();
    const std::size_t n = gamma.rows();
    if (gamma.cols() != h) throw InputError("cache_spm: gamma and m disagree on the number of lists");

    long mt = 0;
    for (int v : m) mt += v;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T half = one / two;

    CacheSpmResult<T> r;
    if (static_cast<long>(n) == mt) {
        // Degenerate saddle: every item is cached, so the capacity equations force
        // every multiplier to infinity and cache_xi_iter cannot converge. Take Z
        // from the exact recursion and report the limit, rather than iterating.
        r.Z = cache_erec(gamma, m);
        r.lZ = log(r.Z);
        r.xi.assign(h, num_traits<T>::from_double(std::numeric_limits<double>::infinity()));
        return r;
    }

    // A list with no capacity has xi=0, which is a boundary of the Laplace integral
    // rather than a direction of it, so it must leave the expansion: kept, its
    // -sum_l log(sqrt(xi_l)) prefactor diverges and Z comes out far too large.
    // Dropping it is exact, since setting z_l=0 in the generating function removes
    // list l from E(m) and prod_l m_l! is unchanged because 0!=1.
    std::vector<std::size_t> keep;
    for (std::size_t l = 0; l < h; ++l)
        if (m[l] > 0) keep.push_back(l);
    const std::size_t hk = keep.size();
    if (hk == 0) {
        // Empty cache: E(0)=1 and prod_l m_l!=1.
        r.Z = one;
        r.lZ = zero;
        r.xi.assign(h, zero);
        return r;
    }

    Matrix<T> gk(n, hk);
    std::vector<int> mk(hk);
    for (std::size_t a = 0; a < hk; ++a) {
        mk[a] = m[keep[a]];
        for (std::size_t k = 0; k < n; ++k) gk(k, a) = gamma(k, keep[a]);
    }

    const std::vector<T> xik = cache_xi_iter(gk, mk);
    std::vector<T> xi(h, zero);  // dropped lists keep xi = 0
    for (std::size_t a = 0; a < hk; ++a) xi[keep[a]] = xik[a];

    std::vector<T> S(n, zero);
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t l = 0; l < hk; ++l) S[k] += gk(k, l) * xik[l];

    T phi = zero;
    for (std::size_t k = 0; k < n; ++k) phi += log(one + S[k]);
    for (std::size_t l = 0; l < hk; ++l)
        phi -= log(xik[l]) * num_traits<T>::from_int(static_cast<long>(mk[l]));

    Matrix<T> C(hk, hk, zero);
    for (std::size_t j = 0; j < hk; ++j) {
        for (std::size_t l = 0; l < hk; ++l) {
            T C1 = zero, C2 = zero;
            for (std::size_t k = 0; k < n; ++k) {
                C1 += gk(k, j) / (one + S[k]);
                C2 += gk(k, j) * gk(k, l) / ((one + S[k]) * (one + S[k]));
            }
            C(j, l) = (j == l ? C1 : zero) - xik[j] * C2;
        }
    }

    const T detC = detail::cache_det(C);
    if (detC == zero) throw NumericError("cache_spm: the saddle-point Hessian is singular");

    const T pi = acos(num_traits<T>::from_int(-1));
    const T s2pi = sqrt(two * pi);

    // log(prod_l m_l!), computed from the exact factorial rather than lgamma.
    T lfact = zero;
    for (std::size_t l = 0; l < hk; ++l)
        lfact += log(num_factorial<T>(static_cast<unsigned>(mk[l])));

    // real() branch: modulus inside every logarithm, as MATLAB's real(lZ) does.
    T lsqrtxi = zero;
    for (std::size_t l = 0; l < hk; ++l) lsqrtxi += half * log(num_abs(xik[l]));

    r.lZ = -num_traits<T>::from_int(static_cast<long>(hk)) * log(s2pi) + phi + lfact - lsqrtxi -
           half * log(num_abs(detC));

    T prodfact = one, prodsqrtxi = one;
    for (std::size_t l = 0; l < hk; ++l) {
        prodfact *= num_factorial<T>(static_cast<unsigned>(mk[l]));
        prodsqrtxi *= sqrt(num_abs(xik[l]));
    }
    r.Z = exp(phi) * num_pow_int(one / s2pi, static_cast<unsigned>(hk)) * prodfact / prodsqrtxi /
          sqrt(num_abs(detC));
    r.xi = xi;
    return r;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_SPM_H
