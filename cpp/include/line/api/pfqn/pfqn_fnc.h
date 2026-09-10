/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_FNC_H
#define LINE_API_PFQN_PFQN_FNC_H

/**
 * Load-dependent rates of the functional server f(n) = n + c.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_fnc.m (Casale, "On Single-Class
 * Load-Dependent Normalizing Constant Equations", QEST 2006). Given the
 * balance-function increments alpha(i,n) of an existing station, the rates
 * mu(i,n) of the functional server are peeled off by
 *
 *   mu(i,1) = alpha(i,1)/(1 + c_i)
 *   mu(i,n) = alphanum(n,n-1) alpha(i,1) / (prod_{k<n} mu(i,k)) / (1 - rho)
 *   rho     = sum_{k=1}^{n-1} (alphanum(n,k) - alphaden(n,k)) / prod_{j<=k} mu(i,j)
 *
 * with alphanum(n,k) = prod_{j=0}^{k-1} alpha(i,n-j) and
 * alphaden(n,k) = prod_{j=0}^{k-1} alpha(i,n-1-j).
 *
 * OFFSET SEARCH. With no c supplied the reference tries c = 0, then
 * c = -1/2, then walks c upwards in steps of 0.05 until every rate is finite
 * or c reaches 2. The port keeps that ladder, including its two documented
 * repairs: the zero-column guard (a caller that shifted a single-column mu
 * gets an M x 0 result rather than an out-of-range index, which used to break
 * SolverNC 'exact' on every closed model of total population one), and the
 * fact that c is a length-M vector rather than a scalar.
 *
 * NON-FINITE RATES. MATLAB maps NaN and any |mu| > 1e15 to Inf and then
 * saturates each row from its first Inf onwards. The port keeps both, since a
 * downstream load-dependent solver reads Inf as "this station cannot hold that
 * many jobs" and would misread a large finite rate as a physical one. The
 * 1e15 threshold is a floating-point guard and is a double constant in every
 * arithmetic, exactly as in the reference.
 *
 * ARITHMETIC. The recursion is a finite sequence of field operations, so it is
 * EXACT in rational arithmetic and is deliberately left ungated -- but the
 * automatic offset ladder tests finiteness, which is a floating-point notion;
 * with an exact type the first ladder step that produces no division by zero
 * is accepted.
 */

#include <cmath>
#include <limits>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_fnc, mirroring [mu, c]. */
template <class T>
struct FncResult {
    Matrix<T> mu;
    std::vector<T> c;
};

namespace detail {

/**
 * MATLAB's retry predicate `~all(isfinite(mu))` used as an `if` condition,
 * reproduced verbatim including its vectorization.
 *
 * REFERENCE DEFECT. `all` reduces along the first non-singleton dimension, so
 * for a multi-station mu (M > 1) it returns one flag PER COLUMN, `~` negates
 * each, and `if` on a vector is true only when EVERY element is true. The
 * retry therefore fires only when every population column contains a
 * non-finite rate, not when any single rate is non-finite. On the two-station
 * example alpha = [1 2 3; 2 2 2] MATLAB accepts c = 0 and returns
 * mu = [1 Inf Inf; 2 2 2], which a plain "any entry non-finite" reading would
 * have rejected and retried. The port reproduces MATLAB, because a caller
 * comparing against it must get the same offset vector.
 */
template <class T>
bool matlab_retry_needed(const Matrix<T>& m) {
    if (m.cols() == 0) return false;
    if (m.rows() == 1) {  // row vector: all() reduces to a scalar
        for (std::size_t j = 0; j < m.cols(); ++j)
            if (is_inf_marker(m(0, j))) return true;
        return false;
    }
    for (std::size_t j = 0; j < m.cols(); ++j) {
        bool colHasNonFinite = false;
        for (std::size_t i = 0; i < m.rows(); ++i)
            if (is_inf_marker(m(i, j))) colHasNonFinite = true;
        if (!colHasNonFinite) return false;
    }
    return true;
}

}  // namespace detail

/** Rates for a given offset vector c (the two-argument MATLAB branch). */
template <class T>
Matrix<T> pfqn_fnc_at(const Matrix<T>& alpha, const std::vector<T>& c) {
    const std::size_t M = alpha.rows(), N = alpha.cols();
    if (c.size() != M) throw InputError("pfqn_fnc: the offset vector must have one entry per station");
    if (N == 0) return Matrix<T>(M, 0);
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T inf = detail::num_inf_marker<T>();

    Matrix<T> mu(M, N, zero);
    for (std::size_t i = 0; i < M; ++i) {
        const T d0 = T(one + c[i]);
        if (d0 == zero) throw NumericError("pfqn_fnc: offset -1 makes the first rate undefined");
        mu(i, 0) = T(alpha(i, 0) / d0);
        // alphanum(n,k) = prod_{j=0}^{k-1} alpha(i, n-j), 1-based n and k.
        Matrix<T> anum(N + 1, N + 1, zero), aden(N + 1, N + 1, zero);
        for (std::size_t n = 2; n <= N; ++n) {
            anum(n, 1) = alpha(i, n - 1);
            aden(n, 1) = alpha(i, n - 2);
            for (std::size_t k = 2; k + 1 <= n; ++k) {
                anum(n, k) = T(anum(n, k - 1) * alpha(i, n - k));
                aden(n, k) = T(aden(n, k - 1) * alpha(i, n - k - 1));
            }
        }
        for (std::size_t n = 2; n <= N; ++n) {
            T rho = zero, muden = one;
            bool bad = false;
            for (std::size_t k = 1; k + 1 <= n; ++k) {
                muden *= mu(i, k - 1);
                if (muden == zero) {
                    bad = true;
                    break;
                }
                rho += T(T(anum(n, k) - aden(n, k)) / muden);
            }
            if (bad || T(one - rho) == zero) {
                mu(i, n - 1) = inf;
                continue;
            }
            T v = T(anum(n, n - 1) * alpha(i, 0) / muden);
            v = T(v / T(one - rho));
            mu(i, n - 1) = v;
        }
    }
    // MATLAB: mu(isnan) = Inf; mu(abs(mu) > 1e15) = Inf; then saturate each
    // row from its first Inf onwards.
    const double big = 1e15;
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t n = 0; n < N; ++n) {
            const double v = num_traits<T>::to_double(mu(i, n));
            if (std::isnan(v) || std::fabs(v) > big) mu(i, n) = inf;  // the reference's Inf
        }
        for (std::size_t n = 0; n < N; ++n) {
            if (detail::is_inf_marker(mu(i, n))) {
                for (std::size_t k = n; k < N; ++k) mu(i, k) = inf;
                break;
            }
        }
    }
    return mu;
}

/** Automatic offset search (the one-argument MATLAB branch). */
template <class T>
FncResult<T> pfqn_fnc(const Matrix<T>& alpha) {
    const std::size_t M = alpha.rows();
    FncResult<T> r;
    if (alpha.cols() == 0) {
        // A caller that shifted a single-column mu has no rate to build.
        r.mu = Matrix<T>(M, 0);
        r.c.assign(M, num_traits<T>::from_int(0));
        return r;
    }
    r.c.assign(M, num_traits<T>::from_int(0));
    r.mu = pfqn_fnc_at(alpha, r.c);
    if (!detail::matlab_retry_needed(r.mu)) return r;
    r.c.assign(M, num_traits<T>::from_rational(-1, 2));
    r.mu = pfqn_fnc_at(alpha, r.c);
    if (!detail::matlab_retry_needed(r.mu)) return r;
    for (int step = 1; step <= 50; ++step) {
        const T v = T(num_traits<T>::from_rational(-1, 2) +
                      num_traits<T>::from_rational(step, 20));  // dt in steps of 0.05
        r.c.assign(M, v);
        r.mu = pfqn_fnc_at(alpha, r.c);
        if (!detail::matlab_retry_needed(r.mu)) return r;
        if (num_traits<T>::to_double(v) >= 2.0) break;
    }
    return r;
}

template <class T>
FncResult<T> pfqn_fnc(const Matrix<T>& alpha, const std::vector<T>& c) {
    FncResult<T> r;
    r.c = c;
    r.mu = pfqn_fnc_at(alpha, c);
    return r;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_FNC_H
