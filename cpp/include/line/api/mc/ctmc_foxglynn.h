/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_FOXGLYNN_H
#define LINE_API_MC_CTMC_FOXGLYNN_H

/**
 * Transient distribution of a CTMC by uniformization with Fox-Glynn Poisson
 * weights.
 *
 * Templated port of matlab/src/api/mc/ctmc_foxglynn.m and
 * jar/src/main/java/jline/api/mc/Ctmc_foxglynn.java. The mixing distribution
 * Poisson(q t) is truncated to a window [L, R] carrying all but tol of its
 * mass, and the weights on that window are built by the two-sided recursion
 *   w(k-1) = w(k) k / lambda,   w(k+1) = w(k) lambda / (k+1)
 * anchored at the mode with w(mode) = 1 and normalized at the end, so neither
 * exp(-lambda) nor lambda^k / k! is ever evaluated and the method is free of
 * the overflow and underflow that force ctmc_uniformization to split long
 * horizons into segments.
 *
 * The window estimate is the Fox-Glynn (1988) one with Jansen's (2011)
 * correction factor 1/(1 - exp(-(2/9) s)) on the right tail, and it is then
 * certified in every regime by tightening or growing R (and L) against the
 * Chernoff exponent lambda h(k/lambda), h(u) = u log u - u + 1. The estimate
 * itself is asymptotic and valid only for lambda >= 25; the certification is
 * what makes the result correct below that.
 *
 * GATED ON TRANSCENDENTAL ARITHMETIC: the truncation window is defined by a
 * logarithmic tail bound and the whole construction is an approximation of
 * exp(Qt) controlled by tol, so there is nothing exact to preserve. The window
 * is located in double precision, exactly as the reference does -- it is an
 * integer pair, and computing it in extended precision would move it by
 * nothing -- while the weights and the vector-matrix recursion run in T, which
 * is what a high-precision instantiation buys.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_uniformization.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct FoxGlynnResult {
    std::vector<T> pi;  ///< distribution at time t
    long left;          ///< left truncation point
    long right;         ///< right truncation point
    std::vector<T> w;   ///< normalized Poisson weights on [left, right]
};

namespace detail {

/** pi, spelled out so the header does not depend on the POSIX M_PI extension. */
constexpr double FOXGLYNN_PI = 3.14159265358979323846;

/** Chernoff exponent lambda h(k/lambda); exp(-e) dominates the Poisson tail. */
inline double foxglynn_chernoff(double lambda, double k) {
    if (k <= 0.0) return lambda;
    return lambda - k + k * std::log(k / lambda);
}

/** Right truncation point R with P{X > R} <= tol/2. */
inline long foxglynn_right(double lambda, double tol) {
    const double target = std::log(2.0 / tol);
    const double m = std::floor(lambda);
    double r = m;
    if (lambda >= 25.0) {
        const double a = (1.0 + 1.0 / lambda) * std::exp(1.0 / 16.0) * std::sqrt(2.0);
        const double spread = std::sqrt(2.0 * lambda);
        for (int k = 1; k <= 64; ++k) {
            const double shift = k * spread + 1.5;
            const double d = 1.0 / (1.0 - std::exp(-(2.0 / 9.0) * shift));
            const double bound = a * d * std::exp(-0.5 * k * k) / (k * std::sqrt(2.0 * FOXGLYNN_PI));
            if (bound <= 0.5 * tol) {
                r = m + std::ceil(shift);
                break;
            }
        }
    }
    while (r > m && foxglynn_chernoff(lambda, r) >= target) r -= 1.0;
    while (foxglynn_chernoff(lambda, r + 1.0) < target) r += 1.0;
    return static_cast<long>(r);
}

/** Left truncation point L with P{X < L} <= tol/2, zero when none is admissible. */
inline long foxglynn_left(double lambda, double tol) {
    const double target = std::log(2.0 / tol);
    const double m = std::floor(lambda);
    if (foxglynn_chernoff(lambda, 0.0) < target) return 0;
    double l = 0.0;
    if (lambda >= 25.0) {
        const double b = (1.0 + 1.0 / lambda) * std::exp(1.0 / (8.0 * lambda));
        const double spread = std::sqrt(lambda);
        for (int k = 1; k <= 64; ++k) {
            const double bound = b * std::exp(-0.5 * k * k) / (k * std::sqrt(2.0 * FOXGLYNN_PI));
            if (bound <= 0.5 * tol) {
                l = m - std::floor(k * spread + 1.5);
                break;
            }
        }
        if (l < 0.0) l = 0.0;
    }
    while (l > 0.0 && foxglynn_chernoff(lambda, l - 1.0) < target) l -= 1.0;
    while (l < m && foxglynn_chernoff(lambda, l) >= target) l += 1.0;
    return static_cast<long>(l);
}

/**
 * Poisson weights on [left, right], by the two-sided recursion. With normalize
 * set the sum is accumulated in increasing order of magnitude, as the reference
 * does with sort(w), and the window is rescaled to one, so the truncated tails
 * are redistributed over it. Cleared, the anchor is instead scaled by the true
 * mode probability, evaluated once in double through a log-gamma, so the
 * returned values are the Poisson probabilities themselves and 1 - sum(w) is
 * the discarded tail rather than being absorbed; ctmc_fau needs that, its error
 * being reported as missing mass rather than as a bound.
 */
template <class T>
std::vector<T> foxglynn_poisson(const T& lambda, long left, long right, double lambdaDouble,
                                bool normalize = true) {
    const std::size_t len = static_cast<std::size_t>(right - left + 1);
    std::vector<T> w(len, num_traits<T>::from_int(0));
    long m = static_cast<long>(std::floor(lambdaDouble));
    if (m < left) m = left;
    if (m > right) m = right;
    w[static_cast<std::size_t>(m - left)] = num_traits<T>::from_int(1);
    for (long k = m; k >= left + 1; --k)
        w[static_cast<std::size_t>(k - 1 - left)] =
            w[static_cast<std::size_t>(k - left)] * num_traits<T>::from_int(k) / lambda;
    for (long k = m; k <= right - 1; ++k)
        w[static_cast<std::size_t>(k + 1 - left)] =
            w[static_cast<std::size_t>(k - left)] * lambda / num_traits<T>::from_int(k + 1);
    if (!normalize) {
        const double logMode = -lambdaDouble + static_cast<double>(m) * std::log(lambdaDouble) -
                               std::lgamma(static_cast<double>(m) + 1.0);
        const T scale = num_traits<T>::from_double(std::exp(logMode));
        for (T& v : w) v *= scale;
        return w;
    }
    std::vector<T> sorted = w;
    std::sort(sorted.begin(), sorted.end());
    T s = num_traits<T>::from_int(0);
    for (const T& v : sorted) s += v;
    if (s == num_traits<T>::from_int(0)) throw NumericError("ctmc_foxglynn: Poisson weights vanish");
    for (T& v : w) v /= s;
    return w;
}

}  // namespace detail

/**
 * @param pi0 initial distribution (row vector)
 * @param Q   generator
 * @param t   time horizon
 * @param tol total Poisson tail mass discarded (MATLAB default 1e-12)
 * @param maxiter cap on the right truncation point; <= 0 leaves it uncapped
 */
template <class T>
FoxGlynnResult<T> ctmc_foxglynn(const std::vector<T>& pi0, const Matrix<T>& Q, const T& t,
                                double tol = 1e-12, long maxiter = -1) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_foxglynn requires transcendental arithmetic: the truncation window is "
                  "defined by a logarithmic Poisson tail bound, and the result is an "
                  "approximation of exp(Qt) controlled by tol rather than an exact quantity");
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_foxglynn: generator is not square");
    if (pi0.size() != n) throw InputError("ctmc_foxglynn: pi0 has the wrong length");
    if (tol <= 0.0) tol = 1e-12;

    // q = 1.1 max |q_ii|, the rate the reference uses.
    T qmax = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        const T a = num_abs(T(Q(i, i)));
        if (a > qmax) qmax = a;
    }
    const T q = qmax * num_traits<T>::from_rational(11, 10);
    const T lambda = q * t;
    const double lambdaDouble = num_traits<T>::to_double(lambda);
    if (!(num_traits<T>::to_double(q) > 0.0) || !(lambdaDouble > 0.0)) {
        FoxGlynnResult<T> r;
        r.pi = pi0;
        r.left = 0;
        r.right = 0;
        r.w.assign(1, num_traits<T>::from_int(1));
        return r;
    }

    long left = detail::foxglynn_left(lambdaDouble, tol);
    long right = detail::foxglynn_right(lambdaDouble, tol);
    if (maxiter > 0 && right > maxiter) {
        right = maxiter;
        left = std::min(left, right);
    }

    FoxGlynnResult<T> r;
    r.left = left;
    r.right = right;
    r.w = detail::foxglynn_poisson(lambda, left, right, lambdaDouble);

    const Matrix<T> Qs = detail::uniformized_matrix(Q, q);
    r.pi.assign(n, num_traits<T>::from_int(0));
    std::vector<T> P = pi0;
    for (long k = 0; k <= right; ++k) {
        if (k >= left) {
            const T& wk = r.w[static_cast<std::size_t>(k - left)];
            for (std::size_t i = 0; i < n; ++i) r.pi[i] += wk * P[i];
        }
        if (k < right) P = detail::vecmat(P, Qs);
    }
    return r;
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_FOXGLYNN_H
