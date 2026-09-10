/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_DNC_H
#define LINE_API_PFQN_DNC_H

/**
 * Distinct-load Normalizing Constant (DNC) at a nonintegral population.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_dnc.m.
 *
 * Normalizing constant and throughput of a single-class closed product-form
 * network at a REAL-VALUED population, by partial-fraction inversion of the
 * network generating function (Dowdy and Gordon 1984). With distinct loads
 * x_1..x_G of multiplicities m_1..m_G the generating function
 * prod_g (1 - x_g u)^{-m_g} expands as
 *
 *   G(n) = sum_g sum_{j = 1..m_g} A_gj C(n+j-1, j-1) x_g^n,
 *
 * every term of which is an analytic function of n. Evaluating it at a real n
 * therefore interpolates the integral normalizing constants exactly (it
 * reproduces them at every integer) and gives a smooth throughput curve
 * X(N) = G(N-1)/G(N) through the integral points, rather than the rounding or
 * linear interpolation the paper compares against. For all-distinct loads the
 * coefficients have the closed form A_g = prod_{l != g} x_g / (x_g - x_l), used
 * directly; with repeated loads they are recovered from the M integral
 * constants G(0..M-1), which determine them uniquely.
 *
 * Only the queueing part admits this continuation: the delay sequence Z^n / n!
 * is entire and has no partial-fraction expansion, so a think time is rejected
 * here. Use pfqn_nintmva for nonintegral populations with a delay.
 *
 * Reference: L. W. Dowdy, K. D. Gordon, "Algorithms for Nonintegral Degrees of
 * Multiprogramming in Closed Queuing Networks", Performance Evaluation
 * 4(1):19-28, 1984.
 *
 * Arithmetic: TRANSCENDENTAL. The partial-fraction series is evaluated in the
 * log domain through gamma functions, so the exact backend is refused at
 * compile time.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Normalizing constant and throughput at a real-valued population. */
template <class T>
struct DncResult {
    T X;   ///< throughput G(N-1)/G(N); NaN where the continuation does not apply
    T G;   ///< normalizing constant at population N
    T lG;  ///< log of the normalizing constant
};

namespace detail {

/**
 * Partial-fraction series evaluated at a real population n. The continuation is
 * analytic for n > -1; below that the binomial factor changes sign and the
 * log-domain evaluation would lose it, so it is not extended there.
 */
template <class T>
T dnc_eval(const T& n, const std::vector<T>& A, const std::vector<T>& u,
           const std::vector<std::size_t>& node, const std::vector<long>& j) {
    using std::exp;
    using std::log;
    const T minus_one = num_traits<T>::from_int(-1);
    if (n <= minus_one) return std::numeric_limits<T>::quiet_NaN();
    const T one = num_traits<T>::from_int(1);
    T g = num_traits<T>::from_int(0);
    for (std::size_t k = 0; k < A.size(); ++k) {
        const T jT = num_traits<T>::from_int(j[k]);
        const T e = T(num_lgamma<T>(T(n + jT)) - num_lgamma<T>(jT) - num_lgamma<T>(T(n + one)) +
                      n * log(u[node[k]]));
        g += A[k] * exp(e);
    }
    return g;
}

}  // namespace detail

/**
 * @param L (M) service demands of the queueing stations
 * @param N population, real and nonnegative (may be fractional)
 */
template <class T>
DncResult<T> pfqn_dnc(const std::vector<T>& L, const T& N) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_dnc needs logarithms and is not available in exact arithmetic");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<T> Lp;
    for (std::size_t i = 0; i < L.size(); ++i)
        if (L[i] > zero) Lp.push_back(L[i]);
    if (Lp.empty())
        throw InputError("pfqn_dnc requires at least one station with positive demand");
    if (N < zero) throw InputError("pfqn_dnc requires a nonnegative population");

    const std::size_t M = Lp.size();
    T xmax = Lp[0];
    for (std::size_t i = 1; i < M; ++i)
        if (Lp[i] > xmax) xmax = Lp[i];
    std::vector<T> y(M);
    for (std::size_t i = 0; i < M; ++i) y[i] = T(Lp[i] / xmax);

    // Distinct loads and multiplicities, merged under a relative tolerance so
    // that numerically coincident loads take the multiplicity branch rather
    // than a near-singular partial-fraction denominator.
    std::vector<T> ys = y;
    std::sort(ys.begin(), ys.end(), [](const T& a, const T& b) { return b < a; });
    const T near = num_traits<T>::from_double(1.0 - 1e-9);
    std::vector<T> u;
    std::vector<long> mult;
    u.push_back(ys[0]);
    mult.push_back(1);
    for (std::size_t i = 1; i < M; ++i) {
        if (ys[i] > T(u.back() * near)) {
            ++mult.back();
        } else {
            u.push_back(ys[i]);
            mult.push_back(1);
        }
    }
    const std::size_t Gd = u.size();

    std::vector<T> A;
    std::vector<std::size_t> node;
    std::vector<long> j;
    bool all_simple = true;
    for (std::size_t g = 0; g < Gd; ++g)
        if (mult[g] != 1) all_simple = false;

    if (all_simple) {
        A.assign(Gd, one);
        node.resize(Gd);
        j.assign(Gd, 1);
        for (std::size_t g = 0; g < Gd; ++g) {
            node[g] = g;
            T prod = one;
            for (std::size_t l = 0; l < Gd; ++l) {
                if (l == g) continue;
                const T den = T(u[g] - u[l]);
                if (den == zero) throw NumericError("pfqn_dnc: coincident loads in the simple branch");
                prod *= T(u[g] / den);
            }
            A[g] = prod;
        }
    } else {
        // Repeated loads: recover the coefficients from G(0..M-1), computed by
        // convolution on the scaled loads (bounded by construction, max u = 1).
        std::vector<T> gint(M, zero);
        gint[0] = one;
        std::size_t len = 1;
        for (std::size_t i = 0; i < M; ++i) {
            std::vector<T> gi(M);
            gi[0] = one;
            for (std::size_t k = 1; k < M; ++k) gi[k] = T(gi[k - 1] * y[i]);
            std::vector<T> out(M, zero);
            const std::size_t newlen = std::min(M, len + M - 1);
            for (std::size_t a = 0; a < len; ++a)
                for (std::size_t b = 0; b + a < M; ++b) out[a + b] += gint[a] * gi[b];
            gint = out;
            len = newlen;
        }
        node.resize(M);
        j.resize(M);
        std::size_t c = 0;
        for (std::size_t g = 0; g < Gd; ++g)
            for (long jj = 1; jj <= mult[g]; ++jj) {
                node[c] = g;
                j[c] = jj;
                ++c;
            }
        Matrix<T> F(M, M);
        for (std::size_t n = 0; n < M; ++n) {
            const T nT = num_traits<T>::from_int(static_cast<long>(n));
            for (std::size_t k = 0; k < M; ++k) {
                const T jT = num_traits<T>::from_int(j[k]);
                F(n, k) = exp(T(detail::num_lgamma<T>(T(nT + jT)) - detail::num_lgamma<T>(jT) -
                                detail::num_lgamma<T>(T(nT + one)) + nT * log(u[node[k]])));
            }
        }
        A = solve(F, gint);
    }

    DncResult<T> res;
    const T GN = detail::dnc_eval<T>(N, A, u, node, j);
    const T GN1 = detail::dnc_eval<T>(T(N - one), A, u, node, j);
    res.lG = T(log(GN) + N * log(xmax));
    res.G = exp(res.lG);
    const bool gn1_nan = !(GN1 == GN1);
    if (N <= zero || GN <= zero || gn1_nan) {
        res.X = std::numeric_limits<T>::quiet_NaN();
    } else {
        res.X = T(T(GN1 / GN) / xmax);
    }
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_DNC_H
