/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_PNT_H
#define LINE_API_MAM_MAP_PNT_H

/**
 * Counting probabilities of a MAP: P_n(t), the matrix whose (i,j) entry is the
 * probability of n arrivals in [0,t) ending in phase j, given phase i at 0.
 *
 * Port of matlab/lib/kpctoolbox/map/map_pntiter.m (and its map_pntbisect
 * helper), mirrored by jline.api.mam.Map_pntiter and the native Python
 * `map_pntiter`.
 *
 * The method is uniformization. With tau = max_i |D0(i,i)|, K = D0/tau + I and
 * K1 = D1/tau, the number of uniformization steps in [0,t) is Poisson(tau t),
 * and conditioning on it gives
 *
 *   V(0,0) = I,  V(0,k) = V(0,k-1) K,
 *   V(n,k) = V(n,k-1) K + V(n-1,k-1) K1,
 *   P_n(t) = sum_{k=0..N} w_k V(n,k),   w_k = e^{-tau t} (tau t)^k / k!,
 *
 * truncated at the N for which the Poisson tail falls below machine epsilon.
 * `map_pntiter` evaluates this at t/2^M and then squares M times through the
 * discrete convolution P_n <- sum_{j=0..n} P_j P_{n-j}, which is exact because
 * the counts over disjoint intervals add and the phase at the split is summed
 * over by the matrix product.
 *
 * REFERENCE DEFECT, NOT REPRODUCED (found 2026-08-01, verified in MATLAB
 * R2025a). `map_pntbisect.m` weights V(n,k) by w_n instead of w_k -- the
 * Poisson weight of the ARRIVAL count rather than of the uniformization step
 * count -- and never propagates V(0,k), leaving it zero for every k >= 1. Both
 * are wrong and neither is visible on a Poisson process, because there K = 0
 * collapses V(n,k) to delta(n,k) and the two weights coincide on the only
 * surviving term. On an Erlang-2 MAP of mean 1 at t = 1.3 the reference
 * returns, against the identities any counting law must satisfy:
 *
 *   P_0(t) vs exp(D0 t)                    max |diff| = 1.93e-1
 *   sum_n P_n(t) vs exp((D0+D1) t)         max |diff| = 4.97e-1
 *   sum_n n pie P_n(t) e vs lambda t       0 against 1.3
 *
 * The version here satisfies all three to round-off, and the tests assert them.
 * `map_pntiter` and `map_pntquad` have NO CALLER in any codebase, so the defect
 * is latent and changes no published result; fixing MATLAB, the JAR and Python
 * is tracked separately.
 *
 * ARITHMETIC: transcendental, for the Poisson weights.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/util/ode.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace detail {

/** Poisson weight e^{-tau t} (tau t)^r / r!, formed in logs. */
inline double pnt_weight(double tau, double t, std::size_t r) {
    const double x = tau * t;
    if (!(x > 0.0)) return r == 0 ? 1.0 : 0.0;
    double logw = -x + static_cast<double>(r) * std::log(x);
    for (std::size_t k = 2; k <= r; ++k) logw -= std::log(static_cast<double>(k));
    return std::exp(logw);
}

/** Smallest N whose Poisson tail beyond it is below machine epsilon. */
inline std::size_t pnt_truncation(double tau, double t) {
    const std::size_t kMax = 2000;
    double acc = 0.0;
    for (std::size_t N = 1; N <= kMax; ++N) {
        acc += pnt_weight(tau, t, N);
        if (1.0 - acc < 2.3e-16) return N;
    }
    return kMax;
}

}  // namespace detail

/**
 * P_0(t) .. P_na(t) by uniformization on one interval, without the squaring.
 *
 * @param m  the MAP
 * @param na highest arrival count to return
 * @param t  interval length
 */
template <class T>
std::vector<Matrix<T>> map_pntbisect(const Map<T>& m, std::size_t na, const T& t) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_pntbisect weights the terms by a Poisson law");
    const std::size_t n = m.order();
    if (n == 0 || m.D1.rows() != n) throw InputError("map_pntbisect: D0 and D1 disagree");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    double tau = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        tau = std::max(tau, -num_traits<T>::to_double(m.D0(i, i)));
    const double tv = num_traits<T>::to_double(t);
    std::vector<Matrix<T>> P(na + 1, Matrix<T>(n, n, zero));
    if (!(tau > 0.0) || !(tv > 0.0)) {  // no clock, or no time: nothing happens
        for (std::size_t i = 0; i < n; ++i) P[0](i, i) = one;
        return P;
    }

    const T taut = num_traits<T>::from_double(tau);
    Matrix<T> K(n, n, zero), K1(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            K(i, j) = T(m.D0(i, j) / taut + (i == j ? one : zero));
            K1(i, j) = T(m.D1(i, j) / taut);
        }

    const std::size_t N = detail::pnt_truncation(tau, tv);
    // V[n] is the current step's V(n,k); the previous step's is kept in Vprev.
    std::vector<Matrix<T>> V(na + 1, Matrix<T>(n, n, zero));
    for (std::size_t i = 0; i < n; ++i) V[0](i, i) = one;  // V(0,0) = I
    for (std::size_t a = 0; a <= na; ++a) {
        const T w = num_traits<T>::from_double(detail::pnt_weight(tau, tv, 0));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) P[a](i, j) = T(P[a](i, j) + w * V[a](i, j));
    }
    for (std::size_t k = 1; k <= N; ++k) {
        const std::vector<Matrix<T>> Vprev = V;
        for (std::size_t a = 0; a <= na; ++a) {
            Matrix<T> next(n, n, zero);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) {
                    T acc = zero;
                    for (std::size_t q = 0; q < n; ++q) acc += Vprev[a](i, q) * K(q, j);
                    if (a > 0)
                        for (std::size_t q = 0; q < n; ++q) acc += Vprev[a - 1](i, q) * K1(q, j);
                    next(i, j) = acc;
                }
            V[a] = next;
        }
        const T w = num_traits<T>::from_double(detail::pnt_weight(tau, tv, k));
        for (std::size_t a = 0; a <= na; ++a)
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j) P[a](i, j) = T(P[a](i, j) + w * V[a](i, j));
    }
    return P;
}

/**
 * P_0(t) .. P_na(t), evaluated on a short interval and squared up.
 *
 * @param m  the MAP
 * @param na highest arrival count to return
 * @param t  interval length
 * @param M  number of squarings; negative selects the reference's default
 *           ceil(log2(100 t / mean)), and a value below zero after that means
 *           the direct evaluation
 */
template <class T>
std::vector<Matrix<T>> map_pnt(const Map<T>& m, std::size_t na, const T& t, long M = -1) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_pnt weights the terms by a Poisson law");
    long steps = M;
    if (steps < 0) {
        const double mean = num_traits<T>::to_double(map_mean(m));
        const double tv = num_traits<T>::to_double(t);
        if (!(mean > 0.0) || !(tv > 0.0)) return map_pntbisect(m, na, t);
        steps = static_cast<long>(std::ceil(std::log2(tv * 100.0 / mean)));
        if (steps < 0) return map_pntbisect(m, na, t);
    }

    const T half = num_traits<T>::from_double(std::pow(2.0, static_cast<double>(steps)));
    std::vector<Matrix<T>> P = map_pntbisect(m, na, T(t / half));
    const std::size_t n = m.order();
    const T zero = num_traits<T>::from_int(0);
    for (long s = 0; s < steps; ++s) {
        const std::vector<Matrix<T>> Pold = P;
        for (std::size_t a = 0; a <= na; ++a) {
            Matrix<T> acc(n, n, zero);
            for (std::size_t j = 0; j <= a; ++j)
                for (std::size_t i = 0; i < n; ++i)
                    for (std::size_t c = 0; c < n; ++c) {
                        T v = zero;
                        for (std::size_t q = 0; q < n; ++q)
                            v += Pold[j](i, q) * Pold[a - j](q, c);
                        acc(i, c) = T(acc(i, c) + v);
                    }
            P[a] = acc;
        }
    }
    return P;
}

/**
 * The same counting probabilities by NUMERICAL INTEGRATION, `map_pntquad`.
 *
 * Port of matlab/lib/kpctoolbox/map/map_pntquad.m, which integrates the forward
 * equations
 *   dP_0/dt = P_0 D0,   dP_n/dt = P_n D0 + P_{n-1} D1
 * from P_0(0) = I with ode45. The reference stacks all na+1 matrices into one
 * state vector and integrates them together, which is what is done here with
 * `ode_rosenbrock4`.
 *
 * IT IS A SECOND ROUTE TO THE SAME OBJECT, not a different quantity, and that
 * is its value: `map_pnt` reaches P_n(t) by uniformization and this one by
 * quadrature, so the two agreeing is evidence neither is wrong. The
 * uniformization route is the cheaper one and is what callers should use; this
 * exists because the reference exposes it and because it is the independent
 * check the tests apply.
 *
 * @param m  the MAP
 * @param na highest arrival count to return
 * @param t  interval length
 */
template <class T>
std::vector<Matrix<T>> map_pntquad(const Map<T>& m, std::size_t na, const T& t) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_pntquad integrates the forward equations");
    const std::size_t K = m.order();
    if (K == 0 || m.D1.rows() != K) throw InputError("map_pntquad: D0 and D1 disagree");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    const std::size_t blk = K * K, dim = (na + 1) * blk;
    std::vector<T> y0(dim, zero);
    for (std::size_t i = 0; i < K; ++i) y0[i * K + i] = one;  // P_0(0) = I

    auto rhs = [&](const T&, const std::vector<T>& y) {
        std::vector<T> dy(dim, zero);
        for (std::size_t n = 0; n <= na; ++n)
            for (std::size_t i = 0; i < K; ++i)
                for (std::size_t j = 0; j < K; ++j) {
                    T acc = zero;
                    for (std::size_t q = 0; q < K; ++q) acc += y[n * blk + i * K + q] * m.D0(q, j);
                    if (n > 0)
                        for (std::size_t q = 0; q < K; ++q)
                            acc += y[(n - 1) * blk + i * K + q] * m.D1(q, j);
                    dy[n * blk + i * K + j] = acc;
                }
        return dy;
    };

    const std::vector<T> yt = ode_rosenbrock4_endpoint(rhs, zero, t, y0);
    std::vector<Matrix<T>> P(na + 1, Matrix<T>(K, K, zero));
    for (std::size_t n = 0; n <= na; ++n)
        for (std::size_t i = 0; i < K; ++i)
            for (std::size_t j = 0; j < K; ++j) P[n](i, j) = yt[n * blk + i * K + j];
    return P;
}

/** The reference's entry point: only the highest count is returned. */
template <class T>
Matrix<T> map_pntiter(const Map<T>& m, std::size_t na, const T& t, long M = -1) {
    return map_pnt(m, na, t, M)[na];
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_PNT_H
