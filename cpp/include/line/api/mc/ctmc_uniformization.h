/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_UNIFORMIZATION_H
#define LINE_API_MC_CTMC_UNIFORMIZATION_H

/**
 * Transient distribution of a CTMC by uniformization (Jensen's method), and
 * the time-averaged distribution over [0, t].
 *
 * Templated port of matlab/src/api/mc/ctmc_uniformization.m and
 * matlab/src/api/mc/ctmc_timeaverage.m.
 *
 *   pi(t) = pi0 exp(Qt) = sum_j Poisson(j; q t) pi0 P^j,  P = I + Q/q
 *
 * with q = 1.1 max |diag(Q)|. The series is truncated at the first k whose
 * Poisson tail falls below tol, and long horizons are split into segments of
 * q t <= 500 exactly as MATLAB does, because the Poisson weights underflow
 * before that bound is reached.
 *
 * Unlike the steady-state routines in ctmc_solve.h, this one is NOT exact in
 * the rational field: exp(-q t) is transcendental and the truncation itself is
 * an approximation controlled by tol. The static_assert makes that explicit at
 * compile time rather than leaving a caller to discover it at runtime. The
 * high-precision instantiation is still useful: it pushes the underflow of the
 * Poisson weights far out, which is the failure mode that forces the
 * segmentation in the first place.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct UniformizationResult {
    std::vector<T> pi;  ///< distribution at time t
    std::size_t kmax;   ///< number of Poisson terms actually used
};

template <class T>
struct TimeAverageResult {
    std::vector<T> piTimeAvg;  ///< time-averaged distribution over [0, t]
    std::vector<T> piExit;     ///< distribution at time t
    std::size_t kmax;
};

namespace detail {

/** Uniformization rate q = 1.1 max_i |Q(i,i)|, and P = I + Q/q. */
template <class T>
T uniformization_rate(const Matrix<T>& Q) {
    const std::size_t n = Q.rows();
    T qmax = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        const T a = num_abs(Q(i, i));
        if (a > qmax) qmax = a;
    }
    return qmax * num_traits<T>::from_rational(11, 10);
}

template <class T>
Matrix<T> uniformized_matrix(const Matrix<T>& Q, const T& q) {
    const std::size_t n = Q.rows();
    Matrix<T> P(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            P(i, j) = Q(i, j) / q + (i == j ? num_traits<T>::from_int(1) : num_traits<T>::from_int(0));
    return P;
}

/** Row-vector times matrix, v P. */
template <class T>
std::vector<T> vecmat(const std::vector<T>& v, const Matrix<T>& P) {
    const std::size_t n = P.rows();
    std::vector<T> r(P.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) {
        if (v[i] == num_traits<T>::from_int(0)) continue;
        for (std::size_t j = 0; j < P.cols(); ++j) r[j] += v[i] * P(i, j);
    }
    return r;
}

/** Truncation point: first k whose Poisson tail at qt is below tol. */
template <class T>
std::size_t poisson_truncation(double qt, double tol, long maxiter) {
    if (maxiter <= 0)
        maxiter = static_cast<long>(std::max(100.0, std::ceil(qt + 10.0 * std::sqrt(qt) + 20.0)));
    double s = 1.0, r = 1.0;
    std::size_t kmax = 1;
    const double e = std::exp(-qt);
    for (long iter = 0, k = 0; iter < maxiter; ++iter) {
        ++k;
        r = r * qt / static_cast<double>(k);
        s += r;
        kmax = static_cast<std::size_t>(k);
        if ((1.0 - e * s) <= tol) break;
    }
    return kmax;
}

}  // namespace detail

/**
 * @param pi0 initial distribution (row vector)
 * @param Q   generator
 * @param t   time horizon
 * @param tol Poisson-tail truncation tolerance (MATLAB default 1e-12)
 * @param maxiter iteration cap, <= 0 for the MATLAB heuristic
 */
template <class T>
UniformizationResult<T> ctmc_uniformization(const std::vector<T>& pi0, const Matrix<T>& Q, const T& t,
                                            double tol = 1e-12, long maxiter = -1) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_uniformization requires transcendental arithmetic: the Poisson weights "
                  "involve exp(-q t), which is not a rational function of the rates");
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_uniformization: generator is not square");
    if (pi0.size() != n) throw InputError("ctmc_uniformization: pi0 has the wrong length");

    const T q = detail::uniformization_rate(Q);
    if (q == num_traits<T>::from_int(0)) return {pi0, 0};

    // Long horizons are split so that q t stays inside the range where the
    // Poisson weights are representable, as in MATLAB (MAXQT = 500).
    const double qt_full = num_traits<T>::to_double(q) * num_traits<T>::to_double(t);
    const double MAXQT = 500.0;
    if (qt_full > MAXQT) {
        const long nSeg = static_cast<long>(std::ceil(qt_full / MAXQT));
        const T tSeg = t / num_traits<T>::from_int(nSeg);
        UniformizationResult<T> r{pi0, 0};
        for (long s = 0; s < nSeg; ++s) {
            UniformizationResult<T> step = ctmc_uniformization(r.pi, Q, tSeg, tol, maxiter);
            r.pi = step.pi;
            r.kmax = step.kmax;
        }
        return r;
    }

    const Matrix<T> P = detail::uniformized_matrix(Q, q);
    const T qt = q * t;
    const std::size_t kmax = detail::poisson_truncation<T>(num_traits<T>::to_double(qt), tol, maxiter);

    using std::exp;
    T ri = exp(-qt);
    std::vector<T> pi(n);
    for (std::size_t i = 0; i < n; ++i) pi[i] = pi0[i] * ri;
    std::vector<T> Pk = pi0;
    for (std::size_t j = 1; j <= kmax; ++j) {
        Pk = detail::vecmat(Pk, P);
        ri = ri * qt / num_traits<T>::from_int(static_cast<long>(j));
        for (std::size_t i = 0; i < n; ++i) pi[i] += ri * Pk[i];
    }
    return {pi, kmax};
}

/**
 * Time-averaged distribution (1/t) int_0^t pi(u) du, plus pi(t) itself.
 * Port of ctmc_timeaverage.m, which accumulates the Poisson survival weights
 * max(1 - W_j, 0) rather than integrating pi(u) numerically.
 */
template <class T>
TimeAverageResult<T> ctmc_timeaverage(const std::vector<T>& pi0, const Matrix<T>& Q, const T& t,
                                      double tol = 1e-12, long maxiter = -1) {
    static_assert(num_traits<T>::has_transcendental,
                  "ctmc_timeaverage requires transcendental arithmetic (Poisson weights)");
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_timeaverage: generator is not square");
    if (pi0.size() != n) throw InputError("ctmc_timeaverage: pi0 has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    const T q = detail::uniformization_rate(Q);
    if (q == zero) return {pi0, pi0, 0};

    const double qt_full = num_traits<T>::to_double(q) * num_traits<T>::to_double(t);
    const double MAXQT = 500.0;
    if (qt_full > MAXQT) {
        const long nSeg = static_cast<long>(std::ceil(qt_full / MAXQT));
        const T tSeg = t / num_traits<T>::from_int(nSeg);
        std::vector<T> cur = pi0, integral(n, zero);
        std::size_t kmax = 0;
        for (long s = 0; s < nSeg; ++s) {
            TimeAverageResult<T> seg = ctmc_timeaverage(cur, Q, tSeg, tol, maxiter);
            for (std::size_t i = 0; i < n; ++i) integral[i] += tSeg * seg.piTimeAvg[i];
            cur = seg.piExit;
            kmax = seg.kmax;
        }
        for (std::size_t i = 0; i < n; ++i) integral[i] /= t;
        return {integral, cur, kmax};
    }

    const Matrix<T> P = detail::uniformized_matrix(Q, q);
    const T qt = q * t;
    const std::size_t kmax = detail::poisson_truncation<T>(num_traits<T>::to_double(qt), tol, maxiter);

    using std::exp;
    T w = exp(-qt);  // Poisson pmf w_0
    T W = w;         // Poisson cdf W_0
    std::vector<T> Pk = pi0;
    std::vector<T> piExit(n), piIntSum(n);
    const T one = num_traits<T>::from_int(1);
    T tail = one - W;
    if (tail < zero) tail = zero;
    for (std::size_t i = 0; i < n; ++i) {
        piExit[i] = w * Pk[i];
        piIntSum[i] = tail * Pk[i];
    }
    for (std::size_t j = 1; j <= kmax; ++j) {
        Pk = detail::vecmat(Pk, P);
        w = w * qt / num_traits<T>::from_int(static_cast<long>(j));
        W += w;
        tail = one - W;
        if (tail < zero) tail = zero;
        for (std::size_t i = 0; i < n; ++i) {
            piExit[i] += w * Pk[i];
            piIntSum[i] += tail * Pk[i];
        }
    }
    std::vector<T> avg(n);
    for (std::size_t i = 0; i < n; ++i) avg[i] = piIntSum[i] / qt;
    return {avg, piExit, kmax};
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_UNIFORMIZATION_H
