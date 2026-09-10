/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_DMC_H
#define LINE_API_QSYS_QSYS_DMC_H

/**
 * D/M/c: deterministic interarrival times, exponential service.
 *
 * Templated port of matlab/src/api/qsys/qsys_dmc.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_dmc.java.
 *
 * The system is embedded at arrival epochs. Between two arrivals only
 * departures occur, so the sub-generator is the death-only bidiagonal matrix
 * A[m,m-1] = min(m,c) mu, A[m,m] = -min(m,c) mu, and the embedded chain is
 * X_{n+1} = exp(A s)[X_n + 1, :] with s = 1/lambda the interarrival time. The
 * stationary vector at arrival epochs is then converted to time averages by
 * integrating the state-count expectations over one interarrival cycle with
 * the trapezoid rule on quadSteps subintervals.
 *
 * The truncation level is max(200, min(2500, floor(15/(1-rho)) + 200)) unless
 * given, as in MATLAB. That default is expensive here because the exponential
 * and the cycle integration are dense; callers benchmarking against MATLAB on
 * a specific instance should pass the same explicit truncation to both.
 *
 * ARITHMETIC. The matrix exponential and the trapezoid rule are both inexact,
 * so the function is gated on transcendental arithmetic.
 *
 * At c = 1 the D/M/1 mean waiting time must agree with qsys_gm1 evaluated at
 * the root of sigma = exp(-mu(1-sigma)/lambda), which is the check the tests
 * apply.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_quadrature.h"
#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

template <class T>
struct DmcResult {
    T meanQueueLength;   ///< time-average E[N]
    T meanWaitingQueue;  ///< time-average Lq = E[(N-c)+]
    T meanWaitingTime;   ///< Wq = Lq/lambda
    T meanSojournTime;   ///< W = Wq + 1/mu
    T utilization;       ///< rho = lambda/(c mu)
};

namespace detail {

/**
 * Matrix exponential by scaling and squaring around a Taylor series.
 *
 * MATLAB uses expm, i.e. scaling and squaring around a Pade approximant. The
 * two agree to round-off: after scaling the argument to infinity-norm below
 * 1/2 the truncated Taylor series of 30 terms has a remainder below 2^-30/30!,
 * far under any working precision, so the squaring stage is what determines
 * the accuracy in both cases.
 */
template <class T>
Matrix<T> expm(const Matrix<T>& A) {
    static_assert(num_traits<T>::has_transcendental, "expm requires transcendental arithmetic");
    const std::size_t n = A.rows();
    if (A.cols() != n) throw InputError("expm: matrix is not square");
    T nrm = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < n; ++j) s += num_abs(A(i, j));
        if (s > nrm) nrm = s;
    }
    unsigned sq = 0;
    T scaled = nrm;
    const T half = num_traits<T>::from_rational(1, 2);
    const T two = num_traits<T>::from_int(2);
    while (scaled > half) {
        scaled /= two;
        ++sq;
    }
    Matrix<T> B = A;
    const T factor = num_pow_int(half, sq);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) B(i, j) *= factor;

    Matrix<T> E = eye<T>(n);
    Matrix<T> term = eye<T>(n);
    for (unsigned q = 1; q <= 30u; ++q) {
        term = matmul(term, B);
        const T inv = num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(q));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) term(i, j) *= inv;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) E(i, j) += term(i, j);
    }
    for (unsigned r = 0; r < sq; ++r) E = matmul(E, E);
    return E;
}

}  // namespace detail

/**
 * @param lambda    deterministic arrival rate, interarrival time 1/lambda
 * @param mu        exponential service rate of one server
 * @param c         number of servers, c >= 1
 * @param truncation state-space truncation; 0 selects the MATLAB default
 * @param quadSteps trapezoid steps over one interarrival cycle (MATLAB 200)
 */
template <class T>
DmcResult<T> qsys_dmc(const T& lambda, const T& mu, unsigned c, unsigned truncation,
                      unsigned quadSteps) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_dmc requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (lambda <= zero) throw InputError("qsys_dmc: arrival rate must be positive");
    if (mu <= zero) throw InputError("qsys_dmc: service rate must be positive");
    if (c < 1) throw InputError("qsys_dmc: number of servers must be at least 1");
    if (quadSteps < 1) throw InputError("qsys_dmc: quadSteps must be at least 1");
    const T ct = num_traits<T>::from_int(static_cast<long>(c));
    const T rho = lambda / (ct * mu);
    if (rho >= one) throw InputError("qsys_dmc: load rho must be strictly less than 1");

    const T s = one / lambda;
    unsigned nMax;
    if (truncation > 0) {
        nMax = truncation;
    } else {
        const double gap = num_traits<T>::to_double(T(one - rho));
        const long guess = static_cast<long>(std::floor(15.0 / gap)) + 200;
        nMax = static_cast<unsigned>(std::max(200L, std::min(2500L, guess)));
    }
    const std::size_t n = nMax + 1;

    Matrix<T> A(n, n, zero);
    for (std::size_t m = 0; m < n; ++m) {
        const unsigned busy = std::min<unsigned>(static_cast<unsigned>(m), c);
        const T rate = num_traits<T>::from_int(static_cast<long>(busy)) * mu;
        A(m, m) = -rate;
        if (m > 0) A(m, m - 1) = rate;
    }

    Matrix<T> As = A, Adt = A;
    const T dt = s / num_traits<T>::from_int(static_cast<long>(quadSteps));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            As(i, j) *= s;
            Adt(i, j) *= dt;
        }
    const Matrix<T> expAs = detail::expm(As);
    const Matrix<T> expAdt = detail::expm(Adt);

    // Row r represents "number in system before the arrival = r"; after the
    // arrival the state is min(r+1, n-1).
    std::vector<std::size_t> yIdx(n);
    for (std::size_t r = 0; r < n; ++r) yIdx[r] = std::min(r + 1, n - 1);

    // Stationary at arrival epochs: (P' - I) pi = 0 with the last row replaced
    // by the normalization.
    Matrix<T> M(n, n, zero);
    for (std::size_t r = 0; r + 1 < n; ++r)
        for (std::size_t col = 0; col < n; ++col) M(r, col) = expAs(yIdx[col], r);
    for (std::size_t r = 0; r + 1 < n; ++r) M(r, r) -= one;
    for (std::size_t col = 0; col < n; ++col) M(n - 1, col) = one;
    std::vector<T> b(n, zero);
    b[n - 1] = one;
    const std::vector<T> piArr = line::solve(M, b);

    std::vector<T> wLq(n), wN(n);
    for (std::size_t m = 0; m < n; ++m) {
        wN[m] = num_traits<T>::from_int(static_cast<long>(m));
        wLq[m] = m > c ? num_traits<T>::from_int(static_cast<long>(m - c)) : zero;
    }

    std::vector<T> ts(quadSteps + 1);
    for (unsigned q = 0; q <= quadSteps; ++q)
        ts[q] = s * num_traits<T>::from_int(static_cast<long>(q)) /
                num_traits<T>::from_int(static_cast<long>(quadSteps));
    std::vector<std::vector<T>> LqAtT(quadSteps + 1), NAtT(quadSteps + 1);
    Matrix<T> expAt = eye<T>(n);
    for (unsigned q = 0; q <= quadSteps; ++q) {
        if (q > 0) expAt = matmul(expAt, expAdt);
        LqAtT[q] = mulvec(expAt, wLq);
        NAtT[q] = mulvec(expAt, wN);
    }
    std::vector<T> LqInt(n, zero), NInt(n, zero);
    std::vector<T> colLq(quadSteps + 1), colN(quadSteps + 1);
    for (std::size_t r = 0; r < n; ++r) {
        for (unsigned q = 0; q <= quadSteps; ++q) {
            colLq[q] = LqAtT[q][r];
            colN[q] = NAtT[q][r];
        }
        LqInt[r] = detail::num_trapz(ts, colLq) / s;
        NInt[r] = detail::num_trapz(ts, colN) / s;
    }

    T LqTime = zero, NTime = zero;
    for (std::size_t r = 0; r < n; ++r) {
        LqTime += piArr[r] * LqInt[yIdx[r]];
        NTime += piArr[r] * NInt[yIdx[r]];
    }

    DmcResult<T> res;
    res.meanQueueLength = NTime;
    res.meanWaitingQueue = LqTime;
    res.meanWaitingTime = LqTime / lambda;
    res.meanSojournTime = res.meanWaitingTime + one / mu;
    res.utilization = rho;
    return res;
}

/** qsys_dmc with the MATLAB defaults, automatic truncation and 200 steps. */
template <class T>
DmcResult<T> qsys_dmc(const T& lambda, const T& mu, unsigned c) {
    return qsys_dmc(lambda, mu, c, 0u, 200u);
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_DMC_H
