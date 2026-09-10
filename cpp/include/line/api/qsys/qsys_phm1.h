/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_PHM1_H
#define LINE_API_QSYS_QSYS_PHM1_H

/**
 * Exact PH/M/1, the GI/M/1 queue with phase-type interarrival times.
 *
 * Templated port of matlab/src/api/qsys/qsys_phm1.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_phm1.java.
 *
 * With psi_A(s) = alpha (sI - T)^-1 (-T e) the interarrival LST, sigma is the
 * root in (0,1) of the GI/M/1 equation
 *
 *   sigma = psi_A(mu(1 - sigma)),
 *
 * and then L = rho/(1-sigma), Lq = rho sigma/(1-sigma), Wq = Lq/lambda,
 * W = Wq + 1/mu. Since psi_A is completely monotone and psi_A(0) = 1, the
 * function f(sigma) = sigma - psi_A(mu(1-sigma)) is negative just above zero
 * and positive just below one whenever rho < 1, so the bracket [0,1] always
 * contains the root and bisection is unconditional. MATLAB brackets on
 * [1e-12, 1-1e-12] with fzero and falls back to fixed-point iteration; the
 * port bisects on the same bracket, which reaches the same root -- f is
 * strictly increasing there -- without the fallback.
 *
 * ARITHMETIC. sigma is defined by a transcendental equation and reached by a
 * tolerance-driven iteration, so the function is gated.
 *
 * At k = 1 with T = [-lambda] the arrival process is Poisson, psi_A is
 * lambda/(s+lambda), the root is sigma = rho and every metric collapses onto
 * the M/M/1 values. That identity is the sharpest available check on the port.
 */

#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

template <class T>
struct PhM1Result {
    T meanQueueLength;   ///< L
    T meanWaitingQueue;  ///< Lq
    T meanWaitingTime;   ///< Wq
    T meanSojournTime;   ///< W
    T utilization;       ///< rho = lambda/mu
    T sigma;             ///< GI/M/1 root in (0,1)
};

namespace detail {

/** psi_A(s) = alpha (sI - T)^-1 (-T e), the interarrival LST. */
template <class T>
T ph_lst(const std::vector<T>& alpha, const Matrix<T>& Tm, const T& s) {
    const std::size_t k = Tm.rows();
    Matrix<T> A(k, k);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) A(i, j) = (i == j ? s : num_traits<T>::from_int(0)) - Tm(i, j);
    std::vector<T> t_vec(k, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) t_vec[i] -= Tm(i, j);
    const std::vector<T> y = line::solve(A, t_vec);
    T out = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < k; ++i) out += alpha[i] * y[i];
    return out;
}

}  // namespace detail

/**
 * @param alpha PH entry probability vector, length k
 * @param Tm    PH sub-generator, k x k
 * @param mu    exponential service rate
 * @param tol   width of the bisection bracket at which to stop (MATLAB's
 *              fzero stops at the double round-off level, so 1e-16 is the
 *              matching default)
 */
template <class T>
PhM1Result<T> qsys_phm1(const std::vector<T>& alpha, const Matrix<T>& Tm, const T& mu,
                        const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_phm1 requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (mu <= zero) throw InputError("qsys_phm1: service rate mu must be positive");
    const std::size_t k = Tm.rows();
    if (Tm.cols() != k) throw InputError("qsys_phm1: T must be square");
    if (alpha.size() != k) throw InputError("qsys_phm1: alpha length must match T dimension");

    // Mean interarrival time = alpha (-T)^-1 e.
    Matrix<T> negT(k, k);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) negT(i, j) = -Tm(i, j);
    const std::vector<T> m = line::solve(negT, ones<T>(k));
    T mean_ia = zero;
    for (std::size_t i = 0; i < k; ++i) mean_ia += alpha[i] * m[i];
    if (mean_ia <= zero) throw InputError("qsys_phm1: non-positive mean interarrival time");
    const T lambda = one / mean_ia;
    const T rho = lambda / mu;
    if (rho >= one) throw InputError("qsys_phm1: load rho must be strictly less than 1");

    T a = T(num_traits<T>::from_double(1e-12));
    T b = one - T(num_traits<T>::from_double(1e-12));
    for (unsigned it = 0; it < 4000u; ++it) {
        if (b - a <= tol) break;
        const T mid = (a + b) / num_traits<T>::from_int(2);
        const T f = mid - detail::ph_lst(alpha, Tm, T(mu * (one - mid)));
        if (f < zero)
            a = mid;
        else
            b = mid;
    }
    const T sigma = (a + b) / num_traits<T>::from_int(2);

    PhM1Result<T> r;
    r.sigma = sigma;
    r.utilization = rho;
    r.meanQueueLength = rho / (one - sigma);
    r.meanWaitingQueue = rho * sigma / (one - sigma);
    r.meanWaitingTime = r.meanWaitingQueue / lambda;
    r.meanSojournTime = r.meanWaitingTime + one / mu;
    return r;
}

/** qsys_phm1 with the fzero-equivalent default bracket tolerance. */
template <class T>
PhM1Result<T> qsys_phm1(const std::vector<T>& alpha, const Matrix<T>& Tm, const T& mu) {
    return qsys_phm1(alpha, Tm, mu, T(num_traits<T>::from_double(1e-16)));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_PHM1_H
