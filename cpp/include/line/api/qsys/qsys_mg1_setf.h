/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MG1_SETF_H
#define LINE_API_QSYS_QSYS_MG1_SETF_H

/**
 * M/G/1 under SETF (shortest elapsed time first), the non-preemptive
 * counterpart of FB/LAS.
 *
 * Templated port of matlab/src/api/qsys/qsys_mg1_setf.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mg1_setf.java.
 *
 * Each class is represented by the single job size x_k = 1/mu_k and evaluated
 * at that size:
 *
 *   rho_x   = sum_i lambda_i int_0^x Fbar_i(t) dt
 *   num(x)  = sum_i lambda_i int_0^x t Fbar_i(t) dt
 *   E[R]    = sum_i (lambda_i/lambda) (1+cs_i^2)/(2 mu_i)
 *   W_k     = num(x)/(1-rho_x)^2 + x/(1-rho_x) + E[R]/(1-rho_x)
 *
 * the last term being the non-preemptive penalty that separates SETF from FB.
 * For an exponential class the two truncated integrals are closed forms; for
 * any other class MATLAB substitutes the bounded surrogates min(x, 1/mu_i) and
 * min(x^2/2, 1/mu_i^2), which is an approximation and is reproduced verbatim
 * here rather than improved on.
 *
 * Note the branch test: MATLAB compares cs(i) == 1 exactly in this file, not
 * within a tolerance as its FB and PSJF siblings do, so a cs of 1 - 1e-12
 * takes the surrogate branch. That asymmetry is part of the reference
 * behaviour and is preserved.
 *
 * ARITHMETIC. exp appears in the exponential branch, so the function is gated
 * on transcendental arithmetic.
 *
 * The returned rhohat is Q/(1+Q) with Q = sum_k lambda_k W_k, the qsys family
 * convention; the utilization sum_k rho_k is not part of the return value
 * because MATLAB overwrites it.
 */

#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

template <class T>
struct Mg1DisciplineResult {
    std::vector<T> W;  ///< per-class mean response time
    T rhohat;          ///< Q/(1+Q) with Q = sum_k lambda_k W_k
};

namespace detail {

/** Shared argument validation for the M/G/1 discipline family. */
template <class T>
T mg1_discipline_check(const char* fn, const std::vector<T>& lambda, const std::vector<T>& mu,
                       const std::vector<T>& cs) {
    const std::size_t K = lambda.size();
    if (mu.size() != K || cs.size() != K)
        throw InputError(std::string(fn) + ": lambda, mu and cs must have the same length");
    if (K == 0) throw InputError(std::string(fn) + ": at least one class is required");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    T rho = zero;
    for (std::size_t i = 0; i < K; ++i) {
        if (lambda[i] <= zero || mu[i] <= zero)
            throw InputError(std::string(fn) + ": lambda and mu must be positive");
        if (cs[i] < zero) throw InputError(std::string(fn) + ": cs must be non-negative");
        rho += lambda[i] / mu[i];
    }
    if (rho >= one) throw InputError(std::string(fn) + ": system is unstable, rho >= 1");
    return rho;
}

/** rhohat = Q/(1+Q) with Q = sum_k lambda_k W_k. */
template <class T>
T mg1_discipline_rhohat(const std::vector<T>& lambda, const std::vector<T>& W) {
    T Q = num_traits<T>::from_int(0);
    for (std::size_t k = 0; k < W.size(); ++k) Q += lambda[k] * W[k];
    return Q / (num_traits<T>::from_int(1) + Q);
}

}  // namespace detail

/**
 * @param lambda per-class arrival rates
 * @param mu     per-class service rates
 * @param cs     per-class coefficients of variation of the service time
 */
template <class T>
Mg1DisciplineResult<T> qsys_mg1_setf(const std::vector<T>& lambda, const std::vector<T>& mu,
                                     const std::vector<T>& cs) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mg1_setf requires transcendental arithmetic");
    detail::mg1_discipline_check("qsys_mg1_setf", lambda, mu, cs);
    const std::size_t K = lambda.size();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);

    T lambda_total = zero;
    for (const T& v : lambda) lambda_total += v;

    // Mean residual service time of the mixture.
    T E_R = zero;
    for (std::size_t i = 0; i < K; ++i) {
        const T p_i = lambda[i] / lambda_total;
        const T E_S_i = one / mu[i];
        const T E_S2_i = (one + cs[i] * cs[i]) / (mu[i] * mu[i]);
        E_R += p_i * E_S2_i / (two * E_S_i);
    }

    Mg1DisciplineResult<T> r;
    r.W.assign(K, zero);
    for (std::size_t k = 0; k < K; ++k) {
        const T x = one / mu[k];
        T rho_x = zero, numer = zero;
        for (std::size_t i = 0; i < K; ++i) {
            T int_Fbar, int_tFbar;
            if (cs[i] == one) {
                const T e = detail::num_exp(T(-mu[i] * x));
                int_Fbar = (one - e) / mu[i];
                int_tFbar = (one - e * (one + mu[i] * x)) / (mu[i] * mu[i]);
            } else {
                int_Fbar = detail::num_min(x, T(one / mu[i]));
                int_tFbar = detail::num_min(T(x * x / two), T(one / (mu[i] * mu[i])));
            }
            rho_x += lambda[i] * int_Fbar;
            numer += lambda[i] * int_tFbar;
        }
        if (rho_x >= one)
            throw NumericError("qsys_mg1_setf: truncated load reaches one, the mean is infinite");
        r.W[k] = numer / ((one - rho_x) * (one - rho_x)) + x / (one - rho_x) + E_R / (one - rho_x);
    }
    r.rhohat = detail::mg1_discipline_rhohat(lambda, r.W);
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MG1_SETF_H
