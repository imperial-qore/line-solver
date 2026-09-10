/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MG1_PRIO_H
#define LINE_API_QSYS_QSYS_MG1_PRIO_H

/**
 * M/G/1 with non-preemptive head-of-line priorities: per-class mean response
 * times from the Cobham/Kleinrock formula.
 *
 * Templated port of matlab/src/api/qsys/qsys_mg1_prio.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mg1_prio.java (identical, including
 * the rhohat convention below).
 *
 *   rho_i = lambda_i/mu_i
 *   B_0   = (1/2) sum_i lambda_i (1 + cs_i^2)/mu_i^2      the mean work left
 *   Wq_k  = B_0 / [ (1 - sum_{i<k} rho_i)(1 - sum_{i<=k} rho_i) ]
 *   W_k   = Wq_k + 1/mu_k
 *
 * Class 1 is the highest priority. The returned rhohat is Q/(1+Q) with
 * Q = sum_k lambda_k W_k, the qsys family convention, not the utilization:
 * MATLAB computes the utilization into the same output variable and then
 * overwrites it, so the utilization is not observable from the return value.
 *
 * Only squares of the coefficients of variation appear, so the whole
 * computation is a rational function of the inputs and exact for
 * T = Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

template <class T>
struct Mg1PrioResult {
    std::vector<T> W;  ///< per-class mean response time, class 1 first
    T rhohat;          ///< Q/(1+Q) with Q = sum_k lambda_k W_k
};

/**
 * @param lambda per-class arrival rates, highest priority first
 * @param mu     per-class service rates
 * @param cs     per-class coefficients of variation of the service time
 */
template <class T>
Mg1PrioResult<T> qsys_mg1_prio(const std::vector<T>& lambda, const std::vector<T>& mu,
                               const std::vector<T>& cs) {
    const std::size_t K = lambda.size();
    if (mu.size() != K || cs.size() != K)
        throw InputError("qsys_mg1_prio: lambda, mu and cs must have the same length");
    if (K == 0) throw InputError("qsys_mg1_prio: at least one class is required");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    for (std::size_t i = 0; i < K; ++i)
        if (lambda[i] <= zero || mu[i] <= zero || cs[i] <= zero)
            throw InputError("qsys_mg1_prio: lambda, mu and cs must all be positive");

    std::vector<T> rho_i(K);
    T rho = zero;
    for (std::size_t i = 0; i < K; ++i) {
        rho_i[i] = lambda[i] / mu[i];
        rho += rho_i[i];
    }
    if (rho >= one) throw InputError("qsys_mg1_prio: system is unstable, rho >= 1");

    T B_0 = zero;
    for (std::size_t i = 0; i < K; ++i)
        B_0 += lambda[i] * (one + cs[i] * cs[i]) / (mu[i] * mu[i]);
    B_0 /= two;

    Mg1PrioResult<T> r;
    r.W.resize(K);
    T rho_prev = zero;
    for (std::size_t k = 0; k < K; ++k) {
        const T rho_curr = rho_prev + rho_i[k];
        r.W[k] = B_0 / ((one - rho_prev) * (one - rho_curr)) + one / mu[k];
        rho_prev = rho_curr;
    }
    T Q = zero;
    for (std::size_t k = 0; k < K; ++k) Q += lambda[k] * r.W[k];
    r.rhohat = Q / (one + Q);
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MG1_PRIO_H
