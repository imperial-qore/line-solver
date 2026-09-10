/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MG1K_LOSS_MGS_H
#define LINE_API_QSYS_QSYS_MG1K_LOSS_MGS_H

/**
 * MacGregor Smith's closed-form approximation of the M/G/1/K loss probability.
 *
 * Templated port of matlab/src/api/qsys/qsys_mg1k_loss_mgs.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_mg1k_loss_mgs.java
 * (identical).
 *
 *   rho = lambda/mu,  s = sqrt(scv),  r = sqrt(rho),  b = 2 + r s^2 - r
 *   Ploss = rho^((r s^2 - r + 2K)/b) (rho - 1) / ( rho^(2(1 + r s^2 - r + K)/b) - 1 )
 *
 * The exponents interpolate the exact M/M/1/K expression in the service
 * variability: at scv = 1 they are K and K+1 and the formula reduces to
 * qsys_mm1k_loss.
 *
 * ARITHMETIC. Both the square root of rho and the real-valued exponents are
 * transcendental, so the function is gated. There is no exact instantiation
 * even in principle: for scv != 1 the exponents are irrational.
 *
 * Reference: J. MacGregor Smith, Optimal design and performance modelling of
 * M/G/1/K queueing systems.
 */

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

template <class T>
struct Mg1kLossMgsResult {
    T lossProbability;
    T utilization;  ///< rho = lambda/mu
};

/**
 * @param lambda arrival rate
 * @param mu     service rate
 * @param mu_scv squared coefficient of variation of the service time
 * @param K      system capacity, jobs in service included
 */
template <class T>
Mg1kLossMgsResult<T> qsys_mg1k_loss_mgs(const T& lambda, const T& mu, const T& mu_scv,
                                        unsigned K) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mg1k_loss_mgs requires transcendental arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    const T s = detail::num_sqrt(mu_scv);
    const T sqrt_rho = detail::num_sqrt(rho);
    const T Kt = num_traits<T>::from_int(static_cast<long>(K));
    const T b = two + sqrt_rho * s * s - sqrt_rho;
    if (b == num_traits<T>::from_int(0))
        throw NumericError("qsys_mg1k_loss_mgs: degenerate exponent denominator");
    const T num = detail::num_pow(rho, T((sqrt_rho * s * s - sqrt_rho + two * Kt) / b)) * (rho - one);
    const T den = detail::num_pow(rho, T(two * (one + sqrt_rho * s * s - sqrt_rho + Kt) / b)) - one;
    if (den == num_traits<T>::from_int(0))
        throw NumericError("qsys_mg1k_loss_mgs: rho == 1, the closed form is singular");
    Mg1kLossMgsResult<T> r;
    r.utilization = rho;
    r.lossProbability = num / den;
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MG1K_LOSS_MGS_H
