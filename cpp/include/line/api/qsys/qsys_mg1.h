/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_MG1_H
#define LINE_API_QSYS_MG1_H

/**
 * Exact mean response time of the M/G/1 queue (Pollaczek-Khinchine).
 *
 * Templated port of matlab/src/api/qsys/qsys_mg1.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mg1.java (identical).
 *
 *   rho = lambda/mu
 *   Q   = rho + rho^2/(2(1-rho)) + lambda^2 cs^2 / mu^2 / (2(1-rho))
 *   W   = Q/lambda,  rhohat = Q/(1+Q)
 *
 * Only integer powers appear, so the function is exact for T = Rational.
 * Note that cs is the coefficient of variation, squared inside the formula.
 */

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"

namespace line {
namespace qsys {

/**
 * @param lambda arrival rate
 * @param mu     service rate
 * @param cs     coefficient of variation of the service time
 * @return W = mean response time, rhohat = Q/(1+Q)
 */
template <class T>
QsysResult<T> qsys_mg1(const T& lambda, const T& mu, const T& cs) {
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    detail::require_no_pole(T(one - rho), "qsys_mg1");
    const T Q = rho + num_pow_int(rho, 2) / (two * (one - rho)) +
                num_pow_int(lambda, 2) * num_pow_int(cs, 2) / num_pow_int(mu, 2) /
                    (two * (one - rho));
    const T W = Q / lambda;
    return {W, Q / (one + Q)};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_MG1_H
