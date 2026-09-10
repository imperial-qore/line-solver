/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_MMK_H
#define LINE_API_QSYS_MMK_H

/**
 * Exact mean response time of the M/M/k queue (Erlang-C).
 *
 * Templated port of matlab/src/api/qsys/qsys_mmk.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mmk.java. The JAR accumulates the
 * factorials incrementally and writes C = 1/(1+(1-rho)*(C*(C-1)!)/(C rho)^C*S),
 * which is the same k! as MATLAB's factorial(k); the two agree.
 *
 *   rho = lambda/(mu k)
 *   S   = sum_{j=0}^{k-1} (k rho)^j / j!
 *   C   = 1 / (1 + (1-rho) k! / (k rho)^k * S)
 *   Q   = rho/(1-rho) C + k rho,   W = Q/lambda
 *
 * Every exponent is an integer, so this stays in the field of the inputs and
 * is exact for T = Rational.
 */

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace qsys {

namespace detail {

/**
 * Erlang-C: the probability that an arriving customer finds all k servers
 * busy. Argument rho is the per-server utilization lambda/(mu k).
 */
template <class T>
T erlang_c(unsigned k, const T& rho) {
    const T one = num_traits<T>::from_int(1);
    const T krho = num_traits<T>::from_int(static_cast<long>(k)) * rho;
    T S = num_traits<T>::from_int(0);
    for (unsigned j = 0; j < k; ++j) S += num_pow_int(krho, j) / num_factorial<T>(j);
    return one / (one + (one - rho) * num_factorial<T>(k) / num_pow_int(krho, k) * S);
}

}  // namespace detail

/**
 * @param lambda arrival rate
 * @param mu     service rate of a single server
 * @param k      number of servers, k >= 1
 * @return W = mean response time, rhohat = rho = lambda/(mu k)
 */
template <class T>
QsysResult<T> qsys_mmk(const T& lambda, const T& mu, unsigned k) {
    if (k == 0) throw InputError("qsys_mmk: k must be at least 1");
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu / num_traits<T>::from_int(static_cast<long>(k));
    detail::require_no_pole(T(one - rho), "qsys_mmk");
    const T Q = rho / (one - rho) * detail::erlang_c(k, rho) +
                num_traits<T>::from_int(static_cast<long>(k)) * rho;
    const T W = Q / lambda;
    return {W, rho};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_MMK_H
