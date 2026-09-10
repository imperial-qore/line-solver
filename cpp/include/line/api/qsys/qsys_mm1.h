/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_MM1_H
#define LINE_API_QSYS_MM1_H

/**
 * Exact mean response time of the M/M/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_mm1.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_mm1.java (identical).
 *
 *   rho = lambda/mu,  W = rho/(1-rho)/lambda
 *
 * Pure field arithmetic, so the function is exact for T = Rational.
 */

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"

namespace line {
namespace qsys {

/**
 * @param lambda arrival rate
 * @param mu     service rate
 * @return W = mean response time, rhohat = rho (MATLAB returns the true rho here)
 */
template <class T>
QsysResult<T> qsys_mm1(const T& lambda, const T& mu) {
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu;
    detail::require_no_pole(T(one - rho), "qsys_mm1");
    const T W = rho / (one - rho) / lambda;
    return {W, rho};
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_MM1_H
