/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_RESPT_2WAY_H
#define LINE_API_FJ_RESPT_2WAY_H

/**
 * Exact mean response time of a 2-way fork-join system of M/M/1 branches.
 *
 * Templated port of matlab/src/api/fj/fj_respt_2way.m, cross-checked against
 * FJ_respt.fj_respt_2way in jar/src/main/java/jline/api/fj/FJ_respt.java
 * (identical; the JAR inlines qsys_mm1 as 1/(mu-lambda), which is the same
 * value MATLAB's qsys_mm1 returns).
 *
 *   R_2 = (H_2 - rho/8) R(rho) = ((12 - rho)/8) / (mu - lambda)
 *
 * Rational in rho, hence exact in the field. This is the only fork-join
 * response time in the family that is exact rather than an approximation, so
 * it is the natural calibration point for the approximations around it.
 */

#include "line/api/fj/fj_types.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param lambda arrival rate
 * @param mu     per-branch service rate
 * @return       mean 2-way fork-join response time
 */
template <class T>
T fj_respt_2way(const T& lambda, const T& mu) {
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu;
    if (rho >= one) throw NumericError("fj_respt_2way: unstable system, rho = lambda/mu >= 1");
    const T R_rho = qsys::qsys_mm1(lambda, mu).W;
    return (num_traits<T>::from_int(12) - rho) / num_traits<T>::from_int(8) * R_rho;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_RESPT_2WAY_H
