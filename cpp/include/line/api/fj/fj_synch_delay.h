/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_SYNCH_DELAY_H
#define LINE_API_FJ_SYNCH_DELAY_H

/**
 * Mean synchronization delay of a 2-way fork-join system of M/M/1 branches,
 * i.e. the time the first-finishing branch waits at the join.
 *
 * Templated port of matlab/src/api/fj/fj_synch_delay.m, cross-checked against
 * jar/src/main/java/jline/api/fj/FJ_synch_delay.java (identical).
 *
 *   S_2 = (1/2)(1 - rho/4) R(rho)
 *
 * Rational in rho, hence exact in the field. It satisfies
 * R_2 = R + S_2 exactly, which is the identity worth checking: both sides are
 * exact rationals, so any discrepancy is a port error and not a rounding one.
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
 * @return       mean synchronization delay
 */
template <class T>
T fj_synch_delay(const T& lambda, const T& mu) {
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu;
    if (rho >= one) throw NumericError("fj_synch_delay: unstable system, rho = lambda/mu >= 1");
    const T R_rho = qsys::qsys_mm1(lambda, mu).W;
    return num_traits<T>::from_rational(1, 2) * (one - rho / num_traits<T>::from_int(4)) * R_rho;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_SYNCH_DELAY_H
