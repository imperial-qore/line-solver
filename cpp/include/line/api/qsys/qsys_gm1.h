/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_GM1_H
#define LINE_API_QSYS_GM1_H

/**
 * Exact mean response time of the G/M/1 queue.
 *
 * Templated port of matlab/src/api/qsys/qsys_gm1.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gm1.java. The two agree on W; the JAR
 * additionally reports rhohat = 0, whereas MATLAB returns W alone. The port
 * follows MATLAB and returns the single value.
 *
 *   W = 1/(1-sigma)/mu
 *
 * sigma is the load seen at arrival instants, i.e. the root in (0,1) of
 * sigma = A*(mu(1-sigma)) with A* the Laplace transform of the interarrival
 * time. It is an input here, so the function itself is pure field arithmetic
 * and exact for T = Rational.
 */

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"

namespace line {
namespace qsys {

/**
 * @param sigma load at arrival instants, 0 <= sigma < 1
 * @param mu    service rate
 * @return mean response time W
 */
template <class T>
T qsys_gm1(const T& sigma, const T& mu) {
    const T one = num_traits<T>::from_int(1);
    detail::require_no_pole(T(one - sigma), "qsys_gm1");
    return one / (one - sigma) / mu;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_GM1_H
