/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_RMAX_H
#define LINE_API_FJ_RMAX_H

/**
 * Pessimistic (independence) fork-join response time: the expected maximum of
 * K independent M/M/1 response times, each exponential with rate mu - lambda.
 *
 * Templated port of matlab/src/api/fj/fj_rmax.m, cross-checked against
 * FJ_rmax.fj_rmax in jar/src/main/java/jline/api/fj/FJ_rmax.java (identical;
 * the JAR additionally rejects rho >= 1 explicitly, where MATLAB leaves the
 * check to qsys_mm1).
 *
 *   Rmax = H_K R(rho) = H_K / (mu - lambda)
 *
 * Rational, hence exact in the field.
 */

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/num/number.h"

namespace line {
namespace fj {

/**
 * @param K      number of parallel branches, K >= 1
 * @param lambda arrival rate
 * @param mu     per-branch service rate
 * @return       H_K / (mu - lambda)
 */
template <class T>
T fj_rmax(unsigned K, const T& lambda, const T& mu) {
    detail::require_positive_K(K, "fj_rmax");
    return fj_harmonic<T>(K) * qsys::qsys_mm1(lambda, mu).W;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_RMAX_H
