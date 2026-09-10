/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_EXP_H
#define LINE_API_FJ_XMAX_EXP_H

/**
 * Expected maximum of K i.i.d. exponential service times.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_exp.m, cross-checked against
 * FJ_xmax.fj_xmax_exp in jar/src/main/java/jline/api/fj/FJ_xmax.java
 * (identical).
 *
 *   X_K^max = H_K / mu
 *
 * Exact in closed form and rational, so this is the oracle the rest of the
 * xmax family is checked against: fj_xmax_2 at equal rates, fj_xmax_hyperexp
 * at mu1 = mu2 and fj_xmax_erlang at k = 1 must all reproduce it exactly.
 */

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K  number of branches, K >= 1
 * @param mu branch service rate, mu > 0
 * @return   H_K / mu
 */
template <class T>
T fj_xmax_exp(unsigned K, const T& mu) {
    detail::require_positive_K(K, "fj_xmax_exp");
    if (mu <= num_traits<T>::from_int(0)) throw InputError("fj_xmax_exp: the service rate mu must be positive");
    return fj_harmonic<T>(K) / mu;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_EXP_H
