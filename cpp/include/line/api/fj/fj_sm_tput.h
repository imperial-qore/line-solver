/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_SM_TPUT_H
#define LINE_API_FJ_SM_TPUT_H

/**
 * Saturated (single-message) maximum throughput of a K-way fork-join system
 * with exponential branch service.
 *
 * Templated port of matlab/src/api/fj/fj_sm_tput.m, cross-checked against
 * jar/src/main/java/jline/api/fj/FJ_sm_tput.java (identical).
 *
 *   lambda_max = mu / H_K = 1 / X_K^max
 *
 * Rational, hence exact in the field, and the exact reciprocal of fj_xmax_exp
 * by construction.
 */

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K  number of parallel branches, K >= 1
 * @param mu per-branch service rate, mu > 0
 * @return   maximum sustainable arrival rate
 */
template <class T>
T fj_sm_tput(unsigned K, const T& mu) {
    detail::require_positive_K(K, "fj_sm_tput");
    if (mu <= num_traits<T>::from_int(0)) throw InputError("fj_sm_tput: the service rate mu must be positive");
    return mu / fj_harmonic<T>(K);
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_SM_TPUT_H
