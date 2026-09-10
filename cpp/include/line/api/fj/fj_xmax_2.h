/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_2_H
#define LINE_API_FJ_XMAX_2_H

/**
 * Expected maximum of two independent, possibly unequal-rate exponentials.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_2.m, cross-checked against
 * FJ_xmax.fj_xmax_2 in jar/src/main/java/jline/api/fj/FJ_xmax.java
 * (identical).
 *
 *   Y_2^max = 1/lambda1 + 1/lambda2 - 1/(lambda1 + lambda2)
 *
 * Rational, hence exact in the field. At lambda1 = lambda2 = mu it collapses
 * to (3/2)/mu = H_2/mu, which is fj_xmax_exp at K = 2.
 */

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param lambda1 rate of the first branch, > 0
 * @param lambda2 rate of the second branch, > 0
 * @return        expected maximum of the two exponentials
 */
template <class T>
T fj_xmax_2(const T& lambda1, const T& lambda2) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (lambda1 <= zero || lambda2 <= zero) throw InputError("fj_xmax_2: the rates must be positive");
    return one / lambda1 + one / lambda2 - one / (lambda1 + lambda2);
}

/** Equal-rate overload, matching MATLAB's single-argument call. */
template <class T>
T fj_xmax_2(const T& lambda) {
    return fj_xmax_2(lambda, lambda);
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_2_H
