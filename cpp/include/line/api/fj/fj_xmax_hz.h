/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_HZ_H
#define LINE_API_FJ_XMAX_HZ_H

/**
 * Harrison-Zertal approximation of the maximum of i.i.d. variables.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_hz.m.
 *
 *   X_K^max ~ m1 + ( m2 / (2 m1) ) ( H_K - 1 )
 *
 * The correction is the equilibrium mean of the branch law scaled by H_K - 1:
 * one branch, plus the residual work still owed by the branches that finish
 * later. Writing m2/(2 m1) = m1 (1+SCV)/2 shows it is exact for the exponential
 * and reduces to m1 at K = 1 for every branch law.
 */

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [Xmax, resid] of fj_xmax_hz. */
template <class T>
struct FJXmaxHzResult {
    T Xmax;
    T resid;
};

/**
 * @param m1 mean of the branch distribution, m1 > 0
 * @param m2 second moment of the branch distribution, m2 >= m1^2
 * @param K  number of branches, K >= 1
 * @return   the approximate expected maximum and the equilibrium mean used
 */
template <class T>
FJXmaxHzResult<T> fj_xmax_hz(const T& m1, const T& m2, unsigned K) {
    const T zero = num_traits<T>::from_int(0), two = num_traits<T>::from_int(2);
    if (!(m1 > zero)) throw InputError("fj_xmax_hz: the branch mean must be positive");
    if (m2 < m1 * m1)
        throw InputError("fj_xmax_hz: the second moment is below the square of the mean");
    detail::require_positive_K(K, "fj_xmax_hz");
    FJXmaxHzResult<T> out;
    out.resid = m2 / (two * m1);
    out.Xmax = m1 + out.resid * (fj_harmonic<T>(K) - num_traits<T>::from_int(1));
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_HZ_H
