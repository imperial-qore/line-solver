/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_COX_FIT_H
#define LINE_API_FJ_COX_FIT_H

/**
 * Two-stage Coxian fit of a mean and a squared coefficient of variation.
 *
 * Templated port of matlab/src/api/fj/fj_cox_fit.m.
 *
 * Marie's balanced-stage condition 1/mu1 = q/mu2 closes the system of two
 * moment equations in three unknowns and gives
 *
 *   mu1 = 2 mu,   q = 1/(2 c2),   mu2 = 2 mu q = mu/c2,
 *
 * which needs q <= 1, hence c2 >= 0.5. The Erlang stage count representing the
 * same target is bracketed by ceil(1/c2) <= k <= floor(1/c2) + 1.
 */

#include <cmath>

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [mu1, mu2, q, kmin, kmax] of fj_cox_fit. */
template <class T>
struct FJCoxFitResult {
    T mu1;
    T mu2;
    T q;
    unsigned kmin;
    unsigned kmax;
};

/**
 * @param m1 target mean, m1 > 0
 * @param c2 target squared coefficient of variation, c2 >= 0.5
 * @return   the two stage rates, the branching probability and the Erlang bracket
 */
template <class T>
FJCoxFitResult<T> fj_cox_fit(const T& m1, const T& c2) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1),
            two = num_traits<T>::from_int(2);
    if (!(m1 > zero)) throw InputError("fj_cox_fit: the target mean must be positive");
    if (c2 < one / two)
        throw InputError("fj_cox_fit: the balanced-stage Coxian fit needs c2 >= 0.5");

    const T mu = one / m1;
    FJCoxFitResult<T> out;
    // Both stages contribute half of the mean
    out.mu1 = two * mu;
    out.q = one / (two * c2);
    out.mu2 = mu / c2;

    const double inv = 1.0 / num_traits<T>::to_double(c2);
    long lo = static_cast<long>(std::ceil(inv - 1e-12));
    long hi = static_cast<long>(std::floor(inv + 1e-12)) + 1;
    if (lo < 1) lo = 1;
    if (hi < lo) hi = lo;
    out.kmin = static_cast<unsigned>(lo);
    out.kmax = static_cast<unsigned>(hi);
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_COX_FIT_H
