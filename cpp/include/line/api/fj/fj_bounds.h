/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_BOUNDS_H
#define LINE_API_FJ_BOUNDS_H

/**
 * Upper and lower bounds on the mean response time of a K-way fork-join system
 * of M/M/1 branches.
 *
 * Templated port of matlab/src/api/fj/fj_bounds.m, cross-checked against
 * jar/src/main/java/jline/api/fj/FJ_bounds.java (identical).
 *
 *   Rmax = H_K / (mu (1 - rho))                       (Thomasian 2014, Eq. 1)
 *   Rmin = (1/mu) [ H_K + sum_{j=1..K} (1/j) rho/(j - rho) ]          (Eq. 2)
 *
 * Both are rational functions of rho = lambda/mu, so the pair is exact in the
 * field. That matters because the bounds are meant to bracket the true mean:
 * a rounded Rmin can exceed a rounded Rmax when the two are close, which the
 * exact instantiation never does.
 */

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K      number of parallel branches, K >= 1
 * @param lambda arrival rate
 * @param mu     per-branch service rate, mu > lambda for stability
 * @return       Rmax (pessimistic) and Rmin (optimistic) bounds
 */
template <class T>
FJBoundsResult<T> fj_bounds(unsigned K, const T& lambda, const T& mu) {
    detail::require_positive_K(K, "fj_bounds");
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu;
    if (rho >= one) throw NumericError("fj_bounds: unstable system, rho = lambda/mu >= 1");

    const T H_K = fj_harmonic<T>(K);
    const T Rmax = H_K / (mu * (one - rho));

    T S_K = num_traits<T>::from_int(0);
    for (unsigned j = 1; j <= K; ++j) {
        const T jj = num_traits<T>::from_int(static_cast<long>(j));
        S_K += (one / jj) * (rho / (jj - rho));
    }
    const T Rmin = (one / mu) * (H_K + S_K);
    return {Rmax, Rmin};
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_BOUNDS_H
