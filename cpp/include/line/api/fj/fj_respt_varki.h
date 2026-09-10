/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_RESPT_VARKI_H
#define LINE_API_FJ_RESPT_VARKI_H

/**
 * Varki approximation to the mean response time of a K-way fork-join system of
 * M/M/1 branches.
 *
 * Templated port of matlab/src/api/fj/fj_respt_varki.m, cross-checked against
 * FJ_respt.fj_respt_varki in jar/src/main/java/jline/api/fj/FJ_respt.java
 * (identical).
 *
 *   R_K = (1/mu) [ H_K + rho/(2(1-rho)) ( S1 + (1-2 rho) S2 ) ]
 *   S1 = sum_{i=1..K} 1/(i - rho),   S2 = sum_{i=1..K} 1/(i (i - rho))
 *
 * Rational in rho, hence exact in the field. Note that the denominators i - rho
 * are positive for every i >= 1 whenever rho < 1, so the only pole is at
 * rho = 1.
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
 * @param mu     per-branch service rate
 * @return       approximate mean fork-join response time
 */
template <class T>
T fj_respt_varki(unsigned K, const T& lambda, const T& mu) {
    detail::require_positive_K(K, "fj_respt_varki");
    const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
    const T rho = lambda / mu;
    if (rho >= one) throw NumericError("fj_respt_varki: unstable system, rho = lambda/mu >= 1");

    const T H_K = fj_harmonic<T>(K);
    T S1 = num_traits<T>::from_int(0), S2 = num_traits<T>::from_int(0);
    for (unsigned i = 1; i <= K; ++i) {
        const T ii = num_traits<T>::from_int(static_cast<long>(i));
        S1 += one / (ii - rho);
        S2 += one / (ii * (ii - rho));
    }
    return (one / mu) * (H_K + (rho / (two * (one - rho))) * (S1 + (one - two * rho) * S2));
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_RESPT_VARKI_H
