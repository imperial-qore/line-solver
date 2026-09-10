/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_RESPT_VM_H
#define LINE_API_FJ_RESPT_VM_H

/**
 * Varma-Makowski light-traffic interpolation for the mean response time of a
 * K-way fork-join system of M/M/1 branches.
 *
 * Templated port of matlab/src/api/fj/fj_respt_vm.m, cross-checked against
 * FJ_respt.fj_respt_vm in jar/src/main/java/jline/api/fj/FJ_respt.java
 * (identical).
 *
 *   R_K = [ H_K + (A_K - H_K) rho ] / (mu - lambda)
 *   A_K = sum_{i=1..K} C(K,i) (-1)^{i-1} sum_{m=1..i} C(i,m) (m-1)! / i^{m+1}
 *
 * Every term is rational, so the whole interpolation is exact in the field.
 * A_K is an alternating sum of binomial coefficients: at K = 30 the largest
 * term is about 1e8 times the result, so in double roughly eight significant
 * digits are lost to cancellation and by K = 60 nothing is left. The exact
 * instantiation is the only way to evaluate A_K at those K, and the only way
 * to measure the loss in the double one.
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
T fj_respt_vm(unsigned K, const T& lambda, const T& mu) {
    detail::require_positive_K(K, "fj_respt_vm");
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu;
    if (rho >= one) throw NumericError("fj_respt_vm: unstable system, rho = lambda/mu >= 1");

    const T H_K = fj_harmonic<T>(K);
    T A_K = num_traits<T>::from_int(0);
    for (unsigned i = 1; i <= K; ++i) {
        const T ii = num_traits<T>::from_int(static_cast<long>(i));
        T inner = num_traits<T>::from_int(0);
        for (unsigned m = 1; m <= i; ++m)
            inner += detail::fj_binom<T>(i, m) * num_factorial<T>(m - 1) / num_pow_int(ii, m + 1);
        const T term = detail::fj_binom<T>(K, i) * inner;
        if ((i - 1) % 2 == 0) A_K += term;
        else A_K -= term;
    }
    return (H_K + (A_K - H_K) * rho) / (mu - lambda);
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_RESPT_VM_H
