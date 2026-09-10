/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_RESPT_NT_H
#define LINE_API_FJ_RESPT_NT_H

/**
 * Nelson-Tantawi approximation to the mean response time of a K-way fork-join
 * system of M/M/1 branches.
 *
 * Templated port of matlab/src/api/fj/fj_respt_nt.m, cross-checked against
 * FJ_respt.fj_respt_nt in jar/src/main/java/jline/api/fj/FJ_respt.java
 * (identical apart from the K > 32 accuracy warning, which MATLAB emits and
 * neither the JAR nor this port does).
 *
 *   R_K = [ H_K/H_2 + (1 - H_K/H_2) 4 rho/11 ] (3/2 - rho/8) / (mu - lambda)
 *
 * Rational in rho, hence exact in the field. At K = 2 it reduces exactly to
 * fj_respt_2way, an identity that only holds bit-for-bit in exact arithmetic.
 */

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K      number of parallel branches, K >= 2
 * @param lambda arrival rate
 * @param mu     per-branch service rate
 * @return       approximate mean fork-join response time
 */
template <class T>
T fj_respt_nt(unsigned K, const T& lambda, const T& mu) {
    if (K < 2) throw InputError("fj_respt_nt: the Nelson-Tantawi approximation requires K >= 2");
    const T one = num_traits<T>::from_int(1);
    const T rho = lambda / mu;
    if (rho >= one) throw NumericError("fj_respt_nt: unstable system, rho = lambda/mu >= 1");

    const T H_K = fj_harmonic<T>(K);
    const T H_2 = fj_harmonic<T>(2);
    const T ratio = H_K / H_2;
    const T S_K = ratio + (one - ratio) * (num_traits<T>::from_int(4) * rho / num_traits<T>::from_int(11));
    const T R2f = num_traits<T>::from_rational(3, 2) - rho / num_traits<T>::from_int(8);
    return S_K * R2f / (mu - lambda);
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_RESPT_NT_H
