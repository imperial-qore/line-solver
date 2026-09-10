/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_HARMONIC_H
#define LINE_API_FJ_HARMONIC_H

/**
 * Harmonic number H_K = sum_{k=1..K} 1/k.
 *
 * Templated port of matlab/src/api/fj/fj_harmonic.m, cross-checked against
 * jar/src/main/java/jline/api/fj/FJ_harmonic.java (identical).
 *
 * H_K is the expected maximum of K i.i.d. unit-rate exponentials and so is the
 * single most reused quantity in the fork-join family. It is a sum of unit
 * fractions, hence exactly representable in the field: the exact instantiation
 * returns the true rational H_K, which the double one does not (the summation
 * loses the low bits from about K = 10 onwards and the alternating sums in
 * fj_respt_vm and fj_xmax_erlang amplify that loss).
 */

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"

namespace line {
namespace fj {

/**
 * @param K number of parallel branches, K >= 1
 * @return  H_K = 1 + 1/2 + ... + 1/K
 */
template <class T>
T fj_harmonic(unsigned K) {
    detail::require_positive_K(K, "fj_harmonic");
    T H = num_traits<T>::from_int(0);
    for (unsigned k = 1; k <= K; ++k) H += num_traits<T>::from_int(1) / num_traits<T>::from_int(static_cast<long>(k));
    return H;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_HARMONIC_H
