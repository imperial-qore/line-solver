/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_QZGBLOW_H
#define LINE_API_PFQN_PFQN_QZGBLOW_H

/**
 * Geometric-bound lower bound on the queue length at station i.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_qzgblow.m. Single-class model: L is the
 * per-station demand vector, N the population, Z the think time.
 *
 * All operations stay in the field, so the bound is exact in rational
 * arithmetic: a bound computed exactly is worth having, since a bound violated
 * only by rounding is indistinguishable from a real violation.
 */

#include <algorithm>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Qgb = y/(1-y) - y^(N+1)/(1-y) with y = N L_i / (Z + sum(L) + Lmax N). */
template <class T>
T pfqn_qzgblow(const std::vector<T>& L, const T& N, const T& Z, std::size_t i) {
    if (i >= L.size()) throw InputError("pfqn_qzgblow: station index out of range");
    T Ltot = num_traits<T>::from_int(0), Lmax = L[0];
    for (const T& d : L) {
        Ltot += d;
        if (d > Lmax) Lmax = d;
    }
    const T yi = N * L[i] / (Z + Ltot + Lmax * N);
    const T one = num_traits<T>::from_int(1);
    if (yi == one) throw NumericError("pfqn_qzgblow: degenerate geometric ratio y = 1");
    // N enters as an integer exponent, so this stays in the field.
    const long Nl = static_cast<long>(num_traits<T>::to_double(N));
    return yi / (one - yi) - num_pow_int(yi, static_cast<unsigned>(Nl + 1)) / (one - yi);
}

}  // namespace pfqn
}  // namespace line

#endif
