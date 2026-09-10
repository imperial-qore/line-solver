/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_XZABALOW_H
#define LINE_API_PFQN_PFQN_XZABALOW_H

/**
 * Asymptotic-bound-analysis lower bound on throughput.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_xzabalow.m. Single-class model: L is the
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

/** X >= N / (Z + N sum(L)), the ABA population bound. */
template <class T>
T pfqn_xzabalow(const std::vector<T>& L, const T& N, const T& Z) {
    if (L.empty()) throw InputError("pfqn_xzabalow: empty demand vector");
    T Ltot = num_traits<T>::from_int(0);
    for (const T& d : L) Ltot += d;
    return N / (Z + Ltot * N);
}

}  // namespace pfqn
}  // namespace line

#endif
