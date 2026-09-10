/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_XZABAUP_H
#define LINE_API_PFQN_PFQN_XZABAUP_H

/**
 * Asymptotic-bound-analysis upper bound on throughput.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_xzabaup.m. Single-class model: L is the
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

/** X <= min(1/Lmax, N/(sum(L)+Z)): capacity bound and population bound. */
template <class T>
T pfqn_xzabaup(const std::vector<T>& L, const T& N, const T& Z) {
    if (L.empty()) throw InputError("pfqn_xzabaup: empty demand vector");
    T Ltot = num_traits<T>::from_int(0), Lmax = L[0];
    for (const T& d : L) {
        Ltot += d;
        if (d > Lmax) Lmax = d;
    }
    if (Lmax == num_traits<T>::from_int(0)) throw InputError("pfqn_xzabaup: all demands are zero");
    const T cap = num_traits<T>::from_int(1) / Lmax;
    const T pop = N / (Ltot + Z);
    return cap < pop ? cap : pop;
}

}  // namespace pfqn
}  // namespace line

#endif
