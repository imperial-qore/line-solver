/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_QZGBUP_H
#define LINE_API_PFQN_PFQN_QZGBUP_H

/**
 * Geometric-bound upper bound on the queue length at station i.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_qzgbup.m. Single-class model: L is the
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
#include "line/api/pfqn/pfqn_xzabaup.h"

namespace line {
namespace pfqn {

/** As the lower bound, with Y from the ABA upper bound and the sigma term. */
template <class T>
T pfqn_qzgbup(const std::vector<T>& L, const T& N, const T& Z, std::size_t i) {
    if (i >= L.size()) throw InputError("pfqn_qzgbup: station index out of range");
    const T one = num_traits<T>::from_int(1);
    T Ltot = num_traits<T>::from_int(0), L2 = num_traits<T>::from_int(0), Lmax = L[0];
    for (const T& d : L) {
        Ltot += d;
        L2 += d * d;
        if (d > Lmax) Lmax = d;
    }
    const T sigma = L2 / Ltot;
    const T Nm1 = N - one;
    const T xup = pfqn_xzabaup(L, Nm1, Z);
    const T cap = one / Lmax;
    const T alt = N / (Z + Ltot + sigma * (Nm1 - Z * xup));
    const T Yi = L[i] * (cap < alt ? cap : alt);
    if (Yi < one) {
        const long Nl = static_cast<long>(num_traits<T>::to_double(N));
        return Yi / (one - Yi) - num_pow_int(Yi, static_cast<unsigned>(Nl + 1)) / (one - Yi);
    }
    return N;
}

}  // namespace pfqn
}  // namespace line

#endif
