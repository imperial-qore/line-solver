/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_XZGSBUP_H
#define LINE_API_PFQN_PFQN_XZGSBUP_H

/**
 * Geometric-square-root bound, upper bound on throughput.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_xzgsbup.m. Single-class model: L is the
 * per-station demand vector, N the population, Z the think time.
 *
 * Needs a square root, so it is available in double and high-precision
 * arithmetic only; the static_assert makes an exact instantiation a compile
 * error rather than a silent approximation.
 */

#include <algorithm>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/api/pfqn/pfqn_qzgbup.h"

namespace line {
namespace pfqn {

/** X = 2N / (R + sqrt(R^2 - 4 Z Lmax N)), R from the geometric queue bound. */
template <class T>
T pfqn_xzgsbup(const std::vector<T>& L, const T& N, const T& Z) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_xzgsbup requires transcendental arithmetic (square root)");
    const T one = num_traits<T>::from_int(1);
    T Ltot = num_traits<T>::from_int(0), Lmax = L[0];
    for (const T& d : L) {
        Ltot += d;
        if (d > Lmax) Lmax = d;
    }
    T R = Z + Ltot + Lmax * (N - one);
    for (std::size_t i = 0; i < L.size(); ++i)
        if (L[i] < Lmax) R += (L[i] - Lmax) * pfqn_qzgbup(L, T(N - one), Z, i);
    T disc = R * R - num_traits<T>::from_int(4) * Z * Lmax * N;
    if (disc < num_traits<T>::from_int(0)) disc = num_traits<T>::from_int(0);
    using std::sqrt;
    return num_traits<T>::from_int(2) * N / (R + sqrt(disc));
}

}  // namespace pfqn
}  // namespace line

#endif
