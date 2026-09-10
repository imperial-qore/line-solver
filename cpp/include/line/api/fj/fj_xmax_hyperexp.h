/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_XMAX_HYPEREXP_H
#define LINE_API_FJ_XMAX_HYPEREXP_H

/**
 * Expected maximum of K i.i.d. two-phase hyperexponential service times.
 *
 * Templated port of matlab/src/api/fj/fj_xmax_hyperexp.m.
 *
 *   X_K^max = sum_{n=1..K} (-1)^{n+1} sum_{m=0..n}
 *                 C(n,m) p1^m p2^{n-m} / (m mu1 + (n-m) mu2)
 *
 * from the inclusion-exclusion expansion of 1 - F(x)^K. Every term is
 * rational, so the sum is exact in the field -- which is the point, because
 * it alternates: the terms grow like 2^K while the result stays O(log K/mu),
 * so in double the answer is destroyed by cancellation somewhere around
 * K = 25 and is pure noise by K = 40.
 *
 * REFERENCE DEFECT: FJ_xmax.fj_xmax_hyperexp in
 * jar/src/main/java/jline/api/fj/FJ_xmax.java drops the p1^m p2^{n-m} factor
 * from the inner sum, so the JAR computes a different quantity and ignores p1
 * entirely except in its validation. MATLAB is ground truth and is what this
 * port follows.
 */

#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/**
 * @param K   number of branches, K >= 1
 * @param p1  probability of the first phase, 0 < p1 < 1
 * @param mu1 rate of the first phase, > 0
 * @param mu2 rate of the second phase, > 0
 * @return    expected maximum of K hyperexponential samples
 */
template <class T>
T fj_xmax_hyperexp(unsigned K, const T& p1, const T& mu1, const T& mu2) {
    detail::require_positive_K(K, "fj_xmax_hyperexp");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (p1 <= zero || p1 >= one) throw InputError("fj_xmax_hyperexp: p1 must lie in (0,1)");
    if (mu1 <= zero || mu2 <= zero) throw InputError("fj_xmax_hyperexp: the rates mu1 and mu2 must be positive");
    const T p2 = one - p1;

    T Xmax = zero;
    for (unsigned n = 1; n <= K; ++n) {
        T inner = zero;
        for (unsigned m = 0; m <= n; ++m) {
            const T den = num_traits<T>::from_int(static_cast<long>(m)) * mu1 +
                          num_traits<T>::from_int(static_cast<long>(n - m)) * mu2;
            if (den > zero)
                inner += detail::fj_binom<T>(n, m) * num_pow_int(p1, m) * num_pow_int(p2, n - m) / den;
        }
        if ((n + 1) % 2 == 0) Xmax += inner;
        else Xmax -= inner;
    }
    return Xmax;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_XMAX_HYPEREXP_H
