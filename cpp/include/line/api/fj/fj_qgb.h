/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_QGB_H
#define LINE_API_FJ_QGB_H

/**
 * Geometric bound on the queue length of a fork-join subnetwork.
 *
 * Templated port of matlab/src/api/fj/fj_qgb.m.
 *
 *   y_n(M) = D_n M / (Z + sum_j D_j H_{P_j} + Dmax M)
 *   Q_n(M) = H_{P_n} [ y_n/(1-y_n) - y_n^(M+1)/(1-y_n) ]
 *
 * The harmonic weights are what distinguishes this from pfqn_qzgblow: a P-way
 * fork-join subnetwork inflates its own demand by H_P in the denominator and
 * its queue length by H_P in the numerator, and setting every P_n to one
 * recovers the ordinary geometric bound exactly.
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [Q, y] of fj_qgb: the bounded queue lengths and the geometric ratios. */
template <class T>
struct FJQgbResult {
    std::vector<T> Q;
    std::vector<T> y;
};

/**
 * @param D per-visit service demand of each subnetwork
 * @param P fork degree of each subnetwork, P[n] >= 1
 * @param M number of circulating jobs, M >= 1
 * @param Z think time, Z >= 0
 * @return  the bounded queue lengths and the geometric ratios
 */
template <class T>
FJQgbResult<T> fj_qgb(const std::vector<T>& D, const std::vector<unsigned>& P, unsigned M,
                      const T& Z) {
    if (D.size() != P.size())
        throw InputError("fj_qgb: D and P must have the same number of elements");
    if (D.empty()) throw InputError("fj_qgb: at least one subnetwork is required");
    if (M < 1) throw InputError("fj_qgb: M must be a positive integer");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (Z < zero) throw InputError("fj_qgb: the think time must be non-negative");

    const std::size_t N = D.size();
    std::vector<T> H(N);
    T Dtot = zero, Dmax = zero;
    for (std::size_t n = 0; n < N; ++n) {
        if (D[n] < zero) throw InputError("fj_qgb: service demands must be non-negative");
        H[n] = fj_harmonic<T>(P[n]);
        Dtot += D[n] * H[n];
        if (n == 0 || D[n] > Dmax) Dmax = D[n];
    }

    FJQgbResult<T> out;
    out.Q.resize(N);
    out.y.resize(N);
    const T Mt = num_traits<T>::from_int(static_cast<long>(M));
    for (std::size_t n = 0; n < N; ++n) {
        out.y[n] = D[n] * Mt / (Z + Dtot + Dmax * Mt);
        if (out.y[n] < one) {
            T pw = one;
            for (unsigned e = 0; e <= M; ++e) pw *= out.y[n];
            out.Q[n] = H[n] * (out.y[n] / (one - out.y[n]) - pw / (one - out.y[n]));
        } else {
            // Degenerate ratio: the bound collapses onto the full population
            out.Q[n] = Mt;
        }
    }
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_QGB_H
