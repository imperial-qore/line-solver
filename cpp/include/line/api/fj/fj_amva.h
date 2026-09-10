/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_AMVA_H
#define LINE_API_FJ_AMVA_H

/**
 * Mean value analysis of a closed network of fork-join subnetworks.
 *
 * Templated port of matlab/src/api/fj/fj_amva.m.
 *
 *   R_n(m) = D_n [ H_{P_n} + Q_n(m-1) ]
 *   X(m)   = m / (Z + sum_n R_n(m))
 *   Q_n(m) = X(m) R_n(m)
 *
 * started from Q_n(0) = 0. With every P_n equal to one this is the exact
 * single-class mean value analysis, because H_1 = 1; above that it is an
 * approximation whose per-subnetwork residence time is an upper bound in the
 * sense of Varki.
 */

#include <cstddef>
#include <vector>

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [R, Q, X, U] of fj_amva. */
template <class T>
struct FJAmvaResult {
    std::vector<T> R;
    std::vector<T> Q;
    T X;
    std::vector<T> U;
};

/**
 * @param D per-visit service demand of each subnetwork
 * @param P fork degree of each subnetwork, P[n] >= 1
 * @param M number of circulating jobs, M >= 1
 * @param Z think time, Z >= 0
 * @return  residence times, queue lengths, throughput and per-queue utilizations
 */
template <class T>
FJAmvaResult<T> fj_amva(const std::vector<T>& D, const std::vector<unsigned>& P, unsigned M,
                        const T& Z) {
    if (D.size() != P.size())
        throw InputError("fj_amva: D and P must have the same number of elements");
    if (D.empty()) throw InputError("fj_amva: at least one subnetwork is required");
    if (M < 1) throw InputError("fj_amva: M must be a positive integer");
    const T zero = num_traits<T>::from_int(0);
    if (Z < zero) throw InputError("fj_amva: the think time must be non-negative");

    const std::size_t N = D.size();
    std::vector<T> H(N);
    for (std::size_t n = 0; n < N; ++n) {
        if (D[n] < zero) throw InputError("fj_amva: service demands must be non-negative");
        H[n] = fj_harmonic<T>(P[n]);
    }

    FJAmvaResult<T> out;
    out.R.assign(N, zero);
    out.Q.assign(N, zero);
    out.U.assign(N, zero);
    out.X = zero;
    for (unsigned m = 1; m <= M; ++m) {
        T Rtot = zero;
        for (std::size_t n = 0; n < N; ++n) {
            out.R[n] = D[n] * (H[n] + out.Q[n]);
            Rtot += out.R[n];
        }
        if (!(Rtot > zero))
            throw NumericError("fj_amva: the total residence time vanished; every demand is zero");
        out.X = num_traits<T>::from_int(static_cast<long>(m)) / (Z + Rtot);
        for (std::size_t n = 0; n < N; ++n) out.Q[n] = out.X * out.R[n];
    }
    // Each subnetwork holds P(n) queues sharing the demand equally
    for (std::size_t n = 0; n < N; ++n)
        out.U[n] = out.X * D[n] / num_traits<T>::from_int(static_cast<long>(P[n]));
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_AMVA_H
