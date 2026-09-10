/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_NINTMVA_H
#define LINE_API_PFQN_NINTMVA_H

/**
 * Mean value analysis at a nonintegral population (fractional-base aMVA).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_nintmva.m.
 *
 * The exact MVA recursion started from the FRACTIONAL base n_0 = N - floor(N)
 * instead of from the empty network, giving mean performance measures of a
 * single-class closed product-form network at a real-valued population (Dowdy
 * and Gordon 1984, "aMVA"). The recursion is the standard Reiser-Lavenberg one,
 *
 *   R_i(n) = D_i (1 + Q_i(n-1)),  X(n) = n / (Z + sum_i R_i(n)),  Q_i(n) = X R_i,
 *
 * stepped in unit increments from n = n_0 (where the arrival-theorem term
 * Q_i(n_0 - 1) is taken as 0, the network below the base being empty) up to
 * n = N. At integer N the base is 0 and the recursion is bit-identical to exact
 * MVA; at fractional N it interpolates smoothly through the integral points,
 * which is what a nonintegral degree of multiprogramming (a time-average over a
 * measurement window) calls for.
 *
 * Unlike pfqn_dnc this accepts a think time, since the delay enters the
 * recursion and not a partial-fraction continuation. It is single-class: the
 * multiclass recursion has no one-dimensional step. For fractional multiclass
 * populations use pfqn_bs, which accepts them directly.
 *
 * Reference: L. W. Dowdy, K. D. Gordon, "Algorithms for Nonintegral Degrees of
 * Multiprogramming in Closed Queuing Networks", Performance Evaluation
 * 4(1):19-28, 1984.
 *
 * Arithmetic: EXACT-CAPABLE. Only field operations on T; the integer part of N
 * is read through a double, which is lossless for any population a closed model
 * can carry.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Mean performance measures of pfqn_nintmva at the requested population. */
template <class T>
struct NintMvaResult {
    T X;                ///< throughput at population N
    std::vector<T> Q;   ///< (M) mean queue lengths
    std::vector<T> U;   ///< (M) utilizations
    std::vector<T> R;   ///< (M) residence times
};

/**
 * @param L (M) service demands of the queueing stations
 * @param N population, real and nonnegative (may be fractional)
 * @param Z think time
 */
template <class T>
NintMvaResult<T> pfqn_nintmva(const std::vector<T>& L, const T& N, const T& Z) {
    const std::size_t M = L.size();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (N < zero) throw InputError("pfqn_nintmva requires a nonnegative population");

    NintMvaResult<T> res;
    res.Q.assign(M, zero);
    res.U.assign(M, zero);
    res.R.assign(M, zero);
    res.X = zero;
    if (N == zero) return res;

    // n_0 = N - floor(N), with the integral case starting the recursion at 1
    const double Nd = num_traits<T>::to_double(N);
    T n = T(N - num_traits<T>::from_int(static_cast<long>(std::floor(Nd))));
    if (n == zero) n = one;

    // The reference guards the real-valued loop bound with an absolute 1e-12
    const T stop = T(N + num_traits<T>::from_double(1e-12));
    while (n <= stop) {
        T sumR = zero;
        for (std::size_t i = 0; i < M; ++i) {
            res.R[i] = T(L[i] * T(one + res.Q[i]));
            sumR += res.R[i];
        }
        const T den = T(Z + sumR);
        if (den == zero) throw NumericError("pfqn_nintmva: the network has zero total demand");
        res.X = T(n / den);
        for (std::size_t i = 0; i < M; ++i) res.Q[i] = T(res.X * res.R[i]);
        n += one;
    }
    for (std::size_t i = 0; i < M; ++i) res.U[i] = T(res.X * L[i]);
    return res;
}

/** MATLAB default: no think time. */
template <class T>
NintMvaResult<T> pfqn_nintmva(const std::vector<T>& L, const T& N) {
    return pfqn_nintmva(L, N, num_traits<T>::from_int(0));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_NINTMVA_H
