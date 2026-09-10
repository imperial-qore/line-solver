/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_AQL_H
#define LINE_API_PFQN_AQL_H

/**
 * Aggregate Queue Length (AQL) approximate MVA.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_aql.m. AQL solves K+1 coupled
 * populations at once: the full population N and each N - e_s. The arrival
 * estimate is
 *   R(k,s|n) = L(k,s) [1 + (|n| - 1) (Q(k|n)/|n| - gamma(k,s))]
 * with the correction gamma(k,s) = Q(k|N)/|N| - Q(k|N - e_s)/(|N| - 1)
 * refreshed after every sweep. The aggregate queue length Q(k|n) is kept per
 * population rather than per class, which is what distinguishes AQL from
 * Linearizer and makes it cheaper by a factor of the class count.
 *
 * Iterates to a relative tolerance, so exact arithmetic buys nothing: the
 * static_assert records that.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_bs.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param L  (M x K) demands
 * @param N  (K) populations
 * @param Z  (K) think times, empty for none
 * @param tol convergence tolerance
 * @param maxiter iteration cap
 * @return the standard AMVA metrics; AN holds the arrival-instant aggregate
 *         queue lengths Q(k | N - e_s), which callers use for the arrival
 *         theorem diagnostics
 */
template <class T>
AmvaResult<T> pfqn_aql(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                       double tol = 1e-7, std::size_t maxiter = 1000) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_aql requires transcendental arithmetic: it iterates to a relative "
                  "tolerance, so its answer is a fixed point only to within tol");
    const std::size_t M = L.rows(), K = L.cols();
    if (N.size() != K) throw InputError("pfqn_aql: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != K) throw InputError("pfqn_aql: Z has the wrong length");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    T Ntot = zero;
    for (const T& v : N) Ntot += v;

    AmvaResult<T> out;
    out.XN.assign(K, zero);
    out.QN = Matrix<T>(M, K, zero);
    out.UN = Matrix<T>(M, K, zero);
    out.RN = Matrix<T>(M, K, zero);
    if (M == 0 || Ntot == zero) return out;

    // Populations: index 0 is N, index s+1 is N - e_s.
    std::vector<std::vector<T>> pops(K + 1, N);
    for (std::size_t s = 0; s < K; ++s)
        pops[s + 1][s] = (N[s] > zero) ? N[s] - one : zero;

    std::vector<std::vector<T>> Q(K + 1, std::vector<T>(M, zero));
    std::vector<Matrix<T>> R(K + 1, Matrix<T>(M, K, zero));
    std::vector<std::vector<T>> X(K + 1, std::vector<T>(K, zero));
    Matrix<T> gamma(M, K, zero);

    for (std::size_t t = 0; t <= K; ++t)
        for (std::size_t k = 0; k < M; ++k)
            Q[t][k] = Ntot / num_traits<T>::from_int(static_cast<long>(M));

    for (std::size_t it = 1; it <= maxiter; ++it) {
        out.iterations = it;
        const std::vector<T> Qprev = Q[0];

        for (std::size_t t = 0; t <= K; ++t) {
            const std::vector<T>& n = pops[t];
            T ntot = zero;
            for (const T& v : n) ntot += v;
            if (ntot == zero) {
                for (std::size_t k = 0; k < M; ++k) Q[t][k] = zero;
                for (std::size_t s = 0; s < K; ++s) X[t][s] = zero;
                continue;
            }
            for (std::size_t k = 0; k < M; ++k)
                for (std::size_t s = 0; s < K; ++s)
                    R[t](k, s) = L(k, s) * (one + (ntot - one) * (Q[t][k] / ntot - gamma(k, s)));
            for (std::size_t s = 0; s < K; ++s) {
                T denom = Z.empty() ? zero : Z[s];
                for (std::size_t k = 0; k < M; ++k) denom += R[t](k, s);
                X[t][s] = (denom == zero) ? zero : n[s] / denom;
            }
            for (std::size_t k = 0; k < M; ++k) {
                T q = zero;
                for (std::size_t s = 0; s < K; ++s) q += X[t][s] * R[t](k, s);
                Q[t][k] = q;
            }
        }

        if (Ntot > one)
            for (std::size_t k = 0; k < M; ++k)
                for (std::size_t s = 0; s < K; ++s)
                    gamma(k, s) = Q[0][k] / Ntot - Q[s + 1][k] / (Ntot - one);

        double delta = 0.0;
        for (std::size_t k = 0; k < M; ++k) {
            if (Q[0][k] == zero) continue;
            const double d = std::fabs(num_traits<T>::to_double(T((Qprev[k] - Q[0][k]) / Q[0][k])));
            if (d > delta) delta = d;
        }
        if (delta < tol) {
            out.converged = true;
            break;
        }
    }

    out.XN = X[0];
    out.RN = R[0];
    for (std::size_t k = 0; k < M; ++k)
        for (std::size_t s = 0; s < K; ++s) {
            out.UN(k, s) = out.XN[s] * L(k, s);
            out.QN(k, s) = out.UN(k, s) * (one + Q[s + 1][k]);
        }
    return out;
}

template <class T>
AmvaResult<T> pfqn_aql(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_aql(L, N, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_AQL_H
