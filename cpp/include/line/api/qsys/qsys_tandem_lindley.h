/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_TANDEM_LINDLEY_H
#define LINE_API_QSYS_QSYS_TANDEM_LINDLEY_H

/**
 * Sample-path Lindley recursion along a tandem of single-server FCFS queues.
 *
 * Templated port of matlab/src/api/qsys/qsys_tandem_lindley.m. No JAR
 * counterpart. Given the interarrival times A of N customers at the first
 * station and the (N x K) service times S, it replays
 *
 *   W(n+1,k) = max(W(n,k) + S(n,k) - G(n,k), 0)
 *
 * where G(n,k) is the interarrival gap SEEN AT STATION k. The coupling that
 * makes this a tandem rather than K independent queues is the gap update
 *
 *   G(n,k+1) = max(G(n,k) - W(n,k) - S(n,k), 0) + S(n+1,k),
 *
 * the interdeparture time of station k: customer n+1 either catches up with
 * customer n, in which case the gap collapses to the service time of n+1 at
 * station k, or it does not, and the residual gap survives. Nothing here is
 * distributional, so the recursion holds for arbitrary service laws and is
 * the reference any conditional-moment approximation is checked against.
 *
 * The last customer has no successor, so its row of G stays unset; MATLAB
 * leaves it NaN and the port does the same rather than filling it with a
 * value that has no meaning.
 *
 * Reference: S. Palomo, J. Pender, "Learning the Tandem Network Lindley
 * Recursion", Proc. Winter Simulation Conference, 2021. Registered in
 * .citations() as 'tandemlindley'.
 *
 * ARITHMETIC: additions, subtractions and comparisons only, so the exact
 * instantiation replays the same sample path with no rounding at all.
 */

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Mirrors the struct MATLAB returns from qsys_tandem_lindley. */
template <class T>
struct TandemLindleyResult {
    Matrix<T> W;          ///< (N x K) waiting time of each customer at each station
    Matrix<T> G;          ///< (N x K) gap seen at each station, last row unset (NaN)
    Matrix<T> T_;         ///< (N x K) sojourn time W + S
    Matrix<T> departure;  ///< (N x K) departure epoch from each station
    std::string analyzer;
};

/**
 * @param A  (N) interarrival times at the first station, finite nonnegative
 * @param S  (N x K) service times, finite nonnegative
 * @param W0 (K) initial waiting times, one per station; empty means all zero
 */
template <class T>
TandemLindleyResult<T> qsys_tandem_lindley(const std::vector<T>& A, const Matrix<T>& S,
                                           const std::vector<T>& W0 = std::vector<T>()) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t N = A.size();
    if (N < 1) throw InputError("qsys_tandem_lindley: A must hold at least one interarrival time");
    if (S.rows() != N)
        throw InputError("qsys_tandem_lindley: S must have one row per interarrival time");
    const std::size_t K = S.cols();
    if (K < 1) throw InputError("qsys_tandem_lindley: S must have at least one column");
    for (std::size_t n = 0; n < N; ++n) {
        if (A[n] < zero)
            throw InputError("qsys_tandem_lindley: A must hold nonnegative interarrival times");
        for (std::size_t k = 0; k < K; ++k)
            if (S(n, k) < zero)
                throw InputError("qsys_tandem_lindley: S must hold nonnegative service times");
    }
    std::vector<T> w0 = W0;
    if (w0.empty()) w0.assign(K, zero);
    if (w0.size() != K)
        throw InputError("qsys_tandem_lindley: W0 must hold one waiting time per station");
    for (std::size_t k = 0; k < K; ++k)
        if (w0[k] < zero)
            throw InputError("qsys_tandem_lindley: W0 must hold nonnegative waiting times");

    TandemLindleyResult<T> r;
    r.analyzer = "qsys_tandem_lindley";
    r.W = Matrix<T>(N, K, zero);
    // the last customer has no successor, so its gap row is undefined
    r.G = Matrix<T>(N, K, num_traits<T>::from_double(std::numeric_limits<double>::quiet_NaN()));
    for (std::size_t k = 0; k < K; ++k) r.W(0, k) = w0[k];

    for (std::size_t n = 0; n + 1 < N; ++n) {
        T gap = A[n];
        for (std::size_t k = 0; k < K; ++k) {
            r.G(n, k) = gap;
            const T adv = r.W(n, k) + S(n, k);
            r.W(n + 1, k) = (adv > gap) ? T(adv - gap) : zero;
            const T residual = (gap > adv) ? T(gap - adv) : zero;
            gap = residual + S(n + 1, k);
        }
    }

    r.T_ = Matrix<T>(N, K, zero);
    for (std::size_t n = 0; n < N; ++n)
        for (std::size_t k = 0; k < K; ++k) r.T_(n, k) = r.W(n, k) + S(n, k);

    r.departure = Matrix<T>(N, K, zero);
    T epoch = zero;
    for (std::size_t n = 0; n < N; ++n) {
        if (n > 0) epoch += A[n - 1];
        r.departure(n, 0) = epoch + r.T_(n, 0);
        for (std::size_t k = 1; k < K; ++k)
            r.departure(n, k) = r.departure(n, k - 1) + r.T_(n, k);
    }
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_TANDEM_LINDLEY_H
