/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MG1K_LOSS_H
#define LINE_API_QSYS_QSYS_MG1K_LOSS_H

/**
 * Exact M/G/1/K loss probability, via the chain embedded at service-start
 * epochs.
 *
 * Templated port of matlab/src/api/qsys/qsys_mg1k_loss.m, cross-checked
 * against jar/src/main/java/jline/api/qsys/Qsys_mg1k_loss.java.
 *
 * The state is the number waiting immediately after a service start,
 * q in {0,...,K-2}, and a_j is the probability of j Poisson arrivals during
 * one service:
 *
 *   q = 0 : both a_0 and a_1 lead to q' = 0, because after an empty departure
 *           the next service starts with the next arrival; j >= 2 gives
 *           q' = j-1;
 *   q >= 1: q' = q-1+j, with arrivals past the free capacity lost and
 *           aggregated in the last column.
 *
 * The loss probability then follows from renewal-reward,
 *
 *   E[cycle] = E[S] + sigma_0 a_0/lambda,   P_loss = 1 - 1/(rho + sigma_0 a_0),
 *
 * with sigma the stationary vector at service-start epochs.
 *
 * The service law is supplied as its density. Both the mean service time and
 * the a_j are obtained by adaptive quadrature on [0, 1e4/lambda], exactly the
 * truncation and the tolerances MATLAB's integral() is given (RelTol 1e-6,
 * AbsTol 1e-10), and the a_j series is stopped at 1e-12 or after 1000 terms as
 * in the reference. For an exponential density the result must reproduce
 * qsys_mm1k_loss, and that identity is the sharpest available check on the
 * embedded chain.
 *
 * ARITHMETIC. The quadrature and the exp in the Poisson weights make this
 * transcendental. The chain itself -- row normalization and the stationary
 * solve -- is finite field arithmetic, but it is fed by the quadrature, so the
 * gate applies to the whole function.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/dtmc_solve.h"
#include "line/api/qsys/qsys_quadrature.h"
#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

template <class T>
struct Mg1kLossResult {
    T lossProbability;
    T utilization;  ///< rho = lambda/mu with mu the reciprocal mean service time
};

namespace detail {

/**
 * Row-normalize P and repair the diagonal, MATLAB's dtmc_makestochastic. A row
 * that sums to zero is replaced by a self-loop.
 */
template <class T>
void dtmc_makestochastic_inplace(Matrix<T>& P) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < P.rows(); ++i) {
        T s = zero;
        for (std::size_t j = 0; j < P.cols(); ++j) s += P(i, j);
        if (s > zero) {
            for (std::size_t j = 0; j < P.cols(); ++j) P(i, j) /= s;
            T off = zero;
            for (std::size_t j = 0; j < P.cols(); ++j)
                if (j != i) off += P(i, j);
            T d = one - off;
            if (d < zero) d = zero;
            if (d > one) d = one;
            P(i, i) = d;
        } else {
            for (std::size_t j = 0; j < P.cols(); ++j) P(i, j) = zero;
            P(i, i) = one;
        }
    }
}

}  // namespace detail

/**
 * @param lambda  Poisson arrival rate
 * @param density service-time density, callable as density(t) -> T
 * @param K       system capacity, jobs in service included, K >= 2
 */
template <class T, class Density>
Mg1kLossResult<T> qsys_mg1k_loss(const T& lambda, Density&& density, unsigned K) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mg1k_loss requires transcendental arithmetic");
    if (K < 2) throw InputError("qsys_mg1k_loss: K must be at least 2");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T reltol = T(num_traits<T>::from_double(1e-6));
    const T abstol = T(num_traits<T>::from_double(1e-10));
    const T tmax = T(num_traits<T>::from_double(1e4)) / lambda;

    const T meanS = detail::num_integral<T>([&](const T& t) { return t * density(t); }, zero, tmax,
                                            reltol, abstol);
    if (meanS <= zero) throw NumericError("qsys_mg1k_loss: non-positive mean service time");
    const T mu = one / meanS;

    // a_j = P(j Poisson arrivals during one service), j = 0, 1, ...
    std::vector<T> a(K - 1, zero);
    const T stop = T(num_traits<T>::from_double(1e-12));
    for (unsigned j = 0; j <= 1000u; ++j) {
        const T aj = detail::num_integral<T>(
                         [&](const T& t) {
                             return detail::num_exp(T(-lambda * t)) *
                                    num_pow_int(T(lambda * t), j) * density(t);
                         },
                         zero, tmax, reltol, abstol) /
                     num_factorial<T>(j);
        if (j < a.size())
            a[j] = aj;
        else
            a.push_back(aj);
        if (aj < stop) break;
    }

    // Embedded chain at service-start epochs, states q = 0..K-2.
    const std::size_t n = K - 1;
    Matrix<T> P(n, n, zero);
    const std::size_t last = n - 1;
    P(0, 0) = a[0] + (a.size() > 1 ? a[1] : zero);
    for (std::size_t i = 1; i + 2 <= n; ++i) P(0, i) = i + 1 < a.size() ? a[i + 1] : zero;
    if (n >= 2) {
        for (std::size_t i = 0; i + 1 < n; ++i) P(1, i) = i < a.size() ? a[i] : zero;
    }
    for (std::size_t r = 2; r < n; ++r)
        for (std::size_t col = r - 1; col + 1 < n; ++col) {
            const std::size_t j = col - r + 1;
            P(r, col) = j < a.size() ? a[j] : zero;
        }
    // Last column absorbs everything the free capacity cannot take.
    for (std::size_t r = 0; r < n; ++r) {
        T s = zero;
        for (std::size_t col = 0; col + 1 < n; ++col) s += P(r, col);
        P(r, last) = one - s;
    }
    detail::dtmc_makestochastic_inplace(P);
    const std::vector<T> sigma = mc::dtmc_solve(P);

    Mg1kLossResult<T> r;
    r.utilization = lambda / mu;
    const T carried = sigma[0] * a[0] + r.utilization;
    if (carried == zero) throw NumericError("qsys_mg1k_loss: degenerate renewal cycle");
    r.lossProbability = one - one / carried;
    return r;
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MG1K_LOSS_H
