/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_BMAPM1_H
#define LINE_API_QSYS_QSYS_BMAPM1_H

/**
 * BMAP/M/1 by the matrix-analytic (M/G/1-type) method.
 *
 * Port of `matlab/src/api/qsys/qsys_bmapm1.m`. What makes this function
 * unusual, and what `getMAMResult` exists to expose, is that it returns the
 * INTERMEDIATE OBJECTS rather than only the mean measures: the randomized
 * blocks, the matrix G, the drift, the measured decay rate and the level
 * probabilities. Mean values alone hide the objects the method is built on, so
 * a matrix-analytic result cannot otherwise be checked against a published
 * derivation.
 *
 * THE RANDOMIZATION. The continuous-time chain is uniformized by q, chosen as
 * `max_i(-D0(i,i)) + mu` unless supplied, into a discrete-time M/G/1-type chain
 * with blocks A0 = (mu/q)I (a service completion, level down), A1 = (1/q)(D0 -
 * mu I) + I (level unchanged), Bk = (1/q)D_k (level up by k) and the boundary
 * local block B0 = (1/q)D0 + I, which differs from A1 because no service can
 * complete at level zero. A q that does not dominate the outflow would give the
 * randomized chain negative entries, and is refused by name rather than
 * silently clamped.
 *
 * THE DECAY RATE IS MEASURED, NOT DERIVED, and that is deliberate in the
 * reference: it is read as the ratio pi_(n+1)/pi_n at a level where the mass is
 * still numerically meaningful (half way up the usable range), rather than
 * taken from a spectral convention that would have to pick a branch. That is
 * also why this function needs NO eigen-decomposition and is therefore
 * instantiable at `Real`, not double-only. Every returned quantity is a linear
 * solve, a functional iteration or a ratio.
 *
 * ARITHMETIC. Gated on `num_traits<T>::has_transcendental` because the G
 * iteration and the level truncation both terminate on a tolerance. The
 * algebra itself is rational.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/dtmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Everything `qsys_bmapm1` returns, mirroring the MATLAB result struct. */
template <class T>
struct BmapM1Result {
    std::vector<T> theta;   ///< stationary vector of the BMAP phase process
    T lambda;               ///< mean arrival rate, theta (sum_k k D_k) e
    T rho;                  ///< offered load lambda/mu
    T q;                    ///< uniformization constant actually used
    Matrix<T> A0, A1, B0;   ///< randomized blocks: down, local, boundary local
    std::vector<Matrix<T>> Bk;  ///< level up by k, k = 1..K
    Matrix<T> A;            ///< A0 + A1 + sum_k Bk, the phase process
    std::vector<T> alpha;   ///< stationary vector of A
    Matrix<T> G;            ///< minimal non-negative solution of the M/G/1-type equation
    T drift;                ///< stable iff strictly negative
    double decayRate;       ///< measured pi_(n+1)/pi_n; NaN when unmeasurable
    Matrix<T> levelProb;    ///< level probabilities, row n = pi_n
    T pi0;                  ///< probability the system is empty
    T meanQueueLength;
    T utilization;
    T throughput;
    std::size_t truncLevel;
    double truncError;
    bool gConverged = true;
};

namespace bmapm1_detail {

/**
 * Stationary distribution of the level-truncated CTMC.
 *
 * The reference replaces the last column of Q by ones and solves b/Q with
 * b = e_last, which is the standard normalized stationary solve; reproduced
 * here on the transpose because `line::solve` takes a column right-hand side.
 */
template <class T>
Matrix<T> solve_levels(const std::vector<Matrix<T>>& D, const T& mu, std::size_t V, std::size_t K,
                       std::size_t levelMax) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t dim = (levelMax + 1) * V;
    Matrix<T> Q(dim, dim, zero);
    for (std::size_t n = 0; n <= levelMax; ++n) {
        const std::size_t base = n * V;
        // D0 on every level
        for (std::size_t i = 0; i < V; ++i)
            for (std::size_t j = 0; j < V; ++j) Q(base + i, base + j) += D[0](i, j);
        // service: level n -> n-1 for n >= 1
        if (n >= 1)
            for (std::size_t i = 0; i < V; ++i) Q(base + i, base - V + i) += mu;
        // batch arrivals: level n -> n+k
        for (std::size_t k = 1; k <= K; ++k) {
            if (n + k > levelMax) continue;
            for (std::size_t i = 0; i < V; ++i)
                for (std::size_t j = 0; j < V; ++j)
                    Q(base + i, (n + k) * V + j) += D[k](i, j);
        }
    }
    // Row-sum correction: the reference subtracts the row sums from the
    // diagonal AFTER assembling, so D0's own negative diagonal is included in
    // the sum and the result is a proper generator of the truncated chain.
    for (std::size_t i = 0; i < dim; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < dim; ++j) s += Q(i, j);
        Q(i, i) -= s;
    }
    // Replace the last column by ones and solve pi Q = e_last^T.
    for (std::size_t i = 0; i < dim; ++i) Q(i, dim - 1) = one;
    Matrix<T> Qt(dim, dim, zero);
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim; ++j) Qt(i, j) = Q(j, i);
    std::vector<T> b(dim, zero);
    b[dim - 1] = one;
    const std::vector<T> pi = solve(Qt, b);

    Matrix<T> levelProb(levelMax + 1, V, zero);
    T tot = zero;
    for (std::size_t n = 0; n <= levelMax; ++n)
        for (std::size_t i = 0; i < V; ++i) {
            const T v = pi[n * V + i];
            levelProb(n, i) = (v < zero) ? zero : v;
            tot += levelProb(n, i);
        }
    if (tot > zero)
        for (std::size_t n = 0; n <= levelMax; ++n)
            for (std::size_t i = 0; i < V; ++i) levelProb(n, i) /= tot;
    return levelProb;
}

/** Relative contribution the truncated tail would add to the mean level. */
template <class T>
double level_tail_error(const Matrix<T>& levelProb) {
    const std::size_t L = levelProb.rows();
    double meanLevel = 0.0, last = 0.0;
    for (std::size_t n = 0; n < L; ++n) {
        double m = 0.0;
        for (std::size_t i = 0; i < levelProb.cols(); ++i)
            m += num_traits<T>::to_double(levelProb(n, i));
        meanLevel += static_cast<double>(n) * m;
        if (n + 1 == L) last = m;
    }
    const double denom = std::max(meanLevel, std::numeric_limits<double>::min());
    return static_cast<double>(L - 1) * last / denom;
}

}  // namespace bmapm1_detail

/**
 * @param D  the BMAP {D0, D1, ..., DK}; D0 carries hidden transitions, Dk the
 *           transitions releasing a batch of k
 * @param mu exponential service rate
 * @param qParam uniformization constant; <= 0 selects the reference's default
 * @param maxLevelParam explicit level truncation; 0 selects the adaptive search
 * @param maxIter iteration cap for the G matrix
 * @param tol convergence tolerance for the G matrix
 * @param tailTol relative truncation target for the level distribution
 */
template <class T>
BmapM1Result<T> qsys_bmapm1(const std::vector<Matrix<T>>& D, const T& mu, const T& qParam,
                            std::size_t maxLevelParam, unsigned maxIter, const T& tol,
                            double tailTol) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_bmapm1 runs tolerance-terminated iterations (G, and the level truncation)");
    using lang::GlobalConstants;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (D.size() < 2)
        throw InputError(
            "qsys_bmapm1: the BMAP must be given as {D0,D1,...,DK} with at least D0 and D1");
    const std::size_t V = D[0].rows();
    for (std::size_t k = 0; k < D.size(); ++k) {
        if (D[k].rows() != V || D[k].cols() != V)
            throw InputError("qsys_bmapm1: BMAP matrix D{" + std::to_string(k) +
                             "} is not conformable");
        if (k >= 1)
            for (std::size_t i = 0; i < V; ++i)
                for (std::size_t j = 0; j < V; ++j)
                    if (num_traits<T>::to_double(D[k](i, j)) < -GlobalConstants::FineTol)
                        throw InputError("qsys_bmapm1: BMAP arrival matrix D{" +
                                         std::to_string(k) + "} must be non-negative");
    }
    Matrix<T> Dsum(V, V, zero);
    for (const Matrix<T>& Dk : D)
        for (std::size_t i = 0; i < V; ++i)
            for (std::size_t j = 0; j < V; ++j) Dsum(i, j) += Dk(i, j);
    for (std::size_t i = 0; i < V; ++i) {
        T s = zero;
        for (std::size_t j = 0; j < V; ++j) s += Dsum(i, j);
        if (std::fabs(num_traits<T>::to_double(s)) > std::sqrt(GlobalConstants::FineTol))
            throw InputError(
                "qsys_bmapm1: BMAP matrices are inconsistent, sum_k D_k must have zero row sums");
    }
    if (!(num_traits<T>::to_double(mu) > 0.0) || !std::isfinite(num_traits<T>::to_double(mu)))
        throw InputError("qsys_bmapm1: the service rate mu must be a finite positive scalar");

    const std::size_t K = D.size() - 1;
    BmapM1Result<T> r;

    r.theta = mc::ctmc_solve(Dsum);
    Matrix<T> sumKDk(V, V, zero);
    for (std::size_t k = 1; k <= K; ++k)
        for (std::size_t i = 0; i < V; ++i)
            for (std::size_t j = 0; j < V; ++j)
                sumKDk(i, j) += num_traits<T>::from_int(static_cast<int>(k)) * D[k](i, j);
    r.lambda = zero;
    {
        const std::vector<T> t = vecmul(r.theta, sumKDk);
        for (const T& v : t) r.lambda += v;
    }
    r.rho = T(r.lambda / mu);

    T outflow = T(-D[0](0, 0));
    for (std::size_t i = 1; i < V; ++i)
        if (T(-D[0](i, i)) > outflow) outflow = T(-D[0](i, i));
    const T qmin = T(outflow + mu);
    r.q = (num_traits<T>::to_double(qParam) > 0.0) ? qParam : qmin;
    if (num_traits<T>::to_double(r.q) <
        num_traits<T>::to_double(qmin) - GlobalConstants::FineTol)
        throw InputError(
            "qsys_bmapm1: the uniformization constant does not dominate the total outflow rate; "
            "the randomized chain would have negative entries");

    r.A0 = Matrix<T>(V, V, zero);
    r.A1 = Matrix<T>(V, V, zero);
    r.B0 = Matrix<T>(V, V, zero);
    for (std::size_t i = 0; i < V; ++i) {
        r.A0(i, i) = T(mu / r.q);
        for (std::size_t j = 0; j < V; ++j) {
            r.A1(i, j) = T(D[0](i, j) / r.q);
            r.B0(i, j) = T(D[0](i, j) / r.q);
        }
        r.A1(i, i) = T(r.A1(i, i) - mu / r.q + one);
        r.B0(i, i) = T(r.B0(i, i) + one);
    }
    r.Bk.resize(K);
    for (std::size_t k = 1; k <= K; ++k) {
        r.Bk[k - 1] = Matrix<T>(V, V, zero);
        for (std::size_t i = 0; i < V; ++i)
            for (std::size_t j = 0; j < V; ++j) r.Bk[k - 1](i, j) = T(D[k](i, j) / r.q);
    }
    r.A = Matrix<T>(V, V, zero);
    for (std::size_t i = 0; i < V; ++i)
        for (std::size_t j = 0; j < V; ++j) {
            r.A(i, j) = r.A0(i, j) + r.A1(i, j);
            for (std::size_t k = 0; k < K; ++k) r.A(i, j) += r.Bk[k](i, j);
        }
    r.alpha = mc::dtmc_solve(r.A);

    // G = A0 + A1 G + sum_k Bk G^(k+1), by functional iteration from zero.
    r.G = Matrix<T>(V, V, zero);
    r.gConverged = false;
    for (unsigned it = 0; it < maxIter; ++it) {
        Matrix<T> Gpow = r.G;
        Matrix<T> Gnew = matmul(r.A1, r.G);
        for (std::size_t i = 0; i < V; ++i)
            for (std::size_t j = 0; j < V; ++j) Gnew(i, j) += r.A0(i, j);
        for (std::size_t k = 0; k < K; ++k) {
            Gpow = matmul(Gpow, r.G);
            const Matrix<T> add = matmul(r.Bk[k], Gpow);
            for (std::size_t i = 0; i < V; ++i)
                for (std::size_t j = 0; j < V; ++j) Gnew(i, j) += add(i, j);
        }
        double diff = 0.0;
        for (std::size_t i = 0; i < V; ++i)
            for (std::size_t j = 0; j < V; ++j)
                diff = std::max(diff, std::fabs(num_traits<T>::to_double(Gnew(i, j)) -
                                                num_traits<T>::to_double(r.G(i, j))));
        r.G = Gnew;
        if (diff < num_traits<T>::to_double(tol)) {
            r.gConverged = true;
            break;
        }
    }

    Matrix<T> upDrift(V, V, zero);
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t i = 0; i < V; ++i)
            for (std::size_t j = 0; j < V; ++j)
                upDrift(i, j) +=
                    num_traits<T>::from_int(static_cast<int>(k + 1)) * r.Bk[k](i, j);
    r.drift = zero;
    {
        const std::vector<T> up = vecmul(r.alpha, upDrift);
        const std::vector<T> dn = vecmul(r.alpha, r.A0);
        for (std::size_t i = 0; i < V; ++i) r.drift += up[i] - dn[i];
    }

    // Level probabilities, adaptively refined until the tail is negligible.
    std::size_t levelMax;
    if (maxLevelParam > 0) {
        levelMax = maxLevelParam;
        r.levelProb = bmapm1_detail::solve_levels(D, mu, V, K, levelMax);
        r.truncError = bmapm1_detail::level_tail_error(r.levelProb);
    } else {
        const double rd = num_traits<T>::to_double(r.rho);
        const double slack = std::max(1.0 - std::min(rd, 0.999),
                                      std::numeric_limits<double>::epsilon());
        levelMax = std::max<std::size_t>(50, static_cast<std::size_t>(std::ceil(20.0 / slack)));
        for (;;) {
            r.levelProb = bmapm1_detail::solve_levels(D, mu, V, K, levelMax);
            r.truncError = bmapm1_detail::level_tail_error(r.levelProb);
            if (r.truncError <= tailTol || (2 * levelMax + 1) * V > 200000) break;
            levelMax *= 2;
        }
    }

    std::vector<double> levelMass(r.levelProb.rows(), 0.0);
    for (std::size_t n = 0; n < r.levelProb.rows(); ++n)
        for (std::size_t i = 0; i < V; ++i)
            levelMass[n] += num_traits<T>::to_double(r.levelProb(n, i));
    // Read the ratio where the mass is still meaningful, not at the boundary.
    // THE INDEX IS THE REFERENCE'S AND IT IS 1-BASED: `usable` is
    // find(levelMass > 1e-12, 1, 'last'), `ref = max(2, floor(usable/2))`, and
    // the ratio is levelMass(ref+1)/levelMass(ref). Reading it one level higher
    // samples the ratio before it has settled: on the Bolch fixture that moves
    // the answer from 0.413155865252739 to 0.413243742742746, a 2.1e-4 relative
    // shift that no other returned quantity shows.
    std::size_t usable1 = 0;  // 1-based, 0 = none
    for (std::size_t n = 0; n < levelMass.size(); ++n)
        if (levelMass[n] > 1e-12) usable1 = n + 1;
    if (usable1 < 3) {
        r.decayRate = std::numeric_limits<double>::quiet_NaN();
    } else {
        const std::size_t ref0 = std::max<std::size_t>(2, usable1 / 2) - 1;
        r.decayRate = (ref0 + 1 < levelMass.size() && levelMass[ref0] > 0.0)
                          ? levelMass[ref0 + 1] / levelMass[ref0]
                          : std::numeric_limits<double>::quiet_NaN();
    }

    r.meanQueueLength = zero;
    for (std::size_t n = 0; n < levelMass.size(); ++n)
        r.meanQueueLength +=
            num_traits<T>::from_int(static_cast<int>(n)) * num_traits<T>::from_double(levelMass[n]);
    r.pi0 = num_traits<T>::from_double(levelMass[0]);
    r.utilization = r.rho;
    r.throughput = r.lambda;
    r.truncLevel = r.levelProb.rows() - 1;
    return r;
}

/** The reference's defaults: adaptive truncation, 10000 iterations, tol 1e-12. */
template <class T>
BmapM1Result<T> qsys_bmapm1(const std::vector<Matrix<T>>& D, const T& mu) {
    return qsys_bmapm1(D, mu, num_traits<T>::from_int(0), static_cast<std::size_t>(0), 10000u,
                       T(num_traits<T>::from_double(1e-12)), 1e-10);
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_BMAPM1_H
