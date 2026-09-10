/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MOMLIN_H
#define LINE_API_PFQN_MOMLIN_H

/**
 * Moment linearizer: approximate first and second queue-length moments of a
 * large closed product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_momlin.m.
 *
 * Means come from the Schweitzer-Bard AMVA fixed point. Second moments use the
 * exact product-form identity
 *
 *   Cov[n_{i,r}, n_{j,s}] = D_{j,s} dQ_{i,r} / dD_{j,s},
 *
 * with the demand derivatives obtained by differentiating the AMVA fixed point
 * itself, which gives a LINEAR fixed point in dQ that is iterated to
 * convergence per parameter (j, s0). Both moments therefore carry the
 * Schweitzer-Bard error and are exact only where Schweitzer-Bard is; for exact
 * moments on tractable models use the pfqn_sens family, which differentiates
 * the exact MVA or CoMoM recursion instead.
 *
 * Arithmetic: no transcendental is used anywhere, so the routine is left
 * UNGATED and will instantiate at Rational. That is of limited practical use:
 * the fixed point converges to its limit only asymptotically, so at Rational
 * the iterates are exact but their numerators and denominators grow with every
 * sweep and the tolerance test is met only after the same number of sweeps as
 * in floating point, at far greater cost. Use double or Real for this one; the
 * exact instantiation exists so that no call site is refused, not because it
 * buys accuracy the algorithm does not have.
 *
 * REFERENCE DEFECTS: none found. The reference reassigns R (the class count)
 * to the residence-time output on its last line, which is legal MATLAB and is
 * correct because every use of the class count precedes it; the port keeps the
 * two separate.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_momlin, mirroring [Q, X, U, R, QVar, QCov, dQ]. */
template <class T>
struct MomlinResult {
    Matrix<T> Q;        ///< (M x R) mean queue length
    std::vector<T> X;   ///< (R) per-class throughput
    Matrix<T> U;        ///< (M x R) utilization
    Matrix<T> R;        ///< (M x R) residence time
    Matrix<T> QVar;     ///< (M x R) queue-length variance
    /// (M x R x M x R) covariance and demand-derivative tensors, flattened as
    /// ((i*R + r)*M + j)*R + s so that a caller can index them without a
    /// four-dimensional container.
    std::vector<T> QCov;
    std::vector<T> dQ;
    std::size_t M, Rc;

    const T& cov(std::size_t i, std::size_t r, std::size_t j, std::size_t s) const {
        return QCov[((i * Rc + r) * M + j) * Rc + s];
    }
    const T& dq(std::size_t i, std::size_t r, std::size_t j, std::size_t s) const {
        return dQ[((i * Rc + r) * M + j) * Rc + s];
    }
};

/**
 * @param L       (M x R) demand matrix
 * @param N       (R) closed populations
 * @param Z       (R) think times
 * @param tol     convergence tolerance on the fixed points
 * @param maxiter iteration cap
 */
template <class T>
MomlinResult<T> pfqn_momlin(const Matrix<T>& L, const std::vector<int>& N,
                            const std::vector<T>& Z, const T& tol, int maxiter) {
    const std::size_t M = L.rows(), R = L.cols();
    if (M == 0 || R == 0) throw InputError("pfqn_momlin: empty demand matrix");
    if (N.size() != R) throw InputError("pfqn_momlin: L and N disagree on the class count");
    if (Z.size() != R) throw InputError("pfqn_momlin: L and Z disagree on the class count");
    if (maxiter < 1) throw InputError("pfqn_momlin: maxiter must be at least one");
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] < 0) throw InputError("pfqn_momlin: pfqn_momlin supports closed classes only");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    // Schweitzer population-scaling coefficients c(r,s) = (N_s - delta_rs)/N_s.
    Matrix<T> c(R, R, one);
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t s = 0; s < R; ++s) {
            if (N[s] > 0)
                c(r, s) = num_traits<T>::from_int(N[s] - (r == s ? 1 : 0)) /
                          num_traits<T>::from_int(N[s]);
            else
                c(r, s) = zero;
        }

    // ---- Schweitzer-Bard fixed point for the means --------------------------------
    Matrix<T> Q(M, R, zero), Rm(M, R, zero);
    std::vector<T> X(R, zero);
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] > 0)
            for (std::size_t i = 0; i < M; ++i)
                Q(i, r) = num_traits<T>::from_int(N[r]) / num_traits<T>::from_int(static_cast<long>(M));
    for (int it = 0; it < maxiter; ++it) {
        const Matrix<T> Qold = Q;
        for (std::size_t r = 0; r < R; ++r) {
            if (N[r] == 0) {
                X[r] = zero;
                for (std::size_t i = 0; i < M; ++i) Rm(i, r) = zero;
                continue;
            }
            T Rtot = zero;
            for (std::size_t i = 0; i < M; ++i) {
                T acc = zero;
                for (std::size_t s = 0; s < R; ++s) acc += c(r, s) * Q(i, s);
                Rm(i, r) = L(i, r) * (one + acc);
                Rtot += Rm(i, r);
            }
            const T den = Z[r] + Rtot;
            if (den <= zero) throw NumericError("pfqn_momlin: degenerate residence time");
            X[r] = num_traits<T>::from_int(N[r]) / den;
            for (std::size_t i = 0; i < M; ++i) Q(i, r) = X[r] * Rm(i, r);
        }
        T mx = zero;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                const T d = num_abs(T(Q(i, r) - Qold(i, r)));
                if (d > mx) mx = d;
            }
        if (mx < tol) break;
    }

    MomlinResult<T> res;
    res.M = M;
    res.Rc = R;
    res.Q = Q;
    res.X = X;
    res.R = Rm;
    res.U = Matrix<T>(M, R, zero);
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t i = 0; i < M; ++i) res.U(i, r) = X[r] * L(i, r);

    // ---- analytic linearization of the fixed point ---------------------------------
    res.dQ.assign(M * R * M * R, zero);
    for (std::size_t j = 0; j < M; ++j) {
        for (std::size_t s0 = 0; s0 < R; ++s0) {
            if (N[s0] == 0) continue;  // empty class: derivative zero
            Matrix<T> dq(M, R, zero);
            for (int it = 0; it < maxiter; ++it) {
                const Matrix<T> dqold = dq;
                for (std::size_t r = 0; r < R; ++r) {
                    if (N[r] == 0) continue;
                    std::vector<T> dR(M, zero);
                    T dRsum = zero;
                    for (std::size_t i = 0; i < M; ++i) {
                        T acc = zero, dacc = zero;
                        for (std::size_t s = 0; s < R; ++s) {
                            acc += c(r, s) * Q(i, s);
                            dacc += c(r, s) * dq(i, s);
                        }
                        dR[i] = L(i, r) * dacc;
                        if (i == j && r == s0) dR[i] += one + acc;
                        dRsum += dR[i];
                    }
                    const T dXr = -(X[r] * X[r] / num_traits<T>::from_int(N[r])) * dRsum;
                    for (std::size_t i = 0; i < M; ++i) dq(i, r) = dXr * Rm(i, r) + X[r] * dR[i];
                }
                T mx = zero;
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r = 0; r < R; ++r) {
                        const T d = num_abs(T(dq(i, r) - dqold(i, r)));
                        if (d > mx) mx = d;
                    }
                if (mx < tol) break;
            }
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r = 0; r < R; ++r)
                    res.dQ[((i * R + r) * M + j) * R + s0] = dq(i, r);
        }
    }

    // ---- second moments via the product-form covariance identity -------------------
    res.QCov.assign(M * R * M * R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t j = 0; j < M; ++j)
                for (std::size_t s = 0; s < R; ++s) {
                    const std::size_t k = ((i * R + r) * M + j) * R + s;
                    res.QCov[k] = L(j, s) * res.dQ[k];
                }
    res.QVar = Matrix<T>(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) res.QVar(i, r) = res.cov(i, r, i, r);
    return res;
}

/** Overload with the reference's defaults (tol 1e-8, maxiter 1000). */
template <class T>
MomlinResult<T> pfqn_momlin(const Matrix<T>& L, const std::vector<int>& N,
                            const std::vector<T>& Z) {
    return pfqn_momlin(L, N, Z, num_traits<T>::from_double(1e-8), 1000);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MOMLIN_H
