/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_PHMC_H
#define LINE_API_QSYS_QSYS_PHMC_H

/**
 * Exact PH/M/c by Neuts' matrix-geometric method.
 *
 * Templated port of matlab/src/api/qsys/qsys_phmc.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_phmc.java.
 *
 * The level is the number in system and the phase is the arrival PH phase.
 * Above level c the QBD is level independent with
 *
 *   A0 = (-T e) alpha,   A1 = T - c mu I,   A2 = c mu I,
 *
 * and pi_n = pi_c R^(n-c) for n >= c, R the minimal solution of
 * R^2 A2 + R A1 + A0 = 0 reached by the fixed point R <- -A0 (A1 + R A2)^-1.
 * The boundary vectors pi_0..pi_c come from the level balance equations with
 * the last one replaced by the normalization
 *
 *   sum_{n<c} pi_n e + pi_c (I-R)^-1 e = 1,
 *
 * after which Lq = pi_c R (I-R)^-2 e and L = sum_{n<c} n pi_n e +
 * pi_c (c (I-R)^-1 + R (I-R)^-2) e.
 *
 * ARITHMETIC. R is a fixed point driven to a tolerance and never terminates in
 * a finite number of field operations, so the function is gated on
 * transcendental arithmetic, for the same reason qbd_R is (see qbd_r.h).
 *
 * At k = 1 with T = [-lambda] the arrival process is Poisson and every metric
 * must collapse onto M/M/c, which is the check the tests apply.
 */

#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

template <class T>
struct PhMcResult {
    T meanQueueLength;   ///< L
    T meanWaitingQueue;  ///< Lq
    T meanWaitingTime;   ///< Wq
    T meanSojournTime;   ///< W
    T utilization;       ///< rho = lambda/(c mu)
};

namespace detail {

/**
 * Determinant through the LU factorization, standing in for MATLAB's det in
 * the singularity guard of the R iteration. A zero pivot is reported as a zero
 * determinant rather than an error, which is the branch MATLAB takes.
 */
template <class T>
T matrix_det(const Matrix<T>& A) {
    Matrix<T> LU = A;
    std::vector<std::size_t> piv;
    try {
        piv = lu_factor(LU);
    } catch (const NumericError&) {
        return num_traits<T>::from_int(0);
    }
    T d = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < LU.rows(); ++i) d *= LU(i, i);
    for (std::size_t i = 0; i < piv.size(); ++i)
        if (piv[i] != i) d = -d;
    return d;
}

}  // namespace detail

/**
 * @param alpha   PH entry probability vector, length k
 * @param Tm      PH sub-generator, k x k
 * @param mu      exponential service rate of one server
 * @param c       number of servers, c >= 1
 * @param maxIter iteration budget for the R fixed point (MATLAB 50000)
 * @param tol     convergence tolerance on R (MATLAB 1e-14)
 */
template <class T>
PhMcResult<T> qsys_phmc(const std::vector<T>& alpha, const Matrix<T>& Tm, const T& mu, unsigned c,
                        unsigned maxIter, const T& tol) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_phmc requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (mu <= zero) throw InputError("qsys_phmc: service rate mu must be positive");
    if (c < 1) throw InputError("qsys_phmc: c must be a positive integer");
    const std::size_t k = Tm.rows();
    if (Tm.cols() != k || alpha.size() != k)
        throw InputError("qsys_phmc: alpha and T dimensions are inconsistent");
    const T ct = num_traits<T>::from_int(static_cast<long>(c));

    std::vector<T> t_vec(k, zero);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) t_vec[i] -= Tm(i, j);
    Matrix<T> D1(k, k);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) D1(i, j) = t_vec[i] * alpha[j];

    Matrix<T> negT(k, k);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) negT(i, j) = -Tm(i, j);
    const std::vector<T> mvec = line::solve(negT, ones<T>(k));
    T mean_ia = zero;
    for (std::size_t i = 0; i < k; ++i) mean_ia += alpha[i] * mvec[i];
    if (mean_ia <= zero) throw InputError("qsys_phmc: non-positive mean interarrival time");
    const T lambda = one / mean_ia;
    const T rho = lambda / (ct * mu);
    if (rho >= one) throw InputError("qsys_phmc: load rho must be strictly less than 1");

    // R by successive substitution: R <- -A0 (A1 + R A2)^-1.
    const Matrix<T>& A0 = D1;
    Matrix<T> A1(k, k), A2(k, k, zero);
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < k; ++j) A1(i, j) = Tm(i, j);
        A1(i, i) -= ct * mu;
        A2(i, i) = ct * mu;
    }
    Matrix<T> R(k, k, zero);
    const T det_floor = T(num_traits<T>::from_double(1e-30));
    for (unsigned it = 0; it < maxIter; ++it) {
        Matrix<T> M = A1;
        const Matrix<T> RA2 = matmul(R, A2);
        for (std::size_t i = 0; i < k; ++i)
            for (std::size_t j = 0; j < k; ++j) M(i, j) += RA2(i, j);
        if (num_abs(detail::matrix_det(M)) < det_floor) break;
        Matrix<T> Rn = matmul(A0, inverse(M));
        for (std::size_t i = 0; i < k; ++i)
            for (std::size_t j = 0; j < k; ++j) Rn(i, j) = -Rn(i, j);
        T gap = zero;
        for (std::size_t i = 0; i < k; ++i)
            for (std::size_t j = 0; j < k; ++j) {
                const T d = num_abs(T(Rn(i, j) - R(i, j)));
                if (d > gap) gap = d;
            }
        R = Rn;
        if (gap < tol) break;
    }

    // Boundary system on the stacked unknowns [pi_0 ... pi_c].
    const std::size_t nvar = (c + 1) * k;
    Matrix<T> M(nvar, nvar, zero);
    // Level 0: pi_0 T + pi_1 (mu I) = 0.
    for (std::size_t j = 0; j < k; ++j) {
        for (std::size_t i = 0; i < k; ++i) M(j, i) += Tm(i, j);
        M(j, k + j) += mu;
    }
    // Levels 1..c-1: pi_{n-1} D1 + pi_n (T - n mu I) + pi_{n+1} ((n+1) mu I) = 0.
    for (unsigned n = 1; n < c; ++n) {
        const T nt = num_traits<T>::from_int(static_cast<long>(n));
        for (std::size_t j = 0; j < k; ++j) {
            const std::size_t row = n * k + j;
            for (std::size_t i = 0; i < k; ++i) {
                M(row, (n - 1) * k + i) += D1(i, j);
                M(row, n * k + i) += Tm(i, j);
            }
            M(row, n * k + j) -= nt * mu;
            M(row, (n + 1) * k + j) += (nt + one) * mu;
        }
    }
    // Level c: pi_{c-1} D1 + pi_c (T - c mu I + R (c mu I)) = 0.
    const Matrix<T> RA2c = matmul(R, A2);
    for (std::size_t j = 0; j < k; ++j) {
        const std::size_t row = c * k + j;
        for (std::size_t i = 0; i < k; ++i) {
            M(row, (c - 1) * k + i) += D1(i, j);
            M(row, c * k + i) += Tm(i, j) + RA2c(i, j);
        }
        M(row, c * k + j) -= ct * mu;
    }
    // Replace the last equation with the normalization.
    Matrix<T> IR(k, k);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) IR(i, j) = (i == j ? one : zero) - R(i, j);
    const std::vector<T> sum_geom = line::solve(IR, ones<T>(k));
    for (std::size_t col = 0; col < nvar; ++col) M(nvar - 1, col) = zero;
    for (unsigned n = 0; n < c; ++n)
        for (std::size_t i = 0; i < k; ++i) M(nvar - 1, n * k + i) = one;
    for (std::size_t i = 0; i < k; ++i) M(nvar - 1, c * k + i) = sum_geom[i];
    std::vector<T> b(nvar, zero);
    b[nvar - 1] = one;
    const std::vector<T> x = line::solve(M, b);

    std::vector<T> pi_c(k);
    for (std::size_t i = 0; i < k; ++i) pi_c[i] = x[c * k + i];

    const Matrix<T> IRinv = inverse(IR);
    const Matrix<T> IRinv2 = matmul(IRinv, IRinv);
    const Matrix<T> R_IRinv2 = matmul(R, IRinv2);
    const std::vector<T> e = ones<T>(k);
    const std::vector<T> v1 = mulvec(R_IRinv2, e);
    T Lq = zero;
    for (std::size_t i = 0; i < k; ++i) Lq += pi_c[i] * v1[i];
    Matrix<T> bulk(k, k);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) bulk(i, j) = ct * IRinv(i, j) + R_IRinv2(i, j);
    const std::vector<T> v2 = mulvec(bulk, e);
    T L = zero;
    for (std::size_t i = 0; i < k; ++i) L += pi_c[i] * v2[i];
    for (unsigned n = 0; n < c; ++n) {
        T s = zero;
        for (std::size_t i = 0; i < k; ++i) s += x[n * k + i];
        L += num_traits<T>::from_int(static_cast<long>(n)) * s;
    }

    PhMcResult<T> r;
    r.meanQueueLength = L;
    r.meanWaitingQueue = Lq;
    r.meanWaitingTime = Lq / lambda;
    r.meanSojournTime = r.meanWaitingTime + one / mu;
    r.utilization = rho;
    return r;
}

/** qsys_phmc with the MATLAB defaults, 50000 iterations and tolerance 1e-14. */
template <class T>
PhMcResult<T> qsys_phmc(const std::vector<T>& alpha, const Matrix<T>& Tm, const T& mu,
                        unsigned c) {
    return qsys_phmc(alpha, Tm, mu, c, 50000u, T(num_traits<T>::from_double(1e-14)));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_PHMC_H
