/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MAPMC_H
#define LINE_API_QSYS_QSYS_MAPMC_H

/**
 * The MAP/M/c FCFS queue: c identical exponential servers of rate mu fed by a
 * Markovian arrival process.
 *
 * ALGORITHM, AND HOW IT DIFFERS FROM THE MATLAB REFERENCE.
 * matlab/src/api/qsys/qsys_mapmc.m obtains these quantities by calling Q-MAM's
 * Q_CT_MAP_M_C, which is not transcribed here. This port computes the SAME
 * quantities from the port's own quasi-birth-death machinery, by building the
 * QBD directly: the level is the
 * number in system and the phase is the arrival phase, so for levels n >= c
 * the process is level independent with
 *
 *     A0 = D1  (an arrival),  A1 = D0 - c mu I  (local),  A2 = c mu I  (a
 *     completion),
 *
 * while levels 0..c-1 form the boundary, where level k has downward rate
 * k mu instead of c mu. R is the minimal non-negative solution of
 * R^2 A2 + R A1 + A0 = 0, reached by cyclic reduction (qbd_fundmat called with
 * B = A2, L = A1, F = A0), and pi_{c+k} = pi_c R^k. The boundary vectors
 * pi_0..pi_c come from the level balance equations with the last one replaced
 * by the normalization
 *
 *     sum_{n<c} pi_n e + pi_c (I - R)^-1 e = 1,
 *
 * after which the number in system is
 *
 *     L = sum_{n<c} n pi_n e + pi_c ( c (I-R)^-1 + R (I-R)^-2 ) e
 *
 * and the queued jobs are Lq = pi_c R (I-R)^-2 e. This is the same
 * construction as qsys_phmc.h with the PH arrival replaced by a general MAP:
 * a PH arrival is the special case D0 = T, D1 = (-T e) alpha. Utilization is
 * lambda/(c mu), the reference's per-server convention.
 *
 * meanWaitingTime is Lq/lambda by Little's law rather than the mean of the
 * waiting-time PH representation that the reference extracts from Q-MAM; the
 * two are the same quantity and the measurements below confirm they agree.
 *
 * MEASURED AGREEMENT (MATLAB R2025a, T = double). Metric order is
 * meanQueueLength / meanWaitingTime / meanSojournTime / utilization. Every
 * relative difference below is limited by the REFERENCE, whose maxNumComp =
 * 500 level truncation shows at the 1e-9 level; this port sums the geometric
 * tail in closed form.
 *  - M/M/c collapse, qsys_mapmc(D0 = [-1.2], D1 = [1.2], mu = 1, c = 2):
 *    MATLAB 1.874999996296936 / 0.5624999998489539 / 1.562499999848954 / 0.6.
 *    The port returns the exact Erlang-C values 1.875 / 0.5625 / 1.5625 / 0.6
 *    to 1e-15, so the differences against MATLAB are 2.0e-9, 2.7e-10 and
 *    9.7e-11; the port matches the textbook formulas and MATLAB does not.
 *  - Correlated MMPP2 arrivals D0 = [-2.5 0.2; 0.1 -0.7], D1 = diag(2.3, 0.6)
 *    (lambda = 7/6), mu = 1, c = 2: MATLAB 3.121850390835039 /
 *    1.675871772757968 / 2.675871772757968 / 0.5833333333333333; port
 *    3.12185040294165 / 1.67587177394999 / 2.67587177394999 /
 *    0.583333333333333. Relative differences 3.9e-9, 7.1e-10, 4.5e-10.
 *  - Erlang-2 arrivals D0 = [-4 4; 0 -4], D1 = [0 0; 4 0] (lambda = 2),
 *    mu = 1.5, c = 2: MATLAB 2.022056024561071 / 0.3443613471500779 /
 *    1.011028013816745 / 0.6666666666666666; port 2.02205602772974 /
 *    0.344361347198204 / 1.01102801386487 / 0.666666666666667. Relative
 *    differences 1.6e-9, 1.4e-10, 4.8e-11.
 *  - Three servers, correlated arrivals, mu = 0.5, c = 3 (rho = 7/9): MATLAB
 *    9.441571472754669 / 6.092775552492106 / 8.092775552492107 /
 *    0.7777777777777777; port 9.44157149906991 / 6.09277557063135 /
 *    8.09277557063135 / 0.777777777777778. Relative differences 2.8e-9,
 *    3.0e-9, 2.2e-9.
 *
 * That the port and not the reference is the accurate one was checked without
 * MATLAB: at c = 1 the same instance is also a MAP/MAP/1 queue with an
 * exponential service MAP, and qsys_mapmap1 -- a completely different route,
 * cyclic reduction on the Kronecker-product QBD followed by the closed-form
 * factorial moment -- reproduces this function to 5e-16 on both the correlated
 * and the Erlang-2 arrival instances. The test file carries that cross-check.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental: cyclic reduction
 * drives R to a tolerance and never terminates in a finite number of field
 * operations, exactly as in qbd_r.h and qsys_phmc.h. The boundary solve, the
 * geometric-tail sums and the moment formulas that consume R are finite exact
 * matrix algebra and add no error of their own.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/qbd_r.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Return value of qsys_mapmc and qsys_mapm1, mirroring the MATLAB struct. */
template <class T>
struct MapMcResult {
    T meanQueueLength;               ///< L, number in system
    T meanWaitingTime;               ///< Wq, time in queue
    T meanSojournTime;               ///< W = Wq + 1/mu
    T utilization;                   ///< rho = lambda / (c mu), per server
    std::vector<T> queueLengthDist;  ///< P(N = n), n = 0, 1, ...
};

/**
 * MAP/M/c by the matrix-geometric solution.
 *
 * @param arrival    arrival MAP (D0, D1)
 * @param mu         exponential service rate of one server
 * @param c          number of servers, c >= 1
 * @param dist_size  how many entries of queueLengthDist to materialize
 */
template <class T>
MapMcResult<T> qsys_mapmc(const mam::Map<T>& arrival, const T& mu, unsigned c,
                          std::size_t dist_size) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapmc requires transcendental arithmetic");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (mu <= zero) throw InputError("qsys_mapmc: service rate mu must be positive");
    if (c < 1) throw InputError("qsys_mapmc: c must be a positive integer");
    if (dist_size == 0) throw InputError("qsys_mapmc: dist_size must be positive");
    const std::size_t k = arrival.D0.rows();
    if (arrival.D0.cols() != k || arrival.D1.rows() != k || arrival.D1.cols() != k)
        throw InputError("qsys_mapmc: D0 and D1 must be square and of equal order");
    const T ct = num_traits<T>::from_int(static_cast<long>(c));

    const Matrix<T>& D0 = arrival.D0;
    const Matrix<T>& D1 = arrival.D1;
    const T lambda = mam::map_lambda(arrival);
    if (lambda <= zero) throw InputError("qsys_mapmc: non-positive arrival rate");
    const T rho = lambda / (ct * mu);
    if (rho >= one) throw InputError("qsys_mapmc: load rho must be strictly less than 1");

    // cyclic-reduction R rationale: see _kb/03-api-layer.md (cpp port notes: qsys)
    const Matrix<T>& A0 = D1;
    Matrix<T> A1(k, k), A2(k, k, zero);
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < k; ++j) A1(i, j) = D0(i, j);
        A1(i, i) -= ct * mu;
        A2(i, i) = ct * mu;
    }
    const Matrix<T> R = mam::qbd_fundmat(A2, A1, A0).R;

    // Boundary system on the stacked unknowns [pi_0 ... pi_c].
    const std::size_t nvar = (c + 1) * k;
    Matrix<T> M(nvar, nvar, zero);
    // Level 0: pi_0 D0 + pi_1 (mu I) = 0.
    for (std::size_t j = 0; j < k; ++j) {
        for (std::size_t i = 0; i < k; ++i) M(j, i) += D0(i, j);
        M(j, k + j) += mu;
    }
    // Levels 1..c-1: pi_{n-1} D1 + pi_n (D0 - n mu I) + pi_{n+1} ((n+1) mu I) = 0.
    for (unsigned n = 1; n < c; ++n) {
        const T nt = num_traits<T>::from_int(static_cast<long>(n));
        for (std::size_t j = 0; j < k; ++j) {
            const std::size_t row = n * k + j;
            for (std::size_t i = 0; i < k; ++i) {
                M(row, (n - 1) * k + i) += D1(i, j);
                M(row, n * k + i) += D0(i, j);
            }
            M(row, n * k + j) -= nt * mu;
            M(row, (n + 1) * k + j) += (nt + one) * mu;
        }
    }
    // Level c: pi_{c-1} D1 + pi_c (D0 - c mu I + R (c mu I)) = 0.
    const Matrix<T> RA2c = matmul(R, A2);
    for (std::size_t j = 0; j < k; ++j) {
        const std::size_t row = c * k + j;
        for (std::size_t i = 0; i < k; ++i) {
            M(row, (c - 1) * k + i) += D1(i, j);
            M(row, c * k + i) += D0(i, j) + RA2c(i, j);
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

    MapMcResult<T> r;
    r.meanQueueLength = L;
    r.meanWaitingTime = Lq / lambda;
    r.meanSojournTime = r.meanWaitingTime + one / mu;
    r.utilization = rho;

    // Level distribution: pi_n from the boundary for n <= c, pi_c R^(n-c) above.
    r.queueLengthDist.assign(dist_size, zero);
    std::vector<T> cur(k, zero);
    for (std::size_t n = 0; n < dist_size; ++n) {
        if (n <= c) {
            for (std::size_t i = 0; i < k; ++i) cur[i] = x[n * k + i];
        } else {
            cur = vecmul(cur, R);
        }
        T s = zero;
        for (std::size_t i = 0; i < k; ++i) s += cur[i];
        r.queueLengthDist[n] = s;
    }
    return r;
}

/** qsys_mapmc with 100 materialized levels, the reference's maxNumComp scale. */
template <class T>
MapMcResult<T> qsys_mapmc(const mam::Map<T>& arrival, const T& mu, unsigned c) {
    return qsys_mapmc(arrival, mu, c, static_cast<std::size_t>(100));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MAPMC_H
