/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_QBD_MAPMAP1_H
#define LINE_API_MAM_QBD_MAPMAP1_H

/**
 * The MAP/MAP/1 queue solved as a quasi-birth-death process.
 *
 * Templated port of matlab/src/api/mam/qbd_mapmap1.m and qbd_rg.m,
 * cross-checked against jar/src/main/java/jline/api/mam/Qbd_mapmap1.java.
 *
 * The level is the number in system and the phase is the pair (arrival phase,
 * service phase), so the level blocks are
 *
 *     F    = D1^a (x) I_ns          an arrival, level up
 *     L    = D0^a (+) D0^s          no event, level unchanged
 *     B    = I_na (x) D1^s          a service completion, level down
 *     Lbar = D0^a (x) I_ns          level zero: no server is busy
 *
 * with (x) the Kronecker product and (+) the Kronecker sum. R and G come from
 * cyclic reduction (qbd_fundmat), the boundary vector from qbd_pi.
 *
 * ARITHMETIC. qbd_mapmap1 itself is gated on transcendental arithmetic because
 * it calls qbd_fundmat; see qbd_r.h for why the fixed-point iterations cannot
 * be exact. Everything that consumes R afterwards is un-gated and instantiates
 * at Rational: qbd_mapmap1_blocks assembles the level blocks with Kronecker
 * products only, and the moment formulas below are closed forms in R,
 *
 *     E[N (N-1) ... (N-m+1)] = m! pi_0 R^m (I - R)^-(m+1) e,
 *
 * a finite product of exact matrix operations. Those closed forms are also
 * more accurate than the reference: MATLAB and the JAR sum the truncated level
 * distribution k pi_k e until the accumulated mass reaches 1 - 1e-10, which
 * discards the tail mass times its (unbounded) level index, whereas the closed
 * form sums the whole geometric tail. The difference grows as the caudal
 * characteristic approaches one; see the note on QN_truncated below.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/qbd_r.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** The four level blocks of the MAP/MAP/1 QBD. */
template <class T>
struct QbdMapMap1Blocks {
    Matrix<T> B;     ///< A_-1, service completion
    Matrix<T> L;     ///< A_0, local
    Matrix<T> F;     ///< A_1, arrival
    Matrix<T> Lbar;  ///< level-zero local block
};

/**
 * Level blocks of the MAP/MAP/1 QBD from the arrival and service MAPs.
 * Kronecker products only, so exact at Rational.
 */
template <class T>
QbdMapMap1Blocks<T> qbd_mapmap1_blocks(const Map<T>& arrival, const Map<T>& service) {
    const std::size_t na = arrival.order();
    const std::size_t ns = service.order();
    QbdMapMap1Blocks<T> b;
    b.F = kron(arrival.D1, eye<T>(ns));
    b.L = krons(arrival.D0, service.D0);
    b.B = kron(eye<T>(na), service.D1);
    b.Lbar = kron(arrival.D0, eye<T>(ns));
    return b;
}

/** R, G and the level blocks of a MAP/MAP/1 queue (qbd_rg.m). */
template <class T>
struct QbdRg {
    Matrix<T> R;
    Matrix<T> G;
    Matrix<T> B, L, F, Lbar;
    Matrix<T> U;  ///< U = L + R B, the generator of the taboo process at a level
};

/**
 * R and G of the MAP/MAP/1 QBD (qbd_rg.m). Gated: qbd_fundmat is iterative.
 *
 * @param util if positive, the service MAP is first rescaled to mean
 *             util / lambda_arrival, exactly as the optional third argument of
 *             qbd_rg.m and qbd_mapmap1.m does through map_scale
 * @param arrival the arrival MAP
 * @param service_in the service MAP, rescaled to the requested utilization
 */
template <class T>
QbdRg<T> qbd_rg(const Map<T>& arrival, const Map<T>& service_in, const T& util) {
    static_assert(num_traits<T>::has_transcendental, "qbd_rg requires transcendental arithmetic");
    Map<T> service = service_in;
    if (util > num_traits<T>::from_int(0)) service = map_scale(service, T(util / map_lambda(arrival)));
    const QbdMapMap1Blocks<T> blk = qbd_mapmap1_blocks(arrival, service);
    const QbdFundMat<T> fm = qbd_fundmat(blk.B, blk.L, blk.F);
    QbdRg<T> out;
    out.R = fm.R;
    out.G = fm.G;
    out.B = blk.B;
    out.L = blk.L;
    out.F = blk.F;
    out.Lbar = blk.Lbar;
    out.U = qbd_detail::madd(blk.L, matmul(fm.R, blk.B));
    return out;
}

/** qbd_rg without rescaling the service process. */
template <class T>
QbdRg<T> qbd_rg(const Map<T>& arrival, const Map<T>& service) {
    return qbd_rg(arrival, service, T(num_traits<T>::from_int(0)));
}

/**
 * Factorial moment of order m of the number in system, computed in closed form
 * from the boundary vector and R:
 *
 *     E[N(N-1)...(N-m+1)] = m! pi_0 R^m (I - R)^-(m+1) e.
 *
 * Un-gated: given pi_0 and R this is exact matrix algebra. m = 0 returns 1.
 */
template <class T>
T qbd_qlen_factmoment(const std::vector<T>& pi0, const Matrix<T>& R, unsigned m) {
    const std::size_t n = R.rows();
    if (pi0.size() != n) throw InputError("qbd_qlen_factmoment: pi0 length mismatch");
    if (m == 0) return num_traits<T>::from_int(1);
    const Matrix<T> ImRinv = inverse(qbd_detail::msub(eye<T>(n), R));
    const Matrix<T> A = matmul(matpow(R, m), matpow(ImRinv, m + 1));
    const std::vector<T> v = vecmul(pi0, A);
    T s = num_traits<T>::from_int(0);
    for (const T& x : v) s += x;
    return num_factorial<T>(m) * s;
}

/**
 * Raw moment of order m of the number in system, E[N^m], assembled from the
 * factorial moments with the Stirling numbers of the second kind,
 * N^m = sum_j S(m,j) N(N-1)...(N-j+1). Integer combinatorics plus exact matrix
 * algebra, so un-gated and exact at Rational.
 */
template <class T>
T qbd_qlen_moment(const std::vector<T>& pi0, const Matrix<T>& R, unsigned m) {
    if (m == 0) return num_traits<T>::from_int(1);
    // S(i,j) by the recurrence S(i,j) = j S(i-1,j) + S(i-1,j-1).
    std::vector<std::vector<long long>> S(m + 1, std::vector<long long>(m + 1, 0));
    S[0][0] = 1;
    for (unsigned i = 1; i <= m; ++i)
        for (unsigned j = 1; j <= i; ++j)
            S[i][j] = static_cast<long long>(j) * S[i - 1][j] + S[i - 1][j - 1];
    T acc = num_traits<T>::from_int(0);
    for (unsigned j = 1; j <= m; ++j)
        acc += num_traits<T>::from_int(static_cast<long>(S[m][j])) *
               qbd_qlen_factmoment(pi0, R, j);
    return acc;
}

/** Result of qbd_mapmap1, mirroring the MATLAB return list. */
template <class T>
struct QbdMapMap1Result {
    T XN;   ///< throughput, = lambda of the arrival MAP
    T QN;   ///< mean number in system, closed form
    T UN;   ///< utilization, = 1 - sum(pi_0)
    T RN;   ///< mean response time, QN / XN by Little's law
    T eta;  ///< caudal characteristic, sp(R)
    Matrix<T> pqueue;    ///< level distribution, row k = pi_k
    std::vector<T> pi0;  ///< boundary vector, normalized
    Matrix<T> R, G, U;
    Matrix<T> B, L, F, Lbar;
    Map<T> service;  ///< the service MAP actually used (rescaled if util was given)
    unsigned iterations = 0;
};

/**
 * MAP/MAP/1 queue (qbd_mapmap1.m).
 *
 * Differences from the MATLAB reference, both deliberate:
 *   - QN is the closed form pi_0 R (I - R)^-2 e rather than the truncated sum
 *     over the levels that MATLAB and the JAR compute. The truncated value is
 *     available as qbd_mapmap1_qlen_truncated(res) for a like-for-like
 *     comparison.
 *   - eta is bracketed by Collatz-Wielandt on R (qbd_caudal) rather than taken
 *     from an eigendecomposition; the value is the same spectral radius.
 *
 * @param util if positive, rescale the service MAP to mean util / lambda_a
 * @param max_levels how many levels of pqueue to materialize
 * @param arrival the arrival MAP
 * @param service_in the service MAP, rescaled to the requested utilization
 */
template <class T>
QbdMapMap1Result<T> qbd_mapmap1(const Map<T>& arrival, const Map<T>& service_in, const T& util,
                                std::size_t max_levels) {
    static_assert(num_traits<T>::has_transcendental,
                  "qbd_mapmap1 requires transcendental arithmetic");
    Map<T> service = service_in;
    const T lambda_a = map_lambda(arrival);
    if (util > num_traits<T>::from_int(0)) service = map_scale(service, T(util / lambda_a));
    const T lambda_s = map_lambda(service);
    if (lambda_a >= lambda_s)
        throw NumericError("qbd_mapmap1: the queue is not stable, lambda_a >= lambda_s");

    const QbdMapMap1Blocks<T> blk = qbd_mapmap1_blocks(arrival, service);
    const QbdFundMat<T> fm = qbd_fundmat(blk.B, blk.L, blk.F);

    QbdMapMap1Result<T> res;
    res.R = fm.R;
    res.G = fm.G;
    res.iterations = fm.iterations;
    res.B = blk.B;
    res.L = blk.L;
    res.F = blk.F;
    res.Lbar = blk.Lbar;
    res.U = qbd_detail::madd(blk.L, matmul(fm.R, blk.B));
    res.service = service;

    res.pqueue = qbd_pi(blk.B, blk.Lbar, fm.R, max_levels, T(num_traits<T>::from_double(1e-10)));
    res.pi0.assign(res.pqueue.cols(), num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < res.pqueue.cols(); ++j) res.pi0[j] = res.pqueue(0, j);

    T s0 = num_traits<T>::from_int(0);
    for (const T& v : res.pi0) s0 += v;
    res.UN = num_traits<T>::from_int(1) - s0;
    res.QN = qbd_qlen_factmoment(res.pi0, fm.R, 1);
    res.XN = lambda_a;
    res.RN = res.QN / res.XN;
    res.eta = qbd_caudal(fm.R);
    return res;
}

/** qbd_mapmap1 without rescaling and with 20000 materialized levels. */
template <class T>
QbdMapMap1Result<T> qbd_mapmap1(const Map<T>& arrival, const Map<T>& service) {
    return qbd_mapmap1(arrival, service, T(num_traits<T>::from_int(0)),
                       static_cast<std::size_t>(20000));
}

/**
 * The mean number in system computed the way MATLAB's qbd_mapmap1 does it, by
 * summing k over the materialized levels. Provided for a like-for-like
 * comparison against the reference; qbd_mapmap1's QN field is the closed form
 * and is the value to use.
 */
template <class T>
T qbd_mapmap1_qlen_truncated(const QbdMapMap1Result<T>& res) {
    T acc = num_traits<T>::from_int(0);
    for (std::size_t k = 1; k < res.pqueue.rows(); ++k) {
        T lvl = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < res.pqueue.cols(); ++j) lvl += res.pqueue(k, j);
        acc += num_traits<T>::from_int(static_cast<long>(k)) * lvl;
    }
    return acc;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_QBD_MAPMAP1_H
