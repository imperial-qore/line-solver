/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LCP_H
#define LINE_API_PFQN_LCP_H

/**
 * Bard Large Customer Population (LCP) approximate MVA.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_lcp.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/Pfqn_lcp.java. Y. Bard, "Some
 * extensions to multiclass queueing network analysis", in Performance of
 * Computer Systems, North-Holland, 1979: the first approximate MVA algorithm.
 * It estimates the arrival-instant queue length by the time-averaged one
 * WITHOUT removing the arriving customer,
 *
 *   A_k^(c)(N) = Q_k(N - 1_c) ~= Q_k(N) = sum_s Q_ks(N),
 *
 * since with a large population one customer less cannot change the mean queue
 * lengths appreciably. Dropping the Bard-Schweitzer proportional factor
 * (N_r - 1)/N_r from pfqn_bs is exactly this algorithm, so LCP is uniformly
 * more pessimistic than pfqn_bs and is inaccurate at small populations.
 *
 * Arithmetic: sums, products and divisions only, so the iterate is EXACT in
 * rational arithmetic. The stopping rule still selects WHICH iterate is
 * returned, the same caveat pfqn_bs carries.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_bs.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param L    (M x R) demands
 * @param N    (R) populations
 * @param Z    (R) think times, empty for none
 * @param type (M) per-station scheduling, empty for all PS
 * @param tol convergence tolerance
 * @param maxiter iteration cap
 * @param QN0 queue lengths that warm-start the iteration; empty for a cold start
 */
template <class T>
AmvaResult<T> pfqn_lcp(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                       const std::vector<AmvaSched>& type, double tol = 1e-6,
                       std::size_t maxiter = 1000, const Matrix<T>& QN0 = Matrix<T>()) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_lcp: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_lcp: Z has the wrong length");
    if (!type.empty() && type.size() != M) throw InputError("pfqn_lcp: type has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    AmvaResult<T> r;
    r.XN.assign(R, zero);
    r.QN = Matrix<T>(M, R, zero);
    r.UN = Matrix<T>(M, R, zero);
    r.RN = Matrix<T>(M, R, zero);
    Matrix<T> CN(M, R, zero);

    if (M == 0) return r;
    if (!QN0.empty()) {
        if (QN0.rows() != M || QN0.cols() != R)
            throw InputError("pfqn_lcp: QN0 has the wrong shape");
        r.QN = QN0;
    } else {
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t s = 0; s < R; ++s)
                r.QN(i, s) = N[s] / num_traits<T>::from_int(static_cast<long>(M));
    }

    for (std::size_t it = 1; it <= maxiter; ++it) {
        r.iterations = it;
        const Matrix<T> Qprev = r.QN;

        for (std::size_t cls = 0; cls < R; ++cls) {
            if (N[cls] == zero) {
                r.XN[cls] = zero;
                for (std::size_t i = 0; i < M; ++i) {
                    CN(i, cls) = zero;
                    r.QN(i, cls) = zero;
                    r.UN(i, cls) = zero;
                }
                continue;
            }
            T ctot = Z.empty() ? zero : Z[cls];
            for (std::size_t i = 0; i < M; ++i) {
                CN(i, cls) = L(i, cls);
                if (L(i, cls) == zero) continue;
                const bool fcfs = !type.empty() && type[i] == AmvaSched::FCFS;
                for (std::size_t s = 0; s < R; ++s) {
                    // the arriving customer is NOT removed: no (N-1)/N factor
                    if (fcfs && s != cls)
                        CN(i, cls) += L(i, s) * r.QN(i, s);
                    else
                        CN(i, cls) += L(i, cls) * r.QN(i, s);
                }
                ctot += CN(i, cls);
            }
            if (ctot == zero) throw NumericError("pfqn_lcp: zero total residence time");
            r.XN[cls] = N[cls] / ctot;
        }
        for (std::size_t cls = 0; cls < R; ++cls)
            for (std::size_t i = 0; i < M; ++i) {
                r.QN(i, cls) = r.XN[cls] * CN(i, cls);
                r.UN(i, cls) = r.XN[cls] * L(i, cls);
            }

        double delta = 0.0;
        for (std::size_t cls = 0; cls < R; ++cls) {
            if (N[cls] == zero) continue;
            for (std::size_t i = 0; i < M; ++i) {
                if (Qprev(i, cls) == zero) {
                    if (r.QN(i, cls) == zero) continue;  // 0/0, omitted by max
                    delta = std::numeric_limits<double>::infinity();
                    continue;
                }
                const T one = num_traits<T>::from_int(1);
                const double d =
                    std::fabs(num_traits<T>::to_double(T(one - r.QN(i, cls) / Qprev(i, cls))));
                if (d > delta) delta = d;
            }
        }
        if (delta < tol) {
            r.converged = true;
            break;
        }
    }

    for (std::size_t cls = 0; cls < R; ++cls)
        for (std::size_t i = 0; i < M; ++i)
            r.RN(i, cls) = (N[cls] == zero) ? zero : r.QN(i, cls) / r.XN[cls];
    return r;
}

template <class T>
AmvaResult<T> pfqn_lcp(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    return pfqn_lcp(L, N, Z, std::vector<AmvaSched>());
}

template <class T>
AmvaResult<T> pfqn_lcp(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_lcp(L, N, std::vector<T>(), std::vector<AmvaSched>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LCP_H
