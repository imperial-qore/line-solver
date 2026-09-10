/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_BS_H
#define LINE_API_PFQN_BS_H

/**
 * Bard-Schweitzer approximate MVA.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_bs.m. The exact arrival theorem
 * Q(i|n - e_r) is replaced by the proportional estimate
 *   Q(i,r|n - e_r) = Q(i,r|n) (N_r - 1)/N_r,   Q(i,s|n - e_r) = Q(i,s|n),
 * and the resulting fixed point is iterated to a relative tolerance on Q.
 *
 * FCFS stations use the other class's own demand in the queueing term
 * (L(i,s) Q(i,s)), the PS family uses the arriving class's demand
 * (L(i,r) Q(i,s)); that distinction is what makes the FCFS variant sensitive
 * to demand heterogeneity, and it is easy to lose when transcribing.
 *
 * The iteration stops on a tolerance, so the result is a fixed point only to
 * within tol whatever the arithmetic. That is a caveat on what the answer
 * MEANS, not a reason to deny the exact backend: an exact run returns the
 * iterate the stopping rule selected, without rounding error, which is what
 * one wants when separating arithmetic error from algorithmic error.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_cntol.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Station scheduling as far as the AMVA formulas distinguish it. */
/// INF marks a delay centre; pfqn_bs only distinguishes FCFS, pfqn_qsa needs it.
enum class AmvaSched { PS, FCFS, INF };

template <class T>
struct AmvaResult {
    std::vector<T> XN;  ///< (R) throughput
    Matrix<T> QN;       ///< (M x R) queue length
    Matrix<T> UN;       ///< (M x R) utilization
    Matrix<T> RN;       ///< (M x R) residence time
    std::size_t iterations = 0;
    bool converged = false;
};

/**
 * @param L    (M x R) demands
 * @param N    (R) populations
 * @param Z    (R) think times, empty for none
 * @param type (M) per-station scheduling, empty for all PS
 * @param tol convergence tolerance; NaN selects the published Linearizer termination
 *            test of Chandy and Neuse, Commun. ACM 25(2), 1982, i.e. the cutoff
 *            pfqn_cntol(N) applied to max_{i,r}|dQ(i,r)|/N_r instead of the
 *            relative-change metric used by default. This is the test LQNS runs,
 *            since it sets it in SchweitzerCommon.
 * @param maxiter iteration cap
 * @param QN0 queue lengths that warm-start the iteration; empty for a cold start
 */
template <class T>
AmvaResult<T> pfqn_bs(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                      const std::vector<AmvaSched>& type, double tol = 1e-6,
                      std::size_t maxiter = 1000, const Matrix<T>& QN0 = Matrix<T>()) {
    // field-arithmetic exactness caveat: see _kb/03-api-layer.md (cpp port notes: pfqn)
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_bs: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_bs: Z has the wrong length");
    if (!type.empty() && type.size() != M) throw InputError("pfqn_bs: type has the wrong length");

    const bool cntest = is_cntol(tol);
    if (cntest) tol = pfqn_cntol(N);

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    AmvaResult<T> r;
    r.XN.assign(R, zero);
    r.QN = Matrix<T>(M, R, zero);
    r.UN = Matrix<T>(M, R, zero);
    r.RN = Matrix<T>(M, R, zero);
    Matrix<T> CN(M, R, zero);

    // QN0 warm-start rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    if (M == 0) return r;
    if (!QN0.empty()) {
        if (QN0.rows() != M || QN0.cols() != R)
            throw InputError("pfqn_bs: QN0 has the wrong shape");
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
                if (L(i, cls) == zero) {
                    continue;
                }
                const bool fcfs = !type.empty() && type[i] == AmvaSched::FCFS;
                for (std::size_t s = 0; s < R; ++s) {
                    if (s != cls)
                        CN(i, cls) += fcfs ? L(i, s) * r.QN(i, s) : L(i, cls) * r.QN(i, s);
                    else
                        CN(i, cls) += L(i, cls) * r.QN(i, cls) * (N[cls] - one) / N[cls];
                }
                ctot += CN(i, cls);
            }
            if (ctot == zero) throw NumericError("pfqn_bs: zero total residence time");
            r.XN[cls] = N[cls] / ctot;
        }
        for (std::size_t cls = 0; cls < R; ++cls)
            for (std::size_t i = 0; i < M; ++i) {
                r.QN(i, cls) = r.XN[cls] * CN(i, cls);
                r.UN(i, cls) = r.XN[cls] * L(i, cls);
            }

        // 0/0 vs x/0 convergence rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
        double delta = 0.0;
        for (std::size_t cls = 0; cls < R; ++cls) {
            if (N[cls] == zero) continue;
            for (std::size_t i = 0; i < M; ++i) {
                if (cntest) {
                    // Chandy and Neuse (1982), p.129: absolute queue-length change scaled by
                    // the class population, not the relative change.
                    const double d = std::fabs(num_traits<T>::to_double(
                        T((r.QN(i, cls) - Qprev(i, cls)) / N[cls])));
                    if (d > delta) delta = d;
                    continue;
                }
                if (Qprev(i, cls) == zero) {
                    if (r.QN(i, cls) == zero) continue;  // 0/0, omitted by max
                    delta = std::numeric_limits<double>::infinity();
                    continue;
                }
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
AmvaResult<T> pfqn_bs(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    return pfqn_bs(L, N, Z, std::vector<AmvaSched>());
}

template <class T>
AmvaResult<T> pfqn_bs(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_bs(L, N, std::vector<T>(), std::vector<AmvaSched>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_BS_H
