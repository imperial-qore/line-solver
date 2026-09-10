/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CHOW_H
#define LINE_API_PFQN_CHOW_H

/**
 * Chow Second Approximation (SA) approximate MVA.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_chow.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/Pfqn_chow.java. W.-M. Chow,
 * "Approximations for large scale closed queueing networks", Perform. Eval.
 * 3(1), 1983. The arrival-instant queue length is written exactly as
 *
 *   A_k^(c)(N) = Q_k(N - 1_c) = Q_k(N) (1 + theta_ck),
 *   theta_ck   = [Q_k(N - 1_c) - Q_k(N)] / Q_k(N),
 *
 * and the theta-terms are estimated ONCE, off the Bard LCP solution, before
 * the fixed point is run. Two estimators are given: the BACKWARD one uses
 * Qhat(N - 1_c), the FORWARD one Qhat(N + 1_c). Chow reports the forward form
 * to be the more accurate of the two, so it is the default here. Setting every
 * theta to zero recovers pfqn_lcp.
 *
 * Arithmetic: as pfqn_lcp, field operations only, with the same stopping-rule
 * caveat.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_lcp.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Which finite difference of the LCP solution estimates the theta-terms. */
enum class ChowVariant { Forward, Backward };

/**
 * @param L    (M x R) demands
 * @param N    (R) populations
 * @param Z    (R) think times, empty for none
 * @param type (M) per-station scheduling, empty for all PS
 * @param tol convergence tolerance
 * @param maxiter iteration cap
 * @param QN0 warm start for the inner LCP solves and the fixed point; may be empty
 * @param variant estimator of the theta-terms
 */
template <class T>
AmvaResult<T> pfqn_chow(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                        const std::vector<AmvaSched>& type, double tol = 1e-6,
                        std::size_t maxiter = 1000, const Matrix<T>& QN0 = Matrix<T>(),
                        ChowVariant variant = ChowVariant::Forward) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_chow: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_chow: Z has the wrong length");
    if (!type.empty() && type.size() != M) throw InputError("pfqn_chow: type has the wrong length");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    // theta-terms from the LCP solution
    const AmvaResult<T> base = pfqn_lcp(L, N, Z, type, tol, maxiter, QN0);
    std::vector<T> Qtot(M, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t s = 0; s < R; ++s) Qtot[i] += base.QN(i, s);

    Matrix<T> theta(M, R, zero);
    for (std::size_t cls = 0; cls < R; ++cls) {
        if (N[cls] == zero) continue;
        std::vector<T> Nalt(N);
        Nalt[cls] = (variant == ChowVariant::Backward) ? T(N[cls] - one) : T(N[cls] + one);
        const AmvaResult<T> alt = pfqn_lcp(L, Nalt, Z, type, tol, maxiter, QN0);
        for (std::size_t i = 0; i < M; ++i) {
            T qalt = zero;
            for (std::size_t s = 0; s < R; ++s) qalt += alt.QN(i, s);
            const T ref = (variant == ChowVariant::Backward) ? Qtot[i] : qalt;
            const T delta = (variant == ChowVariant::Backward) ? T(qalt - Qtot[i])
                                                               : T(Qtot[i] - qalt);
            if (ref > zero) theta(i, cls) = delta / ref;
        }
    }

    // fixed point with A_k^(c) = Q_k (1 + theta_ck)
    AmvaResult<T> r;
    r.XN.assign(R, zero);
    r.QN = Matrix<T>(M, R, zero);
    r.UN = Matrix<T>(M, R, zero);
    r.RN = Matrix<T>(M, R, zero);
    Matrix<T> CN(M, R, zero);
    if (M == 0) return r;
    if (!QN0.empty()) {
        if (QN0.rows() != M || QN0.cols() != R)
            throw InputError("pfqn_chow: QN0 has the wrong shape");
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
                const T infl = T(one + theta(i, cls));
                for (std::size_t s = 0; s < R; ++s) {
                    if (fcfs && s != cls)
                        CN(i, cls) += L(i, s) * r.QN(i, s) * infl;
                    else
                        CN(i, cls) += L(i, cls) * r.QN(i, s) * infl;
                }
                // a theta below -1 would make the arrival-instant queue negative
                if (CN(i, cls) < L(i, cls)) CN(i, cls) = L(i, cls);
                ctot += CN(i, cls);
            }
            if (ctot == zero) throw NumericError("pfqn_chow: zero total residence time");
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
                    if (r.QN(i, cls) == zero) continue;
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
AmvaResult<T> pfqn_chow(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    return pfqn_chow(L, N, Z, std::vector<AmvaSched>());
}

template <class T>
AmvaResult<T> pfqn_chow(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_chow(L, N, std::vector<T>(), std::vector<AmvaSched>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CHOW_H
