/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_TAY_H
#define LINE_API_PFQN_TAY_H

/**
 * Tay's arrival-instant approximate MVA.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_tay.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/Pfqn_tay.java.
 *
 * WHAT MAKES IT DIFFERENT. Every other AMVA in this directory estimates the
 * arrival-instant queue length by shifting the population: Bard-Schweitzer
 * scales Q by (N-1)/N, Linearizer solves R auxiliary networks, AQL carries an
 * aggregate correction. Tay estimates it from the THROUGHPUT ELASTICITIES
 * instead. With E_mkc = (D_mk/X_c) dX_c/dD_mk the elasticity of the class-c
 * throughput in the class-k demand at station m, and B_ir = 1/(1 + D_ir X_r/N_r),
 * the elasticities satisfy R linear equations
 *
 *   E_mkj sum_t B_tj Q_jt (1+Q_jt)
 *     = -[(delta_jk + Q_jm) B_mk Q_km + sum_{c/=j} E_mkc sum_t B_tc Q_jt Q_ct]
 *
 * and the arrival-instant queue length is then simply Q_km^(r) = Q_km + E_mkr,
 * which closes the recursion R_rm = D_rm (1 + sum_k Q_km^(r)).
 *
 * COST. One R x R solve per (station, class) pair per sweep, so O(M R (R^3 + M R^2))
 * per iteration: more than Bard-Schweitzer, less than Linearizer's R+1 auxiliary
 * networks.
 *
 * DELAY STATIONS enter through Z only. They are "AS" servers in the survey's
 * notation (d_t = 0), contributing Z_j X_j to the DENOMINATOR of the elasticity
 * equations and nothing to the numerator.
 *
 * EMPTY CLASSES are solved out and re-expanded, exactly as pfqn_bs does. Their
 * elasticity denominator is identically zero, so leaving them in makes the R x R
 * system singular rather than merely redundant.
 *
 * Reference: Y. C. Tay and R. Suri, "Error bounds for performance prediction in
 * queueing networks", ACM TOCS 3(4), 1985; Y. C. Tay, "An approach to analyzing
 * the behavior of some queueing networks", Operations Research 40(S2), 1992;
 * P. J. Schweitzer, G. Serazzi and M. Broglia, "A survey of bottleneck analysis
 * in closed queueing networks", Sec. 4.8.2, eqs. 4.8.2-1..3.
 *
 * Iterates to an absolute tolerance on the queue lengths, so exact arithmetic
 * buys nothing and the static_assert records that, as in pfqn_aql.h.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_bs.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param L       (M x R) demands
 * @param N       (R) populations
 * @param Z       (R) think times, empty for none
 * @param tol     absolute tolerance on the queue lengths
 * @param maxiter iteration cap
 * @param QN0     (M x R) initial queue lengths, empty for uniform
 *
 * `AmvaResult::RN` holds the residence times. The arrival-instant queue lengths
 * QNarr(m,k,r) that the method is tabulated on are NOT returned: the reference
 * exposes them as a sixth output for diagnostics only, and no caller in this
 * port reads them. They are the auxiliary quantities of the approximation, not
 * the model re-solved at N - e_r, which is the same object only for an exact
 * solution.
 */
template <class T>
AmvaResult<T> pfqn_tay(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                       double tol = 1e-6, std::size_t maxiter = 1000,
                       const Matrix<T>& QN0 = Matrix<T>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_tay requires transcendental arithmetic: it iterates to a tolerance, so "
                  "its answer is a fixed point only to within tol");
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_tay: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_tay: Z has the wrong length");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    AmvaResult<T> out;
    out.XN.assign(R, zero);
    out.QN = Matrix<T>(M, R, zero);
    out.UN = Matrix<T>(M, R, zero);
    out.RN = Matrix<T>(M, R, zero);
    if (M == 0) return out;

    // Empty classes make the elasticity system singular, not merely redundant:
    // their denominator is identically zero. Solve without them and re-expand.
    std::vector<std::size_t> act;
    for (std::size_t r = 0; r < R; ++r)
        if (N[r] > zero) act.push_back(r);
    if (act.empty()) return out;
    if (act.size() < R) {
        Matrix<T> La(M, act.size(), zero);
        std::vector<T> Na(act.size(), zero), Za(act.size(), zero);
        for (std::size_t a = 0; a < act.size(); ++a) {
            for (std::size_t m = 0; m < M; ++m) La(m, a) = L(m, act[a]);
            Na[a] = N[act[a]];
            Za[a] = Z.empty() ? zero : Z[act[a]];
        }
        const AmvaResult<T> sub = pfqn_tay(La, Na, Za, tol, maxiter);
        for (std::size_t a = 0; a < act.size(); ++a) {
            out.XN[act[a]] = sub.XN[a];
            for (std::size_t m = 0; m < M; ++m) {
                out.QN(m, act[a]) = sub.QN(m, a);
                out.UN(m, act[a]) = sub.UN(m, a);
                out.RN(m, act[a]) = sub.RN(m, a);
            }
        }
        out.iterations = sub.iterations;
        out.converged = sub.converged;
        return out;
    }

    Matrix<T> QN(M, R, zero);
    if (QN0.empty()) {
        const T Md = num_traits<T>::from_int(static_cast<long>(M));
        for (std::size_t m = 0; m < M; ++m)
            for (std::size_t r = 0; r < R; ++r) QN(m, r) = N[r] / Md;
    } else {
        if (QN0.rows() != M || QN0.cols() != R)
            throw InputError("pfqn_tay: QN0 has the wrong shape");
        QN = QN0;
    }
    // XN = N ./ (Z + sum(L,1).*(1+sum(QN,1)))
    for (std::size_t r = 0; r < R; ++r) {
        T Lsum = zero, Qsum = zero;
        for (std::size_t m = 0; m < M; ++m) {
            Lsum += L(m, r);
            Qsum += QN(m, r);
        }
        const T den = (Z.empty() ? zero : Z[r]) + Lsum * (one + Qsum);
        out.XN[r] = (den == zero) ? zero : N[r] / den;
    }

    Matrix<T> Qarr(M, R * R, zero);  // Qarr(m, k*R + r) = Q_km seen by class r
    for (std::size_t it = 1; it <= maxiter; ++it) {
        out.iterations = it;
        const Matrix<T> QN_1 = QN;

        // B(i,r) = 1/(1 + L(i,r) X_r/N_r)
        Matrix<T> B(M, R, zero);
        for (std::size_t m = 0; m < M; ++m)
            for (std::size_t r = 0; r < R; ++r)
                B(m, r) = one / (one + L(m, r) * out.XN[r] / N[r]);

        // den(j): the AS-server term Z_j X_j is the delay contribution, which
        // leaves B = 1 because a delay has d_t = 0.
        std::vector<T> den(R, zero);
        for (std::size_t j = 0; j < R; ++j) {
            T s = zero;
            for (std::size_t m = 0; m < M; ++m) s += B(m, j) * QN(m, j) * (one + QN(m, j));
            den[j] = s + (Z.empty() ? zero : Z[j]) * out.XN[j];
        }

        // C(j,c) = sum_t B_tc Q_jt Q_ct
        Matrix<T> C(R, R, zero);
        for (std::size_t j = 0; j < R; ++j)
            for (std::size_t c = 0; c < R; ++c) {
                T s = zero;
                for (std::size_t m = 0; m < M; ++m) s += B(m, c) * QN(m, j) * QN(m, c);
                C(j, c) = s;
            }

        for (std::size_t m = 0; m < M; ++m)
            for (std::size_t k = 0; k < R; ++k) {
                Matrix<T> A(R, R, zero);
                std::vector<T> b(R, zero);
                for (std::size_t j = 0; j < R; ++j) {
                    A(j, j) = one;
                    if (den[j] == zero)
                        throw NumericError(
                            "pfqn_tay: the elasticity system is singular at station " +
                            std::to_string(m + 1) + ", class " + std::to_string(j + 1) +
                            ": the class carries no queue anywhere and no think time");
                    for (std::size_t c = 0; c < R; ++c)
                        if (c != j) A(j, c) = C(j, c) / den[j];
                    const T delta = (j == k) ? one : zero;
                    b[j] = -((delta + QN(m, j)) * B(m, k) * QN(m, k) / den[j]);
                }
                const std::vector<T> E = line::solve(A, b);
                for (std::size_t r = 0; r < R; ++r) Qarr(m, k * R + r) = QN(m, k) + E[r];
            }

        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t m = 0; m < M; ++m) {
                T s = zero;
                for (std::size_t k = 0; k < R; ++k) s += Qarr(m, k * R + r);
                out.RN(m, r) = L(m, r) * (one + s);
            }
        for (std::size_t r = 0; r < R; ++r) {
            T s = (Z.empty() ? zero : Z[r]);
            for (std::size_t m = 0; m < M; ++m) s += out.RN(m, r);
            out.XN[r] = (s == zero) ? zero : N[r] / s;
        }
        for (std::size_t m = 0; m < M; ++m)
            for (std::size_t r = 0; r < R; ++r) QN(m, r) = out.RN(m, r) * out.XN[r];

        double delta = 0.0;
        for (std::size_t m = 0; m < M; ++m)
            for (std::size_t r = 0; r < R; ++r) {
                const double d =
                    std::fabs(num_traits<T>::to_double(T(QN(m, r) - QN_1(m, r))));
                if (d > delta) delta = d;
            }
        if (delta < tol) {
            out.converged = true;
            break;
        }
    }

    out.QN = QN;
    for (std::size_t m = 0; m < M; ++m)
        for (std::size_t r = 0; r < R; ++r) out.UN(m, r) = L(m, r) * out.XN[r];
    return out;
}

template <class T>
AmvaResult<T> pfqn_tay(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_tay(L, N, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_TAY_H
