/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PAM_H
#define LINE_API_PFQN_PAM_H

/**
 * Hsieh-Lam Proportional Approximation Methods (PAMB / PAMI / PAMT).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_pam.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/Pfqn_pam.java. C. T. Hsieh, S. S. Lam,
 * "PAM - A noniterative approximate solution method for closed multichain
 * queueing networks", ACM SIGMETRICS Perform. Eval. Rev. 16(1), 1988. The
 * three variants are NONITERATIVE: the queue lengths are seeded by the
 * proportion of a class demand that falls at each centre,
 *
 *   E_ck = D_ck / sum_i D_ci,   Q_ck(N) = E_ck N_c,
 *
 * and the MVA equations are then unrolled a fixed number of times. PAMB
 * applies the last MVA step; PAMI additionally scales a class down wherever it
 * would drive a centre past full utilization; PAMT seeds at N - 1_i - 1_j and
 * applies the last TWO MVA steps before that capping.
 *
 * The seed spreads the whole class population over the queueing centres and
 * ignores Z, exactly as published: PAM buys speed, not accuracy.
 *
 * Arithmetic: sums, products, divisions and comparisons only, and there is no
 * fixed point, so the result is EXACT in rational arithmetic.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_bs.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Which of the three proportional approximations to run. */
enum class PamVariant { Basic, Improved, Two };

/**
 * @param L (M x R) demands, @param N (R) populations,
 * @param Z (R) think times (empty for none), @param variant PAMB/PAMI/PAMT
 */
template <class T>
AmvaResult<T> pfqn_pam(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                       PamVariant variant = PamVariant::Basic) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_pam: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_pam: Z has the wrong length");

    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    AmvaResult<T> r;
    r.XN.assign(R, zero);
    r.QN = Matrix<T>(M, R, zero);
    r.UN = Matrix<T>(M, R, zero);
    r.RN = Matrix<T>(M, R, zero);
    r.iterations = 1;
    r.converged = true;
    if (M == 0) return r;

    // E_ck, the share of the class-c demand served at centre k
    Matrix<T> E(M, R, zero);
    for (std::size_t s = 0; s < R; ++s) {
        T tot = zero;
        for (std::size_t i = 0; i < M; ++i) tot += L(i, s);
        if (tot > zero)
            for (std::size_t i = 0; i < M; ++i) E(i, s) = L(i, s) / tot;
    }
    Matrix<T> Q(M, R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t s = 0; s < R; ++s) Q(i, s) = E(i, s) * N[s];

    if (variant == PamVariant::Two) {
        for (std::size_t i = 0; i < R; ++i) {
            Matrix<T> Qmi(M, R, zero);  // Q_jk(N - 1_i)
            for (std::size_t j = 0; j < R; ++j) {
                // Q_ck(N - 1_i - 1_j) = Q_ck(N) - E_ck [(c==i) + (c==j)]
                std::vector<T> Rj(M, zero);
                T rtot = zero;
                for (std::size_t k = 0; k < M; ++k) {
                    T agg = zero;
                    for (std::size_t c = 0; c < R; ++c) {
                        T q = Q(k, c);
                        if (c == i) q -= E(k, c);
                        if (c == j) q -= E(k, c);
                        agg += q;
                    }
                    Rj[k] = L(k, j) * T(one + agg);
                    rtot += Rj[k];
                }
                const T nj = (i == j) ? T(N[j] - one) : N[j];
                const T zj = Z.empty() ? zero : Z[j];
                const T den = T(rtot + zj);
                const T Xj = (nj > zero && den > zero) ? T(nj / den) : zero;
                for (std::size_t k = 0; k < M; ++k) Qmi(k, j) = Xj * Rj[k];
            }
            T rtot = zero;
            for (std::size_t k = 0; k < M; ++k) {
                T agg = zero;
                for (std::size_t c = 0; c < R; ++c) agg += Qmi(k, c);
                r.RN(k, i) = L(k, i) * T(one + agg);
                rtot += r.RN(k, i);
            }
            const T zi = Z.empty() ? zero : Z[i];
            if (N[i] > zero && T(rtot + zi) > zero) r.XN[i] = N[i] / T(rtot + zi);
        }
    } else {
        for (std::size_t s = 0; s < R; ++s) {
            // Q_jk(N - 1_s) = Q_jk(N) - E_jk [j == s]
            T rtot = zero;
            for (std::size_t k = 0; k < M; ++k) {
                T agg = zero;
                for (std::size_t c = 0; c < R; ++c) agg += Q(k, c) - (c == s ? E(k, c) : zero);
                r.RN(k, s) = L(k, s) * T(one + agg);
                rtot += r.RN(k, s);
            }
            const T zs = Z.empty() ? zero : Z[s];
            if (N[s] > zero && T(rtot + zs) > zero) r.XN[s] = N[s] / T(rtot + zs);
        }
    }

    if (variant != PamVariant::Basic) {
        // scale a class down when it would drive a centre it visits past U = 1
        std::vector<T> U(M, zero);
        for (std::size_t k = 0; k < M; ++k)
            for (std::size_t c = 0; c < R; ++c) U[k] += L(k, c) * r.XN[c];
        for (std::size_t s = 0; s < R; ++s) {
            bool any = false;
            T best = zero;
            for (std::size_t k = 0; k < M; ++k) {
                if (L(k, s) == zero) continue;
                if (!any || U[k] > best) {
                    best = U[k];
                    any = true;
                }
            }
            if (any && best > one) r.XN[s] = r.XN[s] / best;
        }
    }

    for (std::size_t k = 0; k < M; ++k)
        for (std::size_t s = 0; s < R; ++s) {
            r.QN(k, s) = r.XN[s] * r.RN(k, s);
            r.UN(k, s) = r.XN[s] * L(k, s);
        }
    return r;
}

template <class T>
AmvaResult<T> pfqn_pam(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_pam(L, N, std::vector<T>(), PamVariant::Basic);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PAM_H
