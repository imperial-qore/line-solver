/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_DMLIN_H
#define LINE_API_PFQN_DMLIN_H

/**
 * de Souza e Silva-Muntz Improved Linearizer (IL).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_dmlin.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/Pfqn_dmlin.java. E. de Souza e Silva,
 * R. R. Muntz, "A note on the computational cost of the Linearizer algorithm
 * for queueing networks", IEEE Trans. Computers 39(6), 1990. Linearizer
 * evaluates the arrival-instant queue length as
 *
 *   A_k^(c)(n) = sum_i (n_i - delta_c^(i)) [Q_ik(n)/n_i + Delta^(i)_ck],
 *
 * re-summing the C Delta-terms at every Core iteration, at every one of the
 * C+1 populations: O(K C^3) per refresh pass. IL splits that sum into the part
 * that moves with the Core iterate and the part that does not,
 *
 *   A_k^(c)(n)     = sum_i (n_i - delta_c^(i)) Q_ik(n)/n_i + xi_ck(n),
 *   xi_ck(N)       = sum_i (N_i - delta_c^(i)) Delta^(i)_ck,
 *   xi_ck(N - 1_j) = xi_ck(N) - Delta^(j)_ck,
 *
 * so the C K aggregates xi are computed ONCE per refresh pass and each Core
 * iteration then costs O(K C) instead of O(K C^2). Time drops to O(K C^2) with
 * the space unchanged at O(K C^2), and, because the split is an IDENTITY and
 * not an approximation, the fixed point is the one Linearizer reaches:
 * pfqn_dmlin and pfqn_linearizer agree to round-off. Transcribing (2.50) of
 * the Wang (1997) survey literally -- xi_ck(N - 1_j) ~= xi_ck(N), dropping the
 * Delta^(j)_ck correction -- breaks that agreement and costs about an order of
 * magnitude of accuracy, so the correction is not optional.
 *
 * Arithmetic: field operations only, so each iterate is EXACT in rational
 * arithmetic; the Core stops on enorm(Q_{k+1} - Q_k) < tol, so the returned
 * value still depends on the stopping rule.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_egflinearizer.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/** Fixed point of the aggregated arrival-instant estimate with the MVA equations. */
template <class T>
struct DmlinCore {
    Matrix<T> Q;
    Matrix<T> W;
    std::vector<T> X;
    int iter;
};

template <class T>
DmlinCore<T> dmlin_core(const Matrix<T>& L, const std::vector<int>& N1, const std::vector<T>& Z,
                        const Matrix<T>& Qin, const Matrix<T>& xi, double tol, int maxiter) {
    const std::size_t M = L.rows(), R = L.cols();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    DmlinCore<T> out;
    out.Q = Qin;
    out.W = Matrix<T>(M, R, zero);
    out.X.assign(R, zero);
    out.iter = 0;
    bool converged = false;
    while (!converged) {
        const Matrix<T> Qlast = out.Q;
        for (std::size_t c = 0; c < R; ++c)
            for (std::size_t i = 0; i < M; ++i) {
                T acc = zero;
                for (std::size_t s = 0; s < R; ++s) {
                    if (N1[s] <= 0) continue;
                    const long nr = static_cast<long>(N1[s]) - (s == c ? 1L : 0L);
                    if (nr <= 0) continue;
                    acc += num_traits<T>::from_int(nr) * out.Q(i, s) /
                           num_traits<T>::from_int(static_cast<long>(N1[s]));
                }
                out.W(i, c) = L(i, c) * T(one + acc + xi(i, c));
            }
        for (std::size_t s = 0; s < R; ++s) {
            T wsum = zero;
            for (std::size_t i = 0; i < M; ++i) wsum += out.W(i, s);
            const T zs = Z.empty() ? zero : Z[s];
            if (N1[s] > 0) {
                const T den = T(zs + wsum);
                if (den == zero) throw NumericError("pfqn_dmlin: zero cycle time");
                out.X[s] = num_traits<T>::from_int(static_cast<long>(N1[s])) / den;
            } else {
                out.X[s] = zero;
            }
            for (std::size_t i = 0; i < M; ++i) out.Q(i, s) = out.X[s] * out.W(i, s);
        }
        // enorm_diff, not enorm: sqrt is not an operation of the rational field,
        // so a T-valued norm cannot be instantiated for the exact backend at all
        // (pfqn_amva_common.h:62-72). line-cli instantiates this whole chain for
        // cpp_rational via `-a normconst`, so enorm(diff) broke the C++ build.
        // The sibling linearizers all use this form; it drops the temporary too.
        if (enorm_diff(out.Q, Qlast) < tol || out.iter > maxiter) converged = true;
        ++out.iter;
    }
    return out;
}

}  // namespace detail

/**
 * @param L       (M x R) service demands
 * @param N       (R) population per class
 * @param Z       (K x R) think times, summed over rows; may be empty
 * @param type    (M) scheduling discipline; accepted and unused, as in
 *                pfqn_linearizer, which treats every station as single-server PS
 * @param tol     convergence tolerance
 * @param maxiter total inner-iteration budget
 * @param QN0     (M x R) warm start of the Bard-Schweitzer seed; may be empty
 * @param npasses number of xi refresh passes (3, the Chandy-Neuse rule)
 */
template <class T>
LinearizerResult<T> pfqn_dmlin(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                               const std::vector<SchedStrategy>& type, double tol, int maxiter,
                               const Matrix<T>& QN0, int npasses = 3) {
    (void)type;
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_dmlin: L and N disagree on the class count");
    const T zero = num_traits<T>::from_int(0);

    std::vector<T> Zv(R, zero);
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_dmlin: Z has the wrong width");
        for (std::size_t s = 0; s < R; ++s)
            for (std::size_t k = 0; k < Z.rows(); ++k) Zv[s] += Z(k, s);
    }

    LinearizerResult<T> res;
    res.Q = Matrix<T>(M, R, zero);
    res.U = Matrix<T>(M, R, zero);
    res.W = Matrix<T>(M, R, zero);
    res.C.assign(R, zero);
    res.X.assign(R, zero);
    res.totiter = 0;

    bool allZero = true;
    for (std::size_t i = 0; i < M && allZero; ++i)
        for (std::size_t s = 0; s < R; ++s)
            if (L(i, s) != zero) {
                allZero = false;
                break;
            }
    if (M == 0 || allZero) {
        for (std::size_t s = 0; s < R; ++s) {
            if (Zv[s] != zero) res.X[s] = num_traits<T>::from_int(N[s]) / Zv[s];
            for (std::size_t i = 0; i < M; ++i) res.U(i, s) = res.X[s] * L(i, s);
        }
        return res;
    }

    // Initialize, as Linearizer does, from Bard-Schweitzer at every population
    std::vector<Matrix<T> > Qs(R + 1);
    for (std::size_t s = 0; s <= R; ++s) {
        const std::vector<int> N1 = oner(N, s);
        std::vector<T> Nt(R, zero);
        for (std::size_t c = 0; c < R; ++c) Nt[c] = num_traits<T>::from_int(N1[c]);
        const AmvaResult<T> seed =
            pfqn_bs(L, Nt, Zv, std::vector<AmvaSched>(), tol, static_cast<std::size_t>(maxiter),
                    QN0);
        Qs[s] = seed.QN;
    }

    // Delta[i][r][c] is the Delta^(r)_c term of station i
    std::vector<std::vector<std::vector<T> > > Delta(
        M, std::vector<std::vector<T> >(R, std::vector<T>(R, zero)));
    Matrix<T> xi(M, R, zero);

    for (int pass = 0; pass < npasses; ++pass) {
        for (std::size_t s = 0; s <= R; ++s) {
            const std::vector<int> N1 = oner(N, s);
            // xi at population N - 1_s, exactly; s == 0 leaves xi at N
            Matrix<T> xis(M, R, zero);
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t c = 0; c < R; ++c)
                    xis(i, c) = (s == 0) ? xi(i, c) : T(xi(i, c) - Delta[i][s - 1][c]);
            const detail::DmlinCore<T> cr =
                detail::dmlin_core(L, N1, Zv, Qs[s], xis, tol, maxiter - res.totiter);
            Qs[s] = cr.Q;
            res.totiter += cr.iter;
        }
        // Refresh the Delta-terms, then aggregate them into xi once per pass
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < R; ++r) {
                if (N[r] == 1) Qs[r + 1](i, r) = zero;
                for (std::size_t s = 1; s <= R; ++s) {
                    const long ns = static_cast<long>(N[r]) - (r == s - 1 ? 1L : 0L);
                    if (N[r] > 0 && ns > 0)
                        Delta[i][r][s - 1] = Qs[s](i, r) / num_traits<T>::from_int(ns) -
                                             Qs[0](i, r) / num_traits<T>::from_int(N[r]);
                    else if (N[r] > 0)
                        Delta[i][r][s - 1] = zero - Qs[0](i, r) / num_traits<T>::from_int(N[r]);
                    else
                        Delta[i][r][s - 1] = zero;
                }
            }
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t c = 0; c < R; ++c) {
                T acc = zero;
                for (std::size_t r = 0; r < R; ++r) {
                    const long w = static_cast<long>(N[r]) - (r == c ? 1L : 0L);
                    if (w > 0) acc += num_traits<T>::from_int(w) * Delta[i][r][c];
                }
                xi(i, c) = acc;
            }
    }

    const detail::DmlinCore<T> fin =
        detail::dmlin_core(L, N, Zv, Qs[0], xi, tol, maxiter - res.totiter);
    res.totiter += fin.iter;
    res.Q = fin.Q;
    res.W = fin.W;
    res.X = fin.X;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t s = 0; s < R; ++s) res.U(i, s) = res.X[s] * L(i, s);
    for (std::size_t s = 0; s < R; ++s)
        res.C[s] = (res.X[s] == zero) ? zero
                                      : T(num_traits<T>::from_int(N[s]) / res.X[s] - Zv[s]);
    return res;
}

template <class T>
LinearizerResult<T> pfqn_dmlin(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z) {
    return pfqn_dmlin(L, N, Z, std::vector<SchedStrategy>(), 1e-8, 1000, Matrix<T>());
}

template <class T>
LinearizerResult<T> pfqn_dmlin(const Matrix<T>& L, const std::vector<int>& N) {
    return pfqn_dmlin(L, N, Matrix<T>(), std::vector<SchedStrategy>(), 1e-8, 1000, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_DMLIN_H
