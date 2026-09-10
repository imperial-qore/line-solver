/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MFQ_LD_SOLVE_H
#define LINE_API_MAM_MFQ_LD_SOLVE_H

/**
 * First- and second-order level-dependent (multi-regime) Markovian fluid
 * queues: the matrix-exponential building blocks of the stationary law.
 *
 * Port of matlab/src/api/mam/mfq_ld_solve.m and the BUTools
 * SecondOrderLevelDependentFluidSolve it wraps. The generator, the drift and
 * optionally the VARIANCE change at threshold levels, giving a piecewise
 * homogeneous first- or second-order (Brownian) fluid queue. Setting every
 * variance cell to zero reduces it to first order. The blocks it returns are
 * exactly what mfq_ld_mean and mfq_ld_distr consume, so this closes that
 * family rather than extending it.
 *
 * METHOD. Per regime, the states with neither drift nor variance are censored
 * out, and the remainder is split into three classes: positive drift with no
 * variance, negative drift with no variance, and any state with variance. The
 * density is anchored at both ends of the regime, forward from its lower
 * threshold and backward from its upper one, and each direction's exponent
 * comes from a QBD. The second-order terms are what force the QBD form: with a
 * variance the balance equation is second order in the level, and multiplying
 * through by a constant c chosen from the spectrum turns it into the matrix
 * QUADRATIC of a discrete QBD, from whose R the exponent is recovered as
 * K = (R - I) c. c is the smallest value that keeps the transformed triple
 * substochastic, taken as the maximum over the drift states of -Q_ii/R_ii and
 * over the variance states of the larger root of the discriminant, floored at 1.
 *
 * WHICH QBD SOLVER, AND WHY NO NEW ONE WAS NEEDED. The reference calls
 * QBD_CR(Bm, Lm, Fm), and the R it returns satisfies
 *
 *     Fm + R Lm + R^2 Bm = 0,
 *
 * which is EXACTLY the contract of the port's existing qbd_R,
 * F + R L + R^2 B = 0, under the direct mapping B = Bm, L = Lm, F = Fm. So the
 * cyclic-reduction routine did not have to be transcribed.
 *
 * That equation was established by MEASUREMENT, not by reading the argument
 * names, and the first reading was wrong. QBD_CR's own error message rejects a
 * triple whose sum is not "(sub)stochastic", which invites the conclusion that
 * it wants the discrete-time form R = A0 + R A1 + R^2 A2 with A0 = Bm; the
 * triples this caller builds sum to a GENERATOR, and the shift that reading
 * implies (L = Lm - I) produces a completely different and wrong K. Feeding the
 * actual regime-1 triple of a two-state first-order instance through MATLAB and
 * evaluating all four candidate residuals settles it in one run: 2.2e-16 for the
 * form above against 1.25, 1.5 and 2.75 for the others. The test asserts that
 * residual directly on the port's own R, so the identity is pinned rather than
 * inferred.
 *
 * The boundary system then couples the K+1 point masses to the 2K density
 * initial vectors through flux conservation at every threshold, the boundary
 * conditions that reflective and absorbing states impose, and continuity of the
 * second-order density across a threshold. One normalization row replaces the
 * first flux equation.
 *
 * REFERENCE DEFECT, reported and NOT worked around. mfq_ld_solve CRASHES ON ITS
 * OWN DEFAULT ARGUMENTS whenever there is more than one regime:
 *
 *     Q = {[-2 2;1 -1],[-3 3;2 -2]}; R = {diag([1 -1]), diag([0.5 -2])};
 *     S = {zeros(2), zeros(2)};  mfq_ld_solve(Q,R,S,[1 3])
 *     Error: The logical indices contain a true value outside of the array
 *     bounds.  SecondOrderLevelDependentFluidSolve line 202.
 *
 * Line 24 defaults boundaryL to zeros(K,N), a K x N MATRIX, while the
 * documented contract (line 8) and every use site want ONE ITEM PER BACKGROUND
 * STATE, a length-N vector; line 202 then indexes ix = 1:N with that K x N
 * logical mask. boundaryU inherits it. K = 1 survives only because zeros(1,N)
 * happens to be the right shape. Passing explicit length-N vectors works and
 * gives a correct answer. This port takes the boundary flags as length-N
 * vectors, which is the documented contract, and defaults them to reflective,
 * so it does not reproduce the crash.
 *
 * ARITHMETIC. Templated on T and gated on num_traits<T>::has_transcendental:
 * the QBD iteration is tolerance-terminated and expm is a Pade approximation.
 * As in mfq_ld_distr, the ONLY eigenvalue computation is the branch test inside
 * the normalization's integral of a matrix exponential, which selects between
 * two algebraically equivalent formulas and never supplies a value that reaches
 * the result, so Real instantiation is honest here. Contrast mfq_multiregime,
 * where the Schur factors ARE the basis of the answer and double is the ceiling.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/mfq_ld_distr.h"
#include "line/api/mam/mfq_ld_mean.h"
#include "line/api/mam/mfq_solve.h"
#include "line/api/mam/qbd_r.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Boundary behaviour of one background state at a reflecting level. */
enum class FluidBoundary { Reflective = 0, Absorbing = 1 };

namespace ld_solve_detail {

/** Rows ri, columns ci of A. */
template <class T>
Matrix<T> pick(const Matrix<T>& A, const std::vector<std::size_t>& ri,
               const std::vector<std::size_t>& ci) {
    Matrix<T> B(ri.size(), ci.size(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < ri.size(); ++i)
        for (std::size_t j = 0; j < ci.size(); ++j) B(i, j) = A(ri[i], ci[j]);
    return B;
}

/** Whether v contains x. */
inline bool has(const std::vector<std::size_t>& v, std::size_t x) {
    for (std::size_t y : v)
        if (y == x) return true;
    return false;
}

/**
 * The R of the reference's QBD_CR(Bm, Lm, Fm), which satisfies
 * Fm + R Lm + R^2 Bm = 0 and is therefore the port's qbd_R under the direct
 * mapping B = Bm, L = Lm, F = Fm. See the file header for how that was
 * established and for the reading it displaced.
 */
template <class T>
Matrix<T> qbd_cr_R(const Matrix<T>& Bm, const Matrix<T>& Lm, const Matrix<T>& Fm, const T& tol) {
    return qbd_R(Bm, Lm, Fm, 200000u, tol);
}

}  // namespace ld_solve_detail

/**
 * Solve a first- or second-order level-dependent fluid queue.
 *
 * @param Q         per-regime generators, K of them
 * @param R         per-regime DIAGONAL drift matrices, K of them
 * @param S         per-regime DIAGONAL variance matrices; all zero = first order
 * @param Thr       the K thresholds
 * @param boundaryL per background state, the behaviour at the lower boundary;
 *                  empty means every state reflective
 * @param boundaryU likewise at the upper boundary; empty means the same as
 *                  boundaryL, as in the reference
 * @param Qt        boundary generators, K+1 of them; empty means
 *                  {Q[0], ..., Q[K-1], Q[K-1]}, a single entry is replicated
 * @param prec      tolerance for the state classification and the QBD solves
 */
template <class T>
LevelDependentFluidBlocks<T> mfq_ld_solve(const std::vector<Matrix<T>>& Q,
                                          const std::vector<Matrix<T>>& R,
                                          const std::vector<Matrix<T>>& S,
                                          const std::vector<T>& Thr,
                                          const std::vector<FluidBoundary>& boundaryL,
                                          const std::vector<FluidBoundary>& boundaryU,
                                          const std::vector<Matrix<T>>& Qt, const T& prec) {
    static_assert(num_traits<T>::has_transcendental,
                  "mfq_ld_solve runs a tolerance-terminated QBD iteration and evaluates expm");
    using namespace ld_solve_detail;
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const std::size_t K = Thr.size();
    if (K == 0) throw InputError("mfq_ld_solve: at least one regime is required");
    if (Q.size() != K || R.size() != K || S.size() != K)
        throw InputError("mfq_ld_solve: Q, R and S must have one entry per regime");
    const std::size_t N = Q[0].rows();
    if (N == 0) throw InputError("mfq_ld_solve: the background chain is empty");

    std::vector<FluidBoundary> bL = boundaryL, bU = boundaryU;
    if (bL.empty()) bL.assign(N, FluidBoundary::Reflective);
    if (bU.empty()) bU = bL;
    if (bL.size() != N || bU.size() != N)
        throw InputError(
            "mfq_ld_solve: the boundary flags are one per BACKGROUND STATE, a vector of length N");
    std::vector<Matrix<T>> Qtk = Qt;
    if (Qtk.empty()) {
        for (std::size_t k = 0; k < K; ++k) Qtk.push_back(Q[k]);
        Qtk.push_back(Q[K - 1]);
    } else if (Qtk.size() == 1) {
        while (Qtk.size() < K + 1) Qtk.push_back(Qtk[0]);
    }
    if (Qtk.size() != K + 1) throw InputError("mfq_ld_solve: expected K+1 boundary generators");

    std::vector<T> Tv(K + 1, zero);
    for (std::size_t k = 0; k < K; ++k) Tv[k + 1] = Thr[k];

    // ---- per-regime forward and backward exponents ----
    std::vector<Matrix<T>> Sh(K);  // S/2, the form every later formula uses
    std::vector<Matrix<T>> KF(K), KB(K), cloF(K), cloB(K);
    std::vector<std::size_t> Np(K, 0), Nn(K, 0), Ns(K, 0), NbF(K, 0), NbB(K, 0);
    std::vector<std::vector<std::size_t>> vixp(K), vixn(K), vix0(K), vixs(K);
    for (std::size_t k = 0; k < K; ++k) {
        if (Q[k].rows() != N || R[k].rows() != N || S[k].rows() != N)
            throw InputError("mfq_ld_solve: a regime has the wrong order");
        Sh[k] = Matrix<T>(N, N, zero);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) Sh[k](i, j) = S[k](i, j) / two;

        std::vector<std::size_t> ix0, ixn0;
        for (std::size_t i = 0; i < N; ++i) {
            if (num_abs(R[k](i, i)) <= prec && Sh[k](i, i) <= prec) ix0.push_back(i);
            else ixn0.push_back(i);
        }
        const std::size_t Nzero = ix0.size(), Nnz = ixn0.size();
        if (Nnz == 0)
            throw InputError("mfq_ld_solve: a regime has neither drift nor variance in any state");

        Matrix<T> Qv = pick(Q[k], ixn0, ixn0);
        Matrix<T> Wz(Nnz, Nzero, zero);  // Qn0 inv(-Q00), reused for the closing matrices
        if (Nzero > 0) {
            Matrix<T> negQ00 = pick(Q[k], ix0, ix0);
            for (std::size_t i = 0; i < Nzero; ++i)
                for (std::size_t j = 0; j < Nzero; ++j) negQ00(i, j) = -negQ00(i, j);
            Wz = matmul(pick(Q[k], ixn0, ix0), inverse(negQ00));
            Qv = mfq_detail::add(Qv, matmul(Wz, pick(Q[k], ix0, ixn0)));
        }
        const Matrix<T> Rv = pick(R[k], ixn0, ixn0);
        const Matrix<T> Sv = pick(Sh[k], ixn0, ixn0);

        std::vector<std::size_t> ixp, ixn, ixs;
        for (std::size_t i = 0; i < Nnz; ++i) {
            if (Sv(i, i) > prec) ixs.push_back(i);
            else if (Rv(i, i) > prec) ixp.push_back(i);
            else if (Rv(i, i) < -prec) ixn.push_back(i);
        }
        Np[k] = ixp.size();
        Nn[k] = ixn.size();
        Ns[k] = ixs.size();

        // The discriminant of the second-order balance equation, shared by both
        // directions: Rv^2 - 2 (2 Sv) Qv on the variance states.
        std::vector<T> discr(ixs.size(), zero);
        for (std::size_t i = 0; i < ixs.size(); ++i) {
            const std::size_t q = ixs[i];
            discr[i] = Rv(q, q) * Rv(q, q) -
                       two * (two * Sv(q, q)) * Qv(q, q);
        }

        // ---- FORWARD ----
        {
            T c = one;
            for (std::size_t q : ixp) {
                const T v = -Qv(q, q) / Rv(q, q);
                if (v > c) c = v;
            }
            for (std::size_t i = 0; i < ixs.size(); ++i)
                if (discr[i] > zero) {
                    const std::size_t q = ixs[i];
                    const T v = (-Rv(q, q) + sqrt(discr[i])) / (two * Sv(q, q));
                    if (v > c) c = v;
                }
            std::vector<std::size_t> ixbF = ixs;
            ixbF.insert(ixbF.end(), ixp.begin(), ixp.end());
            NbF[k] = ixbF.size();
            const std::size_t nb = NbF[k], nn = Nn[k];
            Matrix<T> Bm(nb + nn, nb + nn, zero), Lm(nb + nn, nb + nn, zero),
                Fm(nb + nn, nb + nn, zero);
            const Matrix<T> Sbb = pick(Sv, ixbF, ixbF), Rbb = pick(Rv, ixbF, ixbF);
            const Matrix<T> Qbb = pick(Qv, ixbF, ixbF), Qbn = pick(Qv, ixbF, ixn);
            const Matrix<T> Qnb = pick(Qv, ixn, ixbF), Qnn = pick(Qv, ixn, ixn);
            const Matrix<T> Rnn = pick(Rv, ixn, ixn);
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t j = 0; j < nb; ++j) {
                    Bm(i, j) = c * Sbb(i, j);
                    Lm(i, j) = -Rbb(i, j) - two * c * Sbb(i, j);
                    Fm(i, j) = Qbb(i, j) / c + c * Sbb(i, j) + Rbb(i, j);
                }
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t j = 0; j < nn; ++j) Fm(i, nb + j) = Qbn(i, j) / c;
            for (std::size_t i = 0; i < nn; ++i) {
                for (std::size_t j = 0; j < nb; ++j) Lm(nb + i, j) = Qnb(i, j) / c;
                for (std::size_t j = 0; j < nn; ++j) {
                    Bm(nb + i, nb + j) = -Rnn(i, j);
                    Lm(nb + i, nb + j) = Qnn(i, j) / c + Rnn(i, j);
                }
            }
            const Matrix<T> QR = qbd_cr_R(Bm, Lm, Fm, prec);
            KF[k] = Matrix<T>(nb, nb, zero);
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t j = 0; j < nb; ++j)
                    KF[k](i, j) = (QR(i, j) - (i == j ? one : zero)) * c;
            // clovF has a column per non-zero-drift state: identity on ixbF,
            // Psi on ixn.
            Matrix<T> clov(nb, Nnz, zero);
            for (std::size_t i = 0; i < nb; ++i) clov(i, ixbF[i]) = one;
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t j = 0; j < nn; ++j) clov(i, ixn[j]) = QR(i, nb + j);
            cloF[k] = Matrix<T>(nb, N, zero);
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t j = 0; j < Nnz; ++j) cloF[k](i, ixn0[j]) = clov(i, j);
            if (Nzero > 0) {
                const Matrix<T> z = matmul(clov, Wz);
                for (std::size_t i = 0; i < nb; ++i)
                    for (std::size_t j = 0; j < Nzero; ++j) cloF[k](i, ix0[j]) = z(i, j);
            }
        }

        // ---- BACKWARD ----
        {
            T c = one;
            for (std::size_t q : ixn) {
                const T v = -Qv(q, q) / (-Rv(q, q));
                if (v > c) c = v;
            }
            for (std::size_t i = 0; i < ixs.size(); ++i)
                if (discr[i] > zero) {
                    const std::size_t q = ixs[i];
                    const T v = (Rv(q, q) + sqrt(discr[i])) / (two * Sv(q, q));
                    if (v > c) c = v;
                }
            std::vector<std::size_t> ixbB = ixs;
            ixbB.insert(ixbB.end(), ixn.begin(), ixn.end());
            NbB[k] = ixbB.size();
            const std::size_t nb = NbB[k], np = Np[k];
            Matrix<T> Bm(nb + np, nb + np, zero), Lm(nb + np, nb + np, zero),
                Fm(nb + np, nb + np, zero);
            const Matrix<T> Sbb = pick(Sv, ixbB, ixbB), Rbb = pick(Rv, ixbB, ixbB);
            const Matrix<T> Qbb = pick(Qv, ixbB, ixbB), Qbp = pick(Qv, ixbB, ixp);
            const Matrix<T> Qpb = pick(Qv, ixp, ixbB), Qpp = pick(Qv, ixp, ixp);
            const Matrix<T> Rpp = pick(Rv, ixp, ixp);
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t j = 0; j < nb; ++j) {
                    Bm(i, j) = c * Sbb(i, j);
                    Lm(i, j) = Rbb(i, j) - two * c * Sbb(i, j);
                    Fm(i, j) = Qbb(i, j) / c + c * Sbb(i, j) - Rbb(i, j);
                }
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t j = 0; j < np; ++j) Fm(i, nb + j) = Qbp(i, j) / c;
            for (std::size_t i = 0; i < np; ++i) {
                for (std::size_t j = 0; j < nb; ++j) Lm(nb + i, j) = Qpb(i, j) / c;
                for (std::size_t j = 0; j < np; ++j) {
                    Bm(nb + i, nb + j) = Rpp(i, j);
                    Lm(nb + i, nb + j) = Qpp(i, j) / c - Rpp(i, j);
                }
            }
            const Matrix<T> QR = qbd_cr_R(Bm, Lm, Fm, prec);
            KB[k] = Matrix<T>(nb, nb, zero);
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t j = 0; j < nb; ++j)
                    KB[k](i, j) = (QR(i, j) - (i == j ? one : zero)) * c;
            Matrix<T> clov(nb, Nnz, zero);
            for (std::size_t i = 0; i < nb; ++i) clov(i, ixbB[i]) = one;
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t j = 0; j < np; ++j) clov(i, ixp[j]) = QR(i, nb + j);
            cloB[k] = Matrix<T>(nb, N, zero);
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t j = 0; j < Nnz; ++j) cloB[k](i, ixn0[j]) = clov(i, j);
            if (Nzero > 0) {
                const Matrix<T> z = matmul(clov, Wz);
                for (std::size_t i = 0; i < nb; ++i)
                    for (std::size_t j = 0; j < Nzero; ++j) cloB[k](i, ix0[j]) = z(i, j);
            }
        }

        // Index sets lifted back to the original state numbering.
        for (std::size_t q : ixp) vixp[k].push_back(ixn0[q]);
        for (std::size_t q : ixn) vixn[k].push_back(ixn0[q]);
        for (std::size_t q : ixs) vixs[k].push_back(ixn0[q]);
        vix0[k] = ix0;
    }

    // ---- boundary system ----
    // Unknown layout: masses[0], then per regime iniF, iniB, masses[k+1].
    std::vector<std::size_t> pp;
    pp.push_back(0);
    pp.push_back(N);
    for (std::size_t k = 0; k < K; ++k) {
        pp.push_back(pp.back() + NbF[k]);
        pp.push_back(pp.back() + NbB[k]);
        pp.push_back(pp.back() + N);
    }
    // pp[0] = masses[0], pp[1] = iniF[0], pp[2] = iniB[0], pp[3] = masses[1], ...
    const std::size_t Neq = pp.back();
    Matrix<T> M(Neq, Neq, zero);

    auto rowF = [&](std::size_t k) { return pp[1 + 3 * k]; };
    auto rowB = [&](std::size_t k) { return pp[2 + 3 * k]; };
    auto rowM = [&](std::size_t k) { return k == 0 ? std::size_t(0) : pp[3 * k]; };

    // Flux conservation at every threshold.
    for (std::size_t k = 0; k <= K; ++k) {
        const std::size_t col = k * N;
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) M(rowM(k) + i, col + j) = -Qtk[k](i, j);
        if (k > 0) {
            const std::size_t g = k - 1;
            const Matrix<T> E = expm(KF[g], T(Tv[k] - Tv[k - 1]));
            // expm(KF Tk) (-cloF R + KF cloF S)
            const Matrix<T> aa = mfq_detail::sub(
                matmul(matmul(E, matmul(KF[g], cloF[g])), Sh[g]),
                matmul(matmul(E, cloF[g]), R[g]));
            for (std::size_t i = 0; i < NbF[g]; ++i)
                for (std::size_t j = 0; j < N; ++j) M(rowF(g) + i, col + j) = aa(i, j);
            const Matrix<T> bb = mfq_detail::sub(
                mfq_detail::scale(matmul(cloB[g], R[g]), T(-one)),
                matmul(matmul(KB[g], cloB[g]), Sh[g]));
            for (std::size_t i = 0; i < NbB[g]; ++i)
                for (std::size_t j = 0; j < N; ++j) M(rowB(g) + i, col + j) = bb(i, j);
        }
        if (k < K) {
            const Matrix<T> cc = mfq_detail::sub(matmul(cloF[k], R[k]),
                                                 matmul(matmul(KF[k], cloF[k]), Sh[k]));
            for (std::size_t i = 0; i < NbF[k]; ++i)
                for (std::size_t j = 0; j < N; ++j) M(rowF(k) + i, col + j) = cc(i, j);
            const Matrix<T> E = expm(KB[k], T(Tv[k + 1] - Tv[k]));
            const Matrix<T> dd = matmul(
                E, mfq_detail::add(matmul(cloB[k], R[k]), matmul(matmul(KB[k], cloB[k]), Sh[k])));
            for (std::size_t i = 0; i < NbB[k]; ++i)
                for (std::size_t j = 0; j < N; ++j) M(rowB(k) + i, col + j) = dd(i, j);
        }
    }

    // Boundary and continuity conditions.
    std::size_t col = (K + 1) * N;
    for (std::size_t k = 0; k <= K; ++k) {
        if (k == 0) {
            // No mass in an up-drift state, nor in a reflective variance state.
            std::vector<std::size_t> ixr0;
            for (std::size_t q : vixs[0])
                if (bL[q] == FluidBoundary::Reflective) ixr0.push_back(q);
            std::vector<std::size_t> sel = vixp[0];
            sel.insert(sel.end(), ixr0.begin(), ixr0.end());
            for (std::size_t i = 0; i < sel.size(); ++i) M(rowM(0) + sel[i], col + i) = one;
            col += sel.size();
            // Zero density in an absorbing variance state.
            std::vector<std::size_t> ixa0;
            for (std::size_t q : vixs[0])
                if (bL[q] == FluidBoundary::Absorbing) ixa0.push_back(q);
            const Matrix<T> pdfB = matmul(expm(KB[0], T(Tv[1] - Tv[0])), cloB[0]);
            for (std::size_t i = 0; i < ixa0.size(); ++i) {
                for (std::size_t r = 0; r < NbF[0]; ++r)
                    M(rowF(0) + r, col + i) = cloF[0](r, ixa0[i]);
                for (std::size_t r = 0; r < NbB[0]; ++r)
                    M(rowB(0) + r, col + i) = pdfB(r, ixa0[i]);
            }
            col += ixa0.size();
        } else if (k == K) {
            std::vector<std::size_t> ixrB;
            for (std::size_t q : vixs[K - 1])
                if (bU[q] == FluidBoundary::Reflective) ixrB.push_back(q);
            std::vector<std::size_t> sel = vixn[K - 1];
            sel.insert(sel.end(), ixrB.begin(), ixrB.end());
            for (std::size_t i = 0; i < sel.size(); ++i) M(rowM(K) + sel[i], col + i) = one;
            col += sel.size();
            std::vector<std::size_t> ixaB;
            for (std::size_t q : vixs[K - 1])
                if (bU[q] == FluidBoundary::Absorbing) ixaB.push_back(q);
            const Matrix<T> pdfF =
                matmul(expm(KF[K - 1], T(Tv[K] - Tv[K - 1])), cloF[K - 1]);
            for (std::size_t i = 0; i < ixaB.size(); ++i) {
                for (std::size_t r = 0; r < NbF[K - 1]; ++r)
                    M(rowF(K - 1) + r, col + i) = pdfF(r, ixaB[i]);
                for (std::size_t r = 0; r < NbB[K - 1]; ++r)
                    M(rowB(K - 1) + r, col + i) = cloB[K - 1](r, ixaB[i]);
            }
            col += ixaB.size();
        } else {
            const std::size_t g = k - 1;  // regime below, 0-based
            // No mass except where the drift reverses across the threshold.
            std::vector<std::size_t> st0;
            for (std::size_t q = 0; q < N; ++q) {
                const bool keep = (has(vixp[g], q) && has(vixn[k], q)) || has(vix0[g], q) ||
                                  has(vix0[k], q);
                if (!keep) st0.push_back(q);
            }
            for (std::size_t i = 0; i < st0.size(); ++i) M(rowM(k) + st0[i], col + i) = one;
            col += st0.size();
            // Continuity of the second-order density across the threshold.
            std::vector<std::size_t> sts;
            for (std::size_t q = 0; q < N; ++q) {
                const bool isS = has(vixs[g], q) || has(vixs[k], q);
                const bool excl = has(vixn[k], q) || has(vixp[g], q);
                if (isS && !excl) sts.push_back(q);
            }
            Matrix<T> sqSg(N, N, zero), sqSk(N, N, zero);
            for (std::size_t q = 0; q < N; ++q) {
                sqSg(q, q) = sqrt(Sh[g](q, q));
                sqSk(q, q) = sqrt(Sh[k](q, q));
            }
            const Matrix<T> BelowF = matmul(
                matmul(expm(KF[g], T(Tv[k] - Tv[k - 1])),
                       mfq_detail::scale(cloF[g], T(-one))),
                sqSg);
            const Matrix<T> BelowB =
                matmul(mfq_detail::scale(cloB[g], T(-one)), sqSg);
            const Matrix<T> AboveF = matmul(cloF[k], sqSk);
            const Matrix<T> AboveB =
                matmul(matmul(expm(KB[k], T(Tv[k + 1] - Tv[k])), cloB[k]), sqSk);
            for (std::size_t i = 0; i < sts.size(); ++i) {
                for (std::size_t r = 0; r < NbF[g]; ++r)
                    M(rowF(g) + r, col + i) = BelowF(r, sts[i]);
                for (std::size_t r = 0; r < NbB[g]; ++r)
                    M(rowB(g) + r, col + i) = BelowB(r, sts[i]);
                for (std::size_t r = 0; r < NbF[k]; ++r)
                    M(rowF(k) + r, col + i) = AboveF(r, sts[i]);
                for (std::size_t r = 0; r < NbB[k]; ++r)
                    M(rowB(k) + r, col + i) = AboveB(r, sts[i]);
            }
            col += sts.size();
        }
    }
    if (col != Neq)
        throw NumericError(
            "mfq_ld_solve: the boundary conditions do not close the system; check the drift and "
            "variance pattern across the thresholds");

    // Normalization replaces the first flux equation: total mass plus the
    // integral of every regime's density is one.
    {
        std::vector<T> h(Neq, zero);
        for (std::size_t i = 0; i < N; ++i) h[i] = one;
        for (std::size_t k = 0; k < K; ++k) {
            Matrix<T> sF, sB;
            mfq_ld_detail::integ_exp_pair(KF[k], KB[k], T(Tv[k + 1] - Tv[k]), sF, sB);
            const Matrix<T> gF = matmul(sF, cloF[k]);
            const Matrix<T> gB = matmul(sB, cloB[k]);
            for (std::size_t i = 0; i < NbF[k]; ++i) {
                T s = zero;
                for (std::size_t j = 0; j < N; ++j) s += gF(i, j);
                h[rowF(k) + i] = s;
            }
            for (std::size_t i = 0; i < NbB[k]; ++i) {
                T s = zero;
                for (std::size_t j = 0; j < N; ++j) s += gB(i, j);
                h[rowB(k) + i] = s;
            }
            for (std::size_t i = 0; i < N; ++i) h[rowM(k + 1) + i] = one;
        }
        for (std::size_t i = 0; i < Neq; ++i) M(i, 0) = h[i];
    }

    // b M = rhs, i.e. M^T b^T = rhs^T.
    Matrix<T> Mt(Neq, Neq, zero);
    for (std::size_t i = 0; i < Neq; ++i)
        for (std::size_t j = 0; j < Neq; ++j) Mt(i, j) = M(j, i);
    std::vector<T> rhs(Neq, zero);
    rhs[0] = one;
    const std::vector<T> b = solve(Mt, rhs);

    LevelDependentFluidBlocks<T> out;
    out.Thr = Thr;
    out.masses.assign(K + 1, std::vector<T>(N, zero));
    for (std::size_t k = 0; k <= K; ++k)
        for (std::size_t j = 0; j < N; ++j) out.masses[k][j] = b[rowM(k) + j];
    for (std::size_t k = 0; k < K; ++k) {
        out.iniF.push_back(std::vector<T>(b.begin() + static_cast<long>(rowF(k)),
                                          b.begin() + static_cast<long>(rowF(k) + NbF[k])));
        out.iniB.push_back(std::vector<T>(b.begin() + static_cast<long>(rowB(k)),
                                          b.begin() + static_cast<long>(rowB(k) + NbB[k])));
        out.KF.push_back(KF[k]);
        out.KB.push_back(KB[k]);
        out.cloF.push_back(cloF[k]);
        out.cloB.push_back(cloB[k]);
    }
    return out;
}

/** mfq_ld_solve with reflective boundaries, Qt = Q and the default prec = 1e-14. */
template <class T>
LevelDependentFluidBlocks<T> mfq_ld_solve(const std::vector<Matrix<T>>& Q,
                                          const std::vector<Matrix<T>>& R,
                                          const std::vector<Matrix<T>>& S,
                                          const std::vector<T>& Thr) {
    return mfq_ld_solve(Q, R, S, Thr, std::vector<FluidBoundary>(), std::vector<FluidBoundary>(),
                        std::vector<Matrix<T>>(), T(num_traits<T>::from_double(1e-14)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MFQ_LD_SOLVE_H
