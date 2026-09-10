/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MFQ_MULTIREGIME_H
#define LINE_API_MAM_MFQ_MULTIREGIME_H

/**
 * Multi-regime FEEDBACK Markovian fluid queue: density, density derivative and
 * distribution of the fluid level.
 *
 * Port of matlab/src/api/mam/mfq_multiregime.m and the BUTools multiregime it
 * wraps, which implements H. E. Kankaya and N. Akar, "Solving Multi-Regime
 * Feedback Fluid Queues". The generator and the drift rates are regime
 * dependent, and separate FEEDBACK generators and rates govern the behaviour at
 * each threshold, which is what distinguishes this from the level-dependent
 * model behind mfq_ld_solve: there the thresholds only separate regimes, here
 * the process can behave differently while sitting exactly on one.
 *
 * METHOD. Per regime, the zero-drift states are censored out and the remaining
 * generator is rescaled by the drifts, giving A = Qbar diag(1/R). The spectrum
 * of A splits into zero, negative and positive parts, and the three
 * corresponding invariant subspaces carry the constant, the decaying-upward and
 * the decaying-downward components of the density. A block-triangularizing
 * similarity Y is built from the ordered real Schur form of A plus two
 * Sylvester solves, after which the density in regime k is
 *
 *   pi_k(x) = a0_k L0_k + an_k exp(An_k (x - T_k)) Ln_k
 *                       + ap_k exp(-Ap_k (T_{k+1} - x)) Lp_k.
 *
 * The unknown coefficients, together with the point masses at the K+1
 * thresholds, solve one linear system assembled from the reference's equations
 * (8)-(16): flow balance at each boundary, the boundary conditions that the
 * feedback rates impose on each state, and one normalization row that replaces
 * the first balance equation.
 *
 * WHY THE SCHUR FORM, AND WHY NOT EIGENVECTORS. The eigenvector basis exists
 * only when A is diagonalizable and is complex whenever the spectrum is;
 * substituting it here would be a different algorithm wearing the same name.
 * The real Schur basis is orthogonal and real for every real A, and the
 * ordering by eigenvalue sign class is what isolates the three subspaces. That
 * is why line/util/eig.h grew schur_decomposition and schur_reorder for this
 * function.
 *
 * DOUBLE ONLY, AND HONESTLY SO. Unlike mfq_ld_distr -- where the only
 * eigenvalues merely SELECT a branch and never reach a returned number -- here
 * the Schur factors Z and T are the basis in which An, Ap, L0, Ln and Lp are
 * expressed, so they enter every value the function returns. util/eig.h is
 * double-only by design (the eigenvalues of a rational matrix are algebraic,
 * not rational, and LAPACK cannot supply a multiprecision QR), so this function
 * takes and returns Matrix<double> rather than being templated. Instantiating
 * it at Real50 would advertise a precision it cannot deliver, since everything
 * downstream of the Schur step would be carrying double-accurate inputs. It is
 * registered as {Double} for that reason and no other. When the build has no
 * LAPACK it refuses by name through schur_decomposition rather than guessing a
 * basis.
 *
 * COMPARING THE OUTPUT AGAINST MATLAB. The Cdfm column at level zero is an
 * EXACT ZERO in this port and comes back from MATLAB as -6.26e-17, its solve
 * having rounded. Any comparison of that column must therefore be ABSOLUTE, not
 * relative: a relative test against a value that is exactly zero fails on a
 * difference of one ulp and reports a defect that is not there. The same holds
 * for any state carrying no mass at a threshold. This is stated here rather
 * than only in the test so that a later reader does not "tidy" the absolute
 * comparison into a relative one and manufacture a failure.
 *
 * THE SYLVESTER EQUATIONS ARE NOT THE OBSTACLE. A X + X B = C vectorizes to
 * (I kron A + B^T kron I) vec(X) = vec(C), one ordinary linear solve at these
 * block sizes, needing no Bartels-Stewart and no Schur form of its own. The
 * reference calls MATLAB's sylvester, which is Bartels-Stewart; the two compute
 * the same X, and at the sizes reached here the direct solve is not the slower
 * one.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/mfq_solve.h"
#include "line/num/number.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Return value of mfq_multiregime: one row of N per-state values per point. */
struct MultiRegimeResult {
    std::vector<std::vector<double>> pdf;   ///< density at each pdf point
    std::vector<std::vector<double>> pdfd;  ///< derivative of the density
    std::vector<std::vector<double>> cdf;   ///< P(X < p)
    std::vector<std::vector<double>> cdfm;  ///< P(X <= p)
};

namespace multiregime_detail {

/**
 * Solve the Sylvester equation A X + X B = C for X, with A of order n, B of
 * order m and C of size n x m, by the Kronecker form
 * (I_m kron A + B^T kron I_n) vec(X) = vec(C), vec being column-major.
 */
inline Matrix<double> sylvester(const Matrix<double>& A, const Matrix<double>& B,
                                const Matrix<double>& C) {
    const std::size_t n = A.rows(), m = B.rows();
    if (A.cols() != n || B.cols() != m || C.rows() != n || C.cols() != m)
        throw InputError("mfq_multiregime sylvester: blocks are not conformable");
    Matrix<double> X(n, m, 0.0);
    if (n == 0 || m == 0) return X;
    Matrix<double> K(n * m, n * m, 0.0);
    std::vector<double> rhs(n * m, 0.0);
    for (std::size_t j = 0; j < m; ++j)
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t row = j * n + i;
            for (std::size_t k = 0; k < n; ++k) K(row, j * n + k) += A(i, k);
            for (std::size_t k = 0; k < m; ++k) K(row, k * n + i) += B(k, j);
            rhs[row] = C(i, j);
        }
    const std::vector<double> x = solve(K, rhs);
    for (std::size_t j = 0; j < m; ++j)
        for (std::size_t i = 0; i < n; ++i) X(i, j) = x[j * n + i];
    return X;
}

/**
 * expm of a possibly EMPTY block. A regime whose spectrum has no negative (or
 * no positive) part leaves An (or Ap) at order zero, which is a normal state of
 * this model and not a degenerate one: it simply means the density has no
 * component decaying in that direction. util/expm.h rejects an empty matrix, so
 * the empty case is answered here with the empty matrix rather than by relaxing
 * that check for every other caller.
 */
inline Matrix<double> expm0(const Matrix<double>& A, double t) {
    if (A.rows() == 0) return Matrix<double>(0, 0, 0.0);
    return expm(A, t);
}

/** Sub-block A(r0..r0+nr-1, c0..c0+nc-1); an empty request gives an empty matrix. */
inline Matrix<double> blk(const Matrix<double>& A, std::size_t r0, std::size_t c0, std::size_t nr,
                          std::size_t nc) {
    Matrix<double> B(nr, nc, 0.0);
    for (std::size_t i = 0; i < nr; ++i)
        for (std::size_t j = 0; j < nc; ++j) B(i, j) = A(r0 + i, c0 + j);
    return B;
}

/** Rows of A selected by ri, columns by ci. */
inline Matrix<double> pick(const Matrix<double>& A, const std::vector<std::size_t>& ri,
                           const std::vector<std::size_t>& ci) {
    Matrix<double> B(ri.size(), ci.size(), 0.0);
    for (std::size_t i = 0; i < ri.size(); ++i)
        for (std::size_t j = 0; j < ci.size(); ++j) B(i, j) = A(ri[i], ci[j]);
    return B;
}

/** The per-regime quantities the assembly and the evaluation both need. */
struct Regime {
    Matrix<double> An, Ap;      ///< negative and positive spectral blocks
    Matrix<double> L0, Ln, Lp;  ///< the three closing matrices, each ? x N
    Matrix<double> M0, MT, Mi;  ///< regime transfer matrices at 0, at T and integrated
    std::size_t nzero = 0, nneg = 0, npos = 0;
};

}  // namespace multiregime_detail

/**
 * Multi-regime feedback fluid queue.
 *
 * @param Q         per-regime generators, K of them; a single entry is
 *                  replicated across all regimes, as in the reference
 * @param R         per-regime drift RATE VECTORS of length N (not matrices)
 * @param Qt        boundary generators, K+1 of them; empty means
 *                  {Q[0], Q[0], ..., Q[K-1]}, a single entry is replicated
 * @param Rt        boundary rate vectors, K+1 of them; empty means
 *                  {R[0], R[0], ..., R[K-1]}
 * @param Thr       the K thresholds
 * @param pdfpoints levels at which the density and its derivative are wanted
 * @param cdfpoints levels at which the distribution is wanted
 */
inline MultiRegimeResult mfq_multiregime(const std::vector<Matrix<double>>& Q,
                                         const std::vector<std::vector<double>>& R,
                                         const std::vector<Matrix<double>>& Qt,
                                         const std::vector<std::vector<double>>& Rt,
                                         const std::vector<double>& Thr,
                                         const std::vector<double>& pdfpoints,
                                         const std::vector<double>& cdfpoints) {
    using namespace multiregime_detail;
    const std::size_t K = R.size();
    if (K == 0) throw InputError("mfq_multiregime: at least one regime is required");
    if (Thr.size() != K) throw InputError("mfq_multiregime: one threshold per regime is required");
    const std::size_t N = R[0].size();
    if (N == 0) throw InputError("mfq_multiregime: the background chain is empty");

    // Convenience expansions, exactly the reference's.
    std::vector<Matrix<double>> Qk = Q;
    if (Qk.size() == 1)
        for (std::size_t k = 1; k < K; ++k) Qk.push_back(Qk[0]);
    if (Qk.size() != K) throw InputError("mfq_multiregime: expected one generator per regime");
    std::vector<Matrix<double>> Qtk = Qt;
    if (Qtk.empty()) {
        Qtk.push_back(Qk[0]);
        for (std::size_t k = 0; k < K; ++k) Qtk.push_back(Qk[k]);
    } else if (Qtk.size() == 1) {
        for (std::size_t k = 1; k < K + 1; ++k) Qtk.push_back(Qtk[0]);
    }
    if (Qtk.size() != K + 1)
        throw InputError("mfq_multiregime: expected K+1 boundary generators");
    std::vector<std::vector<double>> Rtk = Rt;
    if (Rtk.empty()) {
        Rtk.push_back(R[0]);
        for (std::size_t k = 0; k < K; ++k) Rtk.push_back(R[k]);
    }
    if (Rtk.size() != K + 1)
        throw InputError("mfq_multiregime: expected K+1 boundary rate vectors");
    for (std::size_t k = 0; k < K; ++k) {
        if (R[k].size() != N || Qk[k].rows() != N || Qk[k].cols() != N)
            throw InputError("mfq_multiregime: a regime's Q or R has the wrong order");
    }

    // Thresholds with the implicit zero in front.
    std::vector<double> Tv(K + 1, 0.0);
    for (std::size_t k = 0; k < K; ++k) Tv[k + 1] = Thr[k];

    // ---- per-regime transfer matrices ----
    std::vector<Regime> reg(K);
    std::vector<std::size_t> Nnz(K, 0), Npos(K + 1, 0);
    for (std::size_t k = 0; k < K; ++k) {
        std::vector<std::size_t> zix, nzix;
        for (std::size_t i = 0; i < N; ++i) (R[k][i] == 0.0 ? zix : nzix).push_back(i);
        const std::size_t Nn = nzix.size(), Nz = zix.size();
        if (Nn == 0)
            throw InputError("mfq_multiregime: a regime has no state with a non-zero drift");

        // Censor the zero-drift states, then rescale by the drifts.
        Matrix<double> Qnk = pick(Qk[k], nzix, nzix);
        Matrix<double> Wz(Nn, Nz, 0.0);  // Q(nz,z) inv(-Q(z,z)), reused for the L blocks
        if (Nz > 0) {
            Matrix<double> negQzz = pick(Qk[k], zix, zix);
            for (std::size_t i = 0; i < Nz; ++i)
                for (std::size_t j = 0; j < Nz; ++j) negQzz(i, j) = -negQzz(i, j);
            Matrix<double> iQzz(Nz, Nz, 0.0);
            try {
                iQzz = inverse(negQzz);
            } catch (const NumericError&) {
                throw NumericError(
                    "mfq_multiregime: the zero-drift states of a regime form a closed set, so the "
                    "fluid can be trapped at a constant level and the regime has no stationary "
                    "law");
            }
            Wz = matmul(pick(Qk[k], nzix, zix), iQzz);
            Qnk = mfq_detail::add(Qnk, matmul(Wz, pick(Qk[k], zix, nzix)));
        }
        Matrix<double> A(Nn, Nn, 0.0);
        for (std::size_t i = 0; i < Nn; ++i)
            for (std::size_t j = 0; j < Nn; ++j) A(i, j) = Qnk(i, j) / R[k][nzix[j]];

        // ordered real Schur form rationale: see _kb/03-api-layer.md (cpp port notes: mam)
        const RealSchur s0 = schur_decomposition(A);
        std::vector<double> key(Nn, 0.0);
        std::size_t nzero = 0;
        for (std::size_t i = 0; i < Nn; ++i) {
            const double d = s0.T(i, i);
            key[i] = (std::fabs(d) < 1e-10 ? 5.0 : 0.0) + (d < 0.0 ? 2.0 : 0.0) +
                     (d > 0.0 ? 1.0 : 0.0);
            if (std::fabs(d) < 1e-10) ++nzero;
        }
        // A 2 x 2 block is one complex pair and must not be split; give both
        // rows the key of the first, as MATLAB requires of a select vector.
        for (std::size_t i = 0; i + 1 < Nn; ++i)
            if (s0.T(i + 1, i) != 0.0) key[i + 1] = key[i];
        const RealSchur s = schur_reorder(s0, key);

        std::size_t nneg = 0, npos = 0;
        for (std::size_t i = nzero; i < Nn; ++i) {
            if (s.T(i, i) < 0.0) ++nneg;
            else if (s.T(i, i) > 0.0) ++npos;
        }
        if (nzero + nneg + npos != Nn)
            throw NumericError(
                "mfq_multiregime: the spectrum of a regime could not be split into zero, negative "
                "and positive parts");

        // X1 solves 0 X1 + X1 (-D22) = D12, X2 solves Dnn X2 - X2 Dpp = Dnp.
        const Matrix<double> D22 = blk(s.T, nzero, nzero, Nn - nzero, Nn - nzero);
        Matrix<double> negD22 = D22;
        for (std::size_t i = 0; i < negD22.rows(); ++i)
            for (std::size_t j = 0; j < negD22.cols(); ++j) negD22(i, j) = -negD22(i, j);
        const Matrix<double> X1 =
            sylvester(Matrix<double>(nzero, nzero, 0.0), negD22,
                      blk(s.T, 0, nzero, nzero, Nn - nzero));
        Matrix<double> negDpp = blk(s.T, nzero + nneg, nzero + nneg, npos, npos);
        for (std::size_t i = 0; i < npos; ++i)
            for (std::size_t j = 0; j < npos; ++j) negDpp(i, j) = -negDpp(i, j);
        const Matrix<double> X2 = sylvester(blk(s.T, nzero, nzero, nneg, nneg), negDpp,
                                            blk(s.T, nzero, nzero + nneg, nneg, npos));

        // Y = Z [I -X1; 0 I] [I 0 0; 0 I -X2; 0 0 I].
        Matrix<double> B1 = eye<double>(Nn);
        for (std::size_t i = 0; i < nzero; ++i)
            for (std::size_t j = 0; j < Nn - nzero; ++j) B1(i, nzero + j) = -X1(i, j);
        Matrix<double> B2 = eye<double>(Nn);
        for (std::size_t i = 0; i < nneg; ++i)
            for (std::size_t j = 0; j < npos; ++j) B2(nzero + i, nzero + nneg + j) = -X2(i, j);
        const Matrix<double> Y = matmul(matmul(s.Z, B1), B2);
        const Matrix<double> iY = inverse(Y);
        const Matrix<double> At = matmul(matmul(iY, A), Y);

        Regime& rg = reg[k];
        rg.nzero = nzero;
        rg.nneg = nneg;
        rg.npos = npos;
        rg.An = blk(At, nzero, nzero, nneg, nneg);
        rg.Ap = blk(At, nzero + nneg, nzero + nneg, npos, npos);

        // L blocks: the non-zero-drift columns are the rows of iY, the
        // zero-drift ones are routed through the censoring operator.
        const Matrix<double> iY0 = blk(iY, 0, 0, nzero, Nn);
        const Matrix<double> iYn = blk(iY, nzero, 0, nneg, Nn);
        const Matrix<double> iYp = blk(iY, nzero + nneg, 0, npos, Nn);
        const Matrix<double> Z0 = (Nz > 0) ? matmul(iY0, Wz) : Matrix<double>(nzero, 0, 0.0);
        const Matrix<double> Zn = (Nz > 0) ? matmul(iYn, Wz) : Matrix<double>(nneg, 0, 0.0);
        const Matrix<double> Zp = (Nz > 0) ? matmul(iYp, Wz) : Matrix<double>(npos, 0, 0.0);
        rg.L0 = Matrix<double>(nzero, N, 0.0);
        rg.Ln = Matrix<double>(nneg, N, 0.0);
        rg.Lp = Matrix<double>(npos, N, 0.0);
        for (std::size_t j = 0; j < Nn; ++j) {
            for (std::size_t i = 0; i < nzero; ++i) rg.L0(i, nzix[j]) = iY0(i, j);
            for (std::size_t i = 0; i < nneg; ++i) rg.Ln(i, nzix[j]) = iYn(i, j);
            for (std::size_t i = 0; i < npos; ++i) rg.Lp(i, nzix[j]) = iYp(i, j);
        }
        for (std::size_t j = 0; j < Nz; ++j) {
            for (std::size_t i = 0; i < nzero; ++i) rg.L0(i, zix[j]) = Z0(i, j);
            for (std::size_t i = 0; i < nneg; ++i) rg.Ln(i, zix[j]) = Zn(i, j);
            for (std::size_t i = 0; i < npos; ++i) rg.Lp(i, zix[j]) = Zp(i, j);
        }

        const double Tk = Tv[k + 1] - Tv[k];
        Matrix<double> negApT = rg.Ap;
        for (std::size_t i = 0; i < npos; ++i)
            for (std::size_t j = 0; j < npos; ++j) negApT(i, j) = -negApT(i, j);
        const Matrix<double> EAn = expm0(rg.An, Tk);
        const Matrix<double> EAp = expm0(negApT, Tk);  // exp(-Ap Tk)

        // M0 = [L0; Ln; exp(-Ap Tk) Lp], MT = [L0; exp(An Tk) Ln; Lp].
        rg.M0 = Matrix<double>(Nn, N, 0.0);
        rg.MT = Matrix<double>(Nn, N, 0.0);
        rg.Mi = Matrix<double>(Nn, N, 0.0);
        const Matrix<double> EApLp = matmul(EAp, rg.Lp);
        const Matrix<double> EAnLn = matmul(EAn, rg.Ln);
        for (std::size_t i = 0; i < nzero; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                rg.M0(i, j) = rg.L0(i, j);
                rg.MT(i, j) = rg.L0(i, j);
                rg.Mi(i, j) = Tk * rg.L0(i, j);
            }
        for (std::size_t i = 0; i < nneg; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                rg.M0(nzero + i, j) = rg.Ln(i, j);
                rg.MT(nzero + i, j) = EAnLn(i, j);
            }
        for (std::size_t i = 0; i < npos; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                rg.M0(nzero + nneg + i, j) = EApLp(i, j);
                rg.MT(nzero + nneg + i, j) = rg.Lp(i, j);
            }
        // Mi = [Tk L0; (-An)^-1 (I - exp(An Tk)) Ln; Ap^-1 (I - exp(-Ap Tk)) Lp].
        if (nneg > 0) {
            Matrix<double> negAn = rg.An;
            for (std::size_t i = 0; i < nneg; ++i)
                for (std::size_t j = 0; j < nneg; ++j) negAn(i, j) = -negAn(i, j);
            Matrix<double> ImE = eye<double>(nneg);
            for (std::size_t i = 0; i < nneg; ++i)
                for (std::size_t j = 0; j < nneg; ++j) ImE(i, j) -= EAn(i, j);
            const Matrix<double> blkn = matmul(matmul(inverse(negAn), ImE), rg.Ln);
            for (std::size_t i = 0; i < nneg; ++i)
                for (std::size_t j = 0; j < N; ++j) rg.Mi(nzero + i, j) = blkn(i, j);
        }
        if (npos > 0) {
            Matrix<double> ImE = eye<double>(npos);
            for (std::size_t i = 0; i < npos; ++i)
                for (std::size_t j = 0; j < npos; ++j) ImE(i, j) -= EAp(i, j);
            const Matrix<double> blkp = matmul(matmul(inverse(rg.Ap), ImE), rg.Lp);
            for (std::size_t i = 0; i < npos; ++i)
                for (std::size_t j = 0; j < N; ++j) rg.Mi(nzero + nneg + i, j) = blkp(i, j);
        }

        Nnz[k] = Nn;
        Npos[k + 1] = Npos[k] + Nn;
    }

    // linear system layout: see _kb/03-api-layer.md (cpp port notes: mam)
    const std::size_t d = (K + 1) * N;
    const std::size_t Neq = d + Npos[K];
    Matrix<double> M(Neq, Neq, 0.0);
    std::size_t p = 0;

    // eq. (12): balance at level zero.
    for (std::size_t j = 0; j < N; ++j)
        for (std::size_t i = 0; i < N; ++i) M(i, p + j) = -Qtk[0](i, j);
    for (std::size_t j = 0; j < N; ++j)
        for (std::size_t i = 0; i < Nnz[0]; ++i)
            M(d + Npos[0] + i, p + j) = reg[0].M0(i, j) * R[0][j];
    p += N;
    // eq. (13): balance at each interior threshold.
    for (std::size_t k = 0; k + 1 < K; ++k) {
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t i = 0; i < N; ++i) M((k + 1) * N + i, p + j) = -Qtk[k + 1](i, j);
        for (std::size_t j = 0; j < N; ++j) {
            for (std::size_t i = 0; i < Nnz[k + 1]; ++i)
                M(d + Npos[k + 1] + i, p + j) = reg[k + 1].M0(i, j) * R[k + 1][j];
            for (std::size_t i = 0; i < Nnz[k]; ++i)
                M(d + Npos[k] + i, p + j) = -reg[k].MT(i, j) * R[k][j];
        }
        p += N;
    }
    // eq. (16): balance at the top threshold.
    for (std::size_t j = 0; j < N; ++j)
        for (std::size_t i = 0; i < N; ++i) M(K * N + i, p + j) = -Qtk[K](i, j);
    for (std::size_t j = 0; j < N; ++j)
        for (std::size_t i = 0; i < Nnz[K - 1]; ++i)
            M(d + Npos[K - 1] + i, p + j) = -reg[K - 1].MT(i, j) * R[K - 1][j];
    p += N;

    // eq. (8): no mass at level zero in an up-drift state.
    for (std::size_t m = 0; m < N; ++m)
        if (R[0][m] > 0.0) M(m, p++) = 1.0;
    // eq. (11): no mass at the top threshold in a down-drift state.
    for (std::size_t m = 0; m < N; ++m)
        if (R[K - 1][m] < 0.0) M(K * N + m, p++) = 1.0;
    // eq. (9): no mass where the drift keeps its sign across a threshold.
    for (std::size_t k = 0; k + 1 < K; ++k)
        for (std::size_t m = 0; m < N; ++m)
            if ((R[k][m] > 0.0 && R[k + 1][m] > 0.0) || (R[k][m] < 0.0 && R[k + 1][m] < 0.0))
                M((k + 1) * N + m, p++) = 1.0;
    // eq. (10): nor where the drift reverses upward through a live boundary.
    for (std::size_t k = 0; k + 1 < K; ++k)
        for (std::size_t m = 0; m < N; ++m)
            if (R[k][m] < 0.0 && R[k + 1][m] > 0.0 && Rtk[k + 1][m] != 0.0)
                M((k + 1) * N + m, p++) = 1.0;
    // eq. (14): continuity of the density from below at a threshold.
    for (std::size_t k = 0; k + 1 < K; ++k)
        for (std::size_t m = 0; m < N; ++m)
            if (R[k][m] < 0.0 && Rtk[k + 1][m] >= 0.0) {
                for (std::size_t i = 0; i < Nnz[k]; ++i)
                    M(d + Npos[k] + i, p) = reg[k].MT(i, m);
                ++p;
            }
    // eq. (15): and from above.
    for (std::size_t k = 0; k + 1 < K; ++k)
        for (std::size_t m = 0; m < N; ++m)
            if (R[k + 1][m] > 0.0 && Rtk[k + 1][m] <= 0.0) {
                for (std::size_t i = 0; i < Nnz[k + 1]; ++i)
                    M(d + Npos[k + 1] + i, p) = reg[k + 1].M0(i, m);
                ++p;
            }
    if (p != Neq)
        throw NumericError(
            "mfq_multiregime: the boundary conditions do not close the system; check that the "
            "feedback rates Rt satisfy the model's continuity requirements");

    // Normalization replaces the FIRST equation: total mass plus the integral
    // of every regime's density is one.
    for (std::size_t i = 0; i < d; ++i) M(i, 0) = 1.0;
    for (std::size_t k = 0; k < K; ++k)
        for (std::size_t i = 0; i < Nnz[k]; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j) s += reg[k].Mi(i, j);
            M(d + Npos[k] + i, 0) = s;
        }

    // sol M = rhs, i.e. M^T sol^T = rhs^T.
    Matrix<double> Mt(Neq, Neq, 0.0);
    for (std::size_t i = 0; i < Neq; ++i)
        for (std::size_t j = 0; j < Neq; ++j) Mt(i, j) = M(j, i);
    std::vector<double> rhs(Neq, 0.0);
    rhs[0] = 1.0;
    const std::vector<double> sol = solve(Mt, rhs);

    // ---- extract the masses and the density coefficients ----
    std::vector<std::vector<double>> masses(K + 1, std::vector<double>(N, 0.0));
    for (std::size_t k = 0; k <= K; ++k)
        for (std::size_t j = 0; j < N; ++j) masses[k][j] = sol[k * N + j];
    std::vector<std::vector<double>> a0(K), an(K), ap(K);
    for (std::size_t k = 0; k < K; ++k) {
        const std::size_t base = d + Npos[k];
        a0[k].assign(sol.begin() + static_cast<long>(base),
                     sol.begin() + static_cast<long>(base + reg[k].nzero));
        an[k].assign(sol.begin() + static_cast<long>(base + reg[k].nzero),
                     sol.begin() + static_cast<long>(base + reg[k].nzero + reg[k].nneg));
        ap[k].assign(sol.begin() + static_cast<long>(base + reg[k].nzero + reg[k].nneg),
                     sol.begin() + static_cast<long>(base + Nnz[k]));
    }

    // The regime holding a point, by the reference's walk.
    auto regime_of = [&](double p_) -> std::size_t {
        std::size_t k = 0;
        while (k < K && p_ >= Tv[k]) ++k;
        return k - 1;
    };

    MultiRegimeResult out;
    for (double pt : pdfpoints) {
        if (pt < 0.0) throw InputError("mfq_multiregime: the evaluation points must be non-negative");
        const std::size_t k = regime_of(pt);
        Matrix<double> negAp = reg[k].Ap;
        for (std::size_t i = 0; i < reg[k].npos; ++i)
            for (std::size_t j = 0; j < reg[k].npos; ++j) negAp(i, j) = -negAp(i, j);
        const Matrix<double> En = expm0(reg[k].An, pt - Tv[k]);
        const Matrix<double> Ep = expm0(negAp, Tv[k + 1] - pt);
        std::vector<double> f(N, 0.0), fd(N, 0.0);
        {
            const std::vector<double> v0 = vecmul(a0[k], reg[k].L0);
            const std::vector<double> vn = vecmul(vecmul(an[k], En), reg[k].Ln);
            const std::vector<double> vp = vecmul(vecmul(ap[k], Ep), reg[k].Lp);
            for (std::size_t j = 0; j < N; ++j) f[j] = v0[j] + vn[j] + vp[j];
            const std::vector<double> dn =
                vecmul(vecmul(vecmul(an[k], reg[k].An), En), reg[k].Ln);
            const std::vector<double> dp =
                vecmul(vecmul(vecmul(ap[k], reg[k].Ap), Ep), reg[k].Lp);
            for (std::size_t j = 0; j < N; ++j) fd[j] = dn[j] + dp[j];
        }
        out.pdf.push_back(f);
        out.pdfd.push_back(fd);
    }

    for (double c : cdfpoints) {
        if (c < 0.0) throw InputError("mfq_multiregime: the evaluation points must be non-negative");
        std::vector<double> cres(N, 0.0), cresm(N, 0.0);
        std::size_t k = 0;
        while (k < K && c >= Tv[k]) {
            if (k > 0) {
                std::vector<double> coef;
                coef.insert(coef.end(), a0[k - 1].begin(), a0[k - 1].end());
                coef.insert(coef.end(), an[k - 1].begin(), an[k - 1].end());
                coef.insert(coef.end(), ap[k - 1].begin(), ap[k - 1].end());
                const std::vector<double> v = vecmul(coef, reg[k - 1].Mi);
                for (std::size_t j = 0; j < N; ++j) {
                    cres[j] += v[j];
                    cresm[j] += v[j];
                }
            }
            for (std::size_t j = 0; j < N; ++j) cresm[j] += masses[k][j];
            if (c > Tv[k])
                for (std::size_t j = 0; j < N; ++j) cres[j] += masses[k][j];
            ++k;
        }
        if (k == K && c == Tv[K])
            for (std::size_t j = 0; j < N; ++j) cresm[j] += masses[K][j];
        const std::size_t kk = k - 1;
        const double crem = c - Tv[kk];
        const double Tk = Tv[kk + 1] - Tv[kk];
        Matrix<double> negAn = reg[kk].An, negAp = reg[kk].Ap;
        for (std::size_t i = 0; i < reg[kk].nneg; ++i)
            for (std::size_t j = 0; j < reg[kk].nneg; ++j) negAn(i, j) = -negAn(i, j);
        for (std::size_t i = 0; i < reg[kk].npos; ++i)
            for (std::size_t j = 0; j < reg[kk].npos; ++j) negAp(i, j) = -negAp(i, j);
        std::vector<double> val(N, 0.0);
        {
            const std::vector<double> v0 = vecmul(a0[kk], reg[kk].L0);
            for (std::size_t j = 0; j < N; ++j) val[j] += v0[j] * crem;
        }
        if (reg[kk].nneg > 0) {
            Matrix<double> ImE = eye<double>(reg[kk].nneg);
            const Matrix<double> E = expm0(reg[kk].An, crem);
            for (std::size_t i = 0; i < reg[kk].nneg; ++i)
                for (std::size_t j = 0; j < reg[kk].nneg; ++j) ImE(i, j) -= E(i, j);
            const std::vector<double> v =
                vecmul(vecmul(vecmul(an[kk], inverse(negAn)), ImE), reg[kk].Ln);
            for (std::size_t j = 0; j < N; ++j) val[j] += v[j];
        }
        if (reg[kk].npos > 0) {
            const Matrix<double> Ea = expm0(negAp, Tk);
            const Matrix<double> Eb = expm0(negAp, Tk - crem);
            Matrix<double> D = Ea;
            for (std::size_t i = 0; i < reg[kk].npos; ++i)
                for (std::size_t j = 0; j < reg[kk].npos; ++j) D(i, j) -= Eb(i, j);
            const std::vector<double> v =
                vecmul(vecmul(vecmul(ap[kk], inverse(negAp)), D), reg[kk].Lp);
            for (std::size_t j = 0; j < N; ++j) val[j] += v[j];
        }
        for (std::size_t j = 0; j < N; ++j) {
            cres[j] += val[j];
            cresm[j] += val[j];
        }
        out.cdf.push_back(cres);
        out.cdfm.push_back(cresm);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MFQ_MULTIREGIME_H
