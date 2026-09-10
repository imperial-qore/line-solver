/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LIB_SMC_ETAQA_H
#define LINE_LIB_SMC_ETAQA_H

/**
 * ETAQA: the aggregated stationary vector and the queue-length moments of an
 * M/G/1-type and of a GI/M/1-type Markov chain. Port of MAMSolver's
 * `MG1_G_ETAQA.m`, `MG1_pi_ETAQA.m`, `MG1_qlen_ETAQA.m`, `GIM1_R_ETAQA.m`,
 * `GIM1_pi_ETAQA.m` and `GIM1_qlen_ETAQA.m` (`matlab/lib/thirdparty/MAMSolver`).
 *
 * WHAT ETAQA IS. A structured infinite chain has a stationary vector
 * pi = (pi_0, pi_1, pi_2, ...) with one block per level. Matrix-geometric
 * methods compute every block; ETAQA instead solves a FINITE linear system for
 * exactly three aggregates -- pi_0, pi_1 and pi* = sum_{j>=2} pi_j -- by
 * replacing the infinitely many balance equations for levels 2 and above by
 * their sum. The aggregate system is (mb + 2m) x (mb + 2m) and is EXACT: no
 * level is truncated and no tail is fitted, which is what separates ETAQA from
 * a level-truncated direct solve. Riska and Smirni, Perform. Eval. 54(2), 2003.
 *
 * The moments come the same way. `..._qlen_ETAQA` never forms pi_j either: it
 * propagates the vectors r^(k) = sum_{j>=2} j^k pi_j through a recurrence whose
 * left-hand side is one fixed m x m system, so the n-th moment costs n solves
 * of the same size regardless of how heavy the tail is.
 *
 * WHY THE LAST COLUMN IS ALWAYS DROPPED. The aggregate balance equations are
 * linearly dependent, exactly as the balance equations of any generator are, so
 * the system is rank deficient by one and the normalization e^T pi = 1 supplies
 * the missing equation. The M/G/1 side finds WHICH column to drop with a rank
 * test (a redundant one exists but is not always the last), the GI/M/1 side
 * simply drops the final column. Both then prepend the column of ones.
 *
 * DOUBLE ONLY: see the note in `lib/smc/mg1.h`. The rank test is an SVD and
 * the fundamental matrices come from LAPACK-backed and FFT-backed iterations.
 *
 * REFERENCE DEFECTS REPRODUCED VERBATIM, because these functions are the
 * numerical reference for the JAR and Python ports and a silent repair here
 * would be a divergence nobody could see:
 *
 * 1. `MG1_G_ETAQA` tests whether its input is a DTMC by setting a flag named
 *    `isdicrete` and reading one named `isdiscrete`, so the flag never changes
 *    and the continuous-to-discrete uniformization ALWAYS runs. For the
 *    generators LINE passes it that is the correct branch anyway; a genuine
 *    stochastic input would be divided by -min(diag(A1)), which is <= 0.
 *
 * 2. `MG1_qlen_ETAQA` builds `F0j` from block dega down to block 2 but indexes
 *    it from 1, so `F0j(:,(j-1)*m+1:j*m)` is the tail sum starting at j+1 and
 *    not at j. The moment is therefore the reference's moment, off-by-one
 *    index and all.
 *
 * 3. `GIM1_qlen_ETAQA` initializes `leftr_part1 = A(3)`, a SCALAR read at
 *    linear index 3 of the stacked A rather than the third block, and MATLAB
 *    then broadcasts it over the m x m accumulator. The port reproduces both
 *    the value and the broadcast, including the case where the loop that would
 *    turn it into a matrix does not run and the scalar survives as a vector.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/lib/smc/mg1.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace smc {

namespace etaqa_detail {

/** MATLAB's `x / A` for a row vector x and a square A: solves x A = b. */
inline std::vector<double> right_divide(const std::vector<double>& b, const Matrix<double>& A) {
    if (A.rows() != A.cols()) throw InputError("ETAQA: right division by a non-square matrix");
    if (b.size() != A.rows()) throw InputError("ETAQA: right division with mismatched length");
    return solve(A.transpose(), b);
}

/** Row vector times matrix. */
inline std::vector<double> vmul(const std::vector<double>& v, const Matrix<double>& A) {
    return vecmul(v, A);
}

/** Sum of the entries of a vector. */
inline double vsum(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x;
    return s;
}

/** v A e, the row vector reduced to a scalar by a matrix and a column of ones. */
inline double vmul_e(const std::vector<double>& v, const Matrix<double>& A) {
    const std::vector<double> r = vmul(v, A);
    return vsum(r);
}

/** The binomial coefficient the reference computes as a ratio of factorials. */
inline double bino(std::size_t n, std::size_t k) {
    double num = 1.0, den = 1.0;
    for (std::size_t i = 1; i <= n; ++i) num *= static_cast<double>(i);
    for (std::size_t i = 1; i <= k; ++i) den *= static_cast<double>(i);
    for (std::size_t i = 1; i <= n - k; ++i) den *= static_cast<double>(i);
    return num / den;
}

/** Horizontal concatenation of matrices with the same number of rows. */
inline Matrix<double> hcat2(const std::vector<Matrix<double>>& parts) {
    std::size_t cols = 0;
    for (const Matrix<double>& p : parts) cols += p.cols();
    const std::size_t rows = parts.empty() ? 0 : parts[0].rows();
    Matrix<double> out(rows, cols, 0.0);
    std::size_t off = 0;
    for (const Matrix<double>& p : parts) {
        for (std::size_t i = 0; i < rows; ++i)
            for (std::size_t j = 0; j < p.cols(); ++j) out(i, off + j) = p(i, j);
        off += p.cols();
    }
    return out;
}

/** Vertical concatenation of matrices with the same number of columns. */
inline Matrix<double> vcat2(const std::vector<Matrix<double>>& parts) {
    std::size_t rows = 0;
    for (const Matrix<double>& p : parts) rows += p.rows();
    const std::size_t cols = parts.empty() ? 0 : parts[0].cols();
    Matrix<double> out(rows, cols, 0.0);
    std::size_t off = 0;
    for (const Matrix<double>& p : parts) {
        for (std::size_t i = 0; i < p.rows(); ++i)
            for (std::size_t j = 0; j < cols; ++j) out(off + i, j) = p(i, j);
        off += p.rows();
    }
    return out;
}

/** Columns `[from, to)` of A. */
inline Matrix<double> cols_of(const Matrix<double>& A, std::size_t from, std::size_t to) {
    Matrix<double> out(A.rows(), to - from, 0.0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = from; j < to; ++j) out(i, j - from) = A(i, j);
    return out;
}

/** Rows `[from, to)` of A. */
inline Matrix<double> rows_of(const Matrix<double>& A, std::size_t from, std::size_t to) {
    Matrix<double> out(to - from, A.cols(), 0.0);
    for (std::size_t i = from; i < to; ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) out(i - from, j) = A(i, j);
    return out;
}

/**
 * Block i (i >= 1) of a GI/M/1-type vertical stack whose boundary occupies the
 * first `mb` rows: `M(mb+(i-1)*m+1 : mb+i*m, :)`. The blocks are m x mb, which
 * is not m x m when the boundary has a different size, so this cannot go
 * through the uniform block splitter.
 */
inline Matrix<double> vblock_after(const Matrix<double>& M, std::size_t mb, std::size_t m,
                                   std::size_t i) {
    return rows_of(M, mb + (i - 1) * m, mb + i * m);
}

/** A with column `drop` removed. */
inline Matrix<double> drop_col(const Matrix<double>& A, std::size_t drop) {
    Matrix<double> out(A.rows(), A.cols() - 1, 0.0);
    for (std::size_t i = 0; i < A.rows(); ++i) {
        std::size_t c = 0;
        for (std::size_t j = 0; j < A.cols(); ++j) {
            if (j == drop) continue;
            out(i, c++) = A(i, j);
        }
    }
    return out;
}

}  // namespace etaqa_detail

// ---------------------------------------------------------------------------
// MG1_G_ETAQA.m
// ---------------------------------------------------------------------------

/**
 * G of an M/G/1-type chain, uniformized first. Port of `MG1_G_ETAQA.m`.
 *
 * `A` is the wide `[A0 A1 ... Amax]`. The generator is turned into the
 * transition matrix of the uniformized chain by dividing through by
 * -min(diag(A1)) and adding the identity back onto A1, which is what cyclic
 * reduction expects; see defect 1 in the header for why that branch is
 * unconditional.
 */
inline Matrix<double> mg1_g_etaqa(const Matrix<double>& A) {
    const std::size_t r = A.rows();
    if (A.cols() % r != 0) throw InputError("MG1_G_ETAQA: A is not a block sequence of A's width");
    Matrix<double> An = A;
    // The reference's isdiscrete flag is write-only (`isdicrete`), so this
    // uniformization always runs.
    double t = An(0, r + 0);
    for (std::size_t i = 1; i < r; ++i) t = std::min(t, An(i, r + i));
    if (t > 0.0)
        throw InputError(
            "MG1_G_ETAQA: this is not a stochastic matrix, neither continuous nor discrete; "
            "every row must sum to 0 or 1");
    for (std::size_t i = 0; i < An.rows(); ++i)
        for (std::size_t j = 0; j < An.cols(); ++j) An(i, j) /= (-t);
    for (std::size_t i = 0; i < r; ++i) An(i, r + i) += 1.0;

    return mg1_cr(blocks_of(An, r));
}

// ---------------------------------------------------------------------------
// MG1_pi_ETAQA.m
// ---------------------------------------------------------------------------

/**
 * Aggregated stationary vector [pi0, pi1, pi2+pi3+...] of an M/G/1-type chain.
 * Port of `MG1_pi_ETAQA.m`.
 *
 * `B` may be empty, in which case the boundary repeats the repetitive blocks.
 * `C0` is the reference's 'Boundary' option, the block that takes level 1 back
 * to a boundary of a different size; pass an empty matrix for the default A0.
 */
inline std::vector<double> mg1_pi_etaqa(const Matrix<double>& Bin, const Matrix<double>& Ain,
                                        const Matrix<double>& G,
                                        const Matrix<double>& C0in = Matrix<double>()) {
    using namespace etaqa_detail;
    Matrix<double> A = Ain;
    const std::size_t m = A.rows();
    const std::size_t dega = A.cols() / m - 1;
    const Matrix<double> A0 = cols_of(A, 0, m);

    Matrix<double> B;
    std::size_t mb = 0, degb = 0;
    if (Bin.empty()) {
        mb = m;
        degb = dega;
        B = A;
    } else {
        B = Bin;
        mb = B.rows();
        if ((B.cols() - mb) % m != 0)
            throw InputError("MG1_pi_ETAQA: matrix B has an incorrect number of columns");
        degb = (B.cols() - mb) / m;
    }
    Matrix<double> C0 = C0in.empty() ? A0 : C0in;
    if (C0in.empty() && mb != m)
        throw InputError(
            "MG1_pi_ETAQA: the Boundary option must be used since a dimension of B0 is not "
            "identical to A0");
    if (!C0in.empty() && (C0in.rows() != m || C0in.cols() != mb))
        throw InputError("MG1_pi_ETAQA: the boundary parameter value has an incorrect dimension");

    // A transition matrix is turned into a generator, as the reference tests it.
    const std::vector<double> brs = rowsums(B);
    double tot = 0.0;
    for (double x : brs) tot += x;
    if (tot > 1e-12 && (tot - static_cast<double>(mb)) < 1e-12) {
        for (std::size_t i = 0; i < mb; ++i) B(i, i) -= 1.0;
        for (std::size_t i = 0; i < m; ++i) A(i, m + i) -= 1.0;
    }

    const Blocks Ab = blocks_of(A, m);
    const Drift d = mg1_drift(Ab);
    if (d.value >= 1.0)
        throw NumericError(
            "MG1_pi_ETAQA: the Markov chain characterized by A is not positive recurrent (drift = " +
            std::to_string(d.value) + ")");

    // Shat(j) = B(j) + B(j+1) G + B(j+2) G^2 + ..., j = 1..degb.
    std::vector<Matrix<double>> Shat;
    Shat.push_back(cols_of(B, mb + (degb - 1) * m, mb + degb * m));
    for (std::size_t i = degb; i-- > 1;) {
        const Matrix<double> temp =
            madd(cols_of(B, mb + (i - 1) * m, mb + i * m), matmul(Shat.front(), G));
        Shat.insert(Shat.begin(), temp);
    }

    // S(j) = A(j) + A(j+1) G + A(j+2) G^2 + ..., j = 1..dega.
    if (dega <= 1)
        throw InputError(
            "MG1_pi_ETAQA: the number of repetitive state blocks is less than 2, this is not an "
            "irreducible Markov chain");
    std::vector<Matrix<double>> S;
    S.push_back(Ab[dega]);
    for (std::size_t i = dega; i-- > 1;) {
        const Matrix<double> temp = madd(Ab[i], matmul(S.front(), G));
        S.insert(S.begin(), temp);
    }

    const Matrix<double> zmb(m, mb, 0.0), zmm(m, m, 0.0);
    Matrix<double> Firstc(mb + 2 * m, 1, 1.0);
    std::vector<Matrix<double>> secondp;
    secondp.push_back(cols_of(B, 0, mb));
    secondp.push_back(C0);
    secondp.push_back(zmb);
    const Matrix<double> Secondc = vcat2(secondp);

    if (Shat.size() < 2) Shat.push_back(Matrix<double>(mb, m, 0.0));
    if (S.size() < 2) S.push_back(Matrix<double>(m, m, 0.0));

    std::vector<Matrix<double>> thirdp;
    thirdp.push_back(madd(cols_of(B, mb, mb + m), matmul(Shat[1], G)));
    thirdp.push_back(madd(Ab[1], matmul(S[1], G)));
    thirdp.push_back(zmm);
    const Matrix<double> Thirdc = vcat2(thirdp);

    Matrix<double> Bsum(mb, m, 0.0), Shat_sum(mb, m, 0.0);
    if (degb <= 2) {
        if (degb == 2) Bsum = cols_of(B, mb + m, mb + 2 * m);
    } else {
        for (std::size_t i = 2; i <= degb - 1; ++i) {
            Bsum = madd(Bsum, cols_of(B, mb + (i - 1) * m, mb + i * m));
            Shat_sum = madd(Shat_sum, Shat[i]);
        }
        Bsum = madd(Bsum, cols_of(B, mb + (degb - 1) * m, B.cols()));
    }

    Matrix<double> Asum(m, m, 0.0), Ssum(m, m, 0.0);
    if (dega >= 3) {
        for (std::size_t i = 2; i <= dega - 1; ++i) {
            Ssum = madd(Ssum, S[i]);
            Asum = madd(Asum, Ab[i]);
        }
        Asum = madd(Asum, Ab[dega]);
    } else if (dega == 2) {
        Asum = madd(Asum, Ab[dega]);
    } else {
        throw InputError(
            "MG1_pi_ETAQA: the number of repetitive state blocks is less than 3, the Markov chain "
            "is reducible");
    }

    std::vector<Matrix<double>> fourthp;
    fourthp.push_back(madd(Bsum, matmul(Shat_sum, G)));
    fourthp.push_back(madd(Asum, matmul(Ssum, G)));
    fourthp.push_back(madd(madd(Asum, Ab[1]), matmul(madd(Ssum, S[1]), G)));
    const Matrix<double> Fourthc = vcat2(fourthp);

    std::vector<Matrix<double>> xparts;
    xparts.push_back(Secondc);
    xparts.push_back(Thirdc);
    xparts.push_back(Fourthc);
    Matrix<double> Xtemp = hcat2(xparts);

    // Drop the first column whose removal leaves the rank unchanged.
    const std::size_t full = matrix_rank(Xtemp);
    const std::size_t n = mb + 2 * m;
    std::size_t drop = n - 1;
    for (std::size_t i = 0; i < n; ++i) {
        if (matrix_rank(drop_col(Xtemp, i)) == full) {
            drop = i;
            break;
        }
    }
    Xtemp = drop_col(Xtemp, drop);

    std::vector<Matrix<double>> xnewp;
    xnewp.push_back(Firstc);
    xnewp.push_back(Xtemp);
    const Matrix<double> Xnew = hcat2(xnewp);

    std::vector<double> rside(n, 0.0);
    rside[0] = 1.0;
    return right_divide(rside, Xnew);
}

// ---------------------------------------------------------------------------
// MG1_qlen_ETAQA.m
// ---------------------------------------------------------------------------

/**
 * n-th moment of the level (the queue length) of an M/G/1-type chain from the
 * ETAQA aggregates. Port of `MG1_qlen_ETAQA.m`.
 */
inline double mg1_qlen_etaqa(const Matrix<double>& Bin, const Matrix<double>& Ain,
                             const std::vector<double>& pi, std::size_t n,
                             const Matrix<double>& C0in = Matrix<double>()) {
    using namespace etaqa_detail;
    const Matrix<double>& A = Ain;
    const std::size_t m = A.rows();
    const std::size_t dega = A.cols() / m - 1;
    const Blocks Ab = blocks_of(A, m);

    Matrix<double> B;
    std::size_t mb = 0, degb = 0;
    if (Bin.empty()) {
        mb = m;
        degb = dega;
        B = A;
    } else {
        B = Bin;
        mb = B.rows();
        if ((B.cols() - mb) % m != 0)
            throw InputError("MG1_qlen_ETAQA: matrix B has an incorrect number of columns");
        degb = (B.cols() - mb) / m;
    }
    if (C0in.empty() && mb != m)
        throw InputError(
            "MG1_qlen_ETAQA: the Boundary option must be used since the column size of B0 is not "
            "identical to A0");
    if (!C0in.empty() && (C0in.rows() != m || C0in.cols() != mb))
        throw InputError("MG1_qlen_ETAQA: the boundary parameter value has an incorrect dimension");

    double mass = 0.0;
    for (double x : pi) mass += x;
    if (std::fabs(mass - 1.0) > 1e-10)
        throw InputError("MG1_qlen_ETAQA: the input probability vector does not sum up to 1");
    if ((pi.size() - mb) % m != 0)
        throw InputError(
            "MG1_qlen_ETAQA: the probability vector has an incorrect number of columns");

    const std::vector<double> pi0(pi.begin(), pi.begin() + mb);
    const std::vector<double> pi1(pi.begin() + mb, pi.begin() + mb + m);
    const std::vector<double> pistar(pi.begin() + mb + m, pi.begin() + mb + 2 * m);

    // lsleft = [(A0+...+Amax) without its last column, (F11 - A0) e].
    Matrix<double> Asum = Ab[0];
    for (std::size_t i = 1; i <= dega; ++i) Asum = madd(Asum, Ab[i]);
    Matrix<double> F11 = Ab[2];
    for (std::size_t i = 3; i <= dega; ++i)
        F11 = madd(F11, mscale(Ab[i], static_cast<double>(i - 1)));
    Matrix<double> lsleft(m, m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j + 1 < m; ++j) lsleft(i, j) = Asum(i, j);
        double s = 0.0;
        for (std::size_t j = 0; j < m; ++j) s += F11(i, j) - Ab[0](i, j);
        lsleft(i, m - 1) = s;
    }

    // Fhat0(j) = sum_{l>=j} B(l), j = 1..degb.
    std::vector<Matrix<double>> Fhat0j;
    Fhat0j.push_back(cols_of(B, mb + (degb - 1) * m, B.cols()));
    for (std::size_t j = degb; j-- > 1;)
        Fhat0j.insert(Fhat0j.begin(),
                      madd(cols_of(B, mb + (j - 1) * m, mb + j * m), Fhat0j.front()));

    // F0(j) = sum_{l>=j+1} A(l): built from block dega down to block 2 but
    // indexed from 1, which is the off-by-one of defect 2 in the header.
    std::vector<Matrix<double>> F0j;
    F0j.push_back(Ab[dega]);
    for (std::size_t j = dega; j-- > 2;) F0j.insert(F0j.begin(), madd(Ab[j], F0j.front()));

    std::vector<std::vector<double>> r;
    r.push_back(pistar);

    // frestsaver(l) = sum_{j>=2} j^l A(j); fcrestsaver(l) = sum_j j^l F0j(j).
    std::vector<Matrix<double>> frestsaver, fcrestsaver;
    for (std::size_t l = 1; l <= n; ++l) {
        Matrix<double> t1(m, m, 0.0);
        for (std::size_t j = 2; j <= dega; ++j)
            t1 = madd(t1, mscale(Ab[j], std::pow(static_cast<double>(j), static_cast<double>(l))));
        frestsaver.push_back(t1);
        Matrix<double> t2(m, m, 0.0);
        for (std::size_t j = 1; j <= dega; ++j)
            if (j <= F0j.size())
                t2 = madd(t2, mscale(F0j[j - 1],
                                     std::pow(static_cast<double>(j), static_cast<double>(l))));
        fcrestsaver.push_back(t2);
    }

    for (std::size_t k = 1; k <= n; ++k) {
        const double dk = static_cast<double>(k);
        Matrix<double> fhatkM(mb, m, 0.0);
        for (std::size_t j = 1; j <= degb; ++j)
            fhatkM = madd(fhatkM, mscale(cols_of(B, mb + (j - 1) * m, mb + j * m),
                                         std::pow(static_cast<double>(j + 1), dk)));
        const std::vector<double> fhatk = vmul(pi0, fhatkM);

        Matrix<double> fkM(m, m, 0.0);
        for (std::size_t j = 2; j <= dega; ++j)
            fkM = madd(fkM, mscale(Ab[j], std::pow(static_cast<double>(j + 1), dk)));
        fkM = madd(mscale(Ab[1], std::pow(2.0, dk)), fkM);
        const std::vector<double> fk = vmul(pi1, fkM);

        std::vector<double> frest(m, 0.0);
        for (std::size_t l = 1; l <= k; ++l) {
            const Matrix<double> M = madd(Ab[1], frestsaver[l - 1]);
            const std::vector<double> t = vmul(r[k - l], M);
            for (std::size_t j = 0; j < m; ++j) frest[j] += bino(k, l) * t[j];
        }

        std::vector<double> bk(m, 0.0);
        for (std::size_t j = 0; j < m; ++j) bk[j] = -fhatk[j] - fk[j] - frest[j];

        Matrix<double> fchatkM(mb, m, 0.0);
        for (std::size_t j = 2; j <= degb; ++j)
            fchatkM = madd(fchatkM, mscale(Fhat0j[j - 1], std::pow(static_cast<double>(j), dk)));
        const double fchatk = vmul_e(pi0, fchatkM);

        Matrix<double> fckM(m, m, 0.0);
        for (std::size_t j = 1; j + 1 <= dega; ++j)
            if (j <= F0j.size())
                fckM = madd(fckM,
                            mscale(F0j[j - 1], std::pow(static_cast<double>(j + 1), dk)));
        const double fck = vmul_e(pi1, fckM);

        double fcrest = 0.0;
        for (std::size_t l = 1; l <= k; ++l)
            fcrest += bino(k, l) * vmul_e(r[k - l], fcrestsaver[l - 1]);

        const double ck = -fchatk - fck - fcrest;

        std::vector<double> rside(m, 0.0);
        for (std::size_t j = 0; j + 1 < m; ++j) rside[j] = bk[j];
        rside[m - 1] = ck;
        r.push_back(right_divide(rside, lsleft));
    }

    return vsum(r.back()) + vsum(pi1);
}

// ---------------------------------------------------------------------------
// GIM1_R_ETAQA.m
// ---------------------------------------------------------------------------

/**
 * R of a GI/M/1-type chain, uniformized first. Port of `GIM1_R_ETAQA.m`.
 *
 * `A` is the VERTICAL stack `[A0; A1; ...; Amax]`, which is how a GI/M/1-type
 * sequence is written; the reference transposes that stack into the horizontal
 * one `GIM1_R` wants, and asks for the automatic dual with functional
 * iterations.
 */
inline Matrix<double> gim1_r_etaqa(const Matrix<double>& A) {
    const std::size_t s = A.cols();
    if (A.rows() % s != 0)
        throw InputError("GIM1_R_ETAQA: A is not a vertical stack of square blocks");
    Matrix<double> An = A;
    double t = An(s + 0, 0);
    for (std::size_t i = 1; i < s; ++i) t = std::min(t, An(s + i, i));
    const bool isdiscrete = t > 0.0;
    if (!isdiscrete) {
        if (t > 0.0)
            throw InputError(
                "GIM1_R_ETAQA: this is not a stochastic matrix, neither continuous nor discrete");
        for (std::size_t i = 0; i < An.rows(); ++i)
            for (std::size_t j = 0; j < An.cols(); ++j) An(i, j) /= (-t);
        for (std::size_t i = 0; i < s; ++i) An(s + i, i) += 1.0;
    }
    return gim1_r(vblocks_of(An, s), "A", "FI");
}

// ---------------------------------------------------------------------------
// GIM1_pi_ETAQA.m
// ---------------------------------------------------------------------------

/**
 * Aggregated stationary vector [pi0, pi1, pi2+pi3+...] of a GI/M/1-type chain.
 * Port of `GIM1_pi_ETAQA.m`. `B` and `A` are vertical stacks; `B0` is the
 * reference's 'Boundary' option (pass empty for the default A0).
 */
inline std::vector<double> gim1_pi_etaqa(const Matrix<double>& Bin, const Matrix<double>& Ain,
                                         const Matrix<double>& R,
                                         const Matrix<double>& B0in = Matrix<double>()) {
    using namespace etaqa_detail;
    Matrix<double> B = Bin, A = Ain;
    const std::size_t m = R.rows();
    const std::size_t mb = B.cols();
    if ((B.rows() - mb) % m != 0)
        throw InputError("GIM1_pi_ETAQA: input matrix B has an incorrect number of rows");
    const std::size_t degb = (B.rows() - mb) / m;
    if (A.rows() % m != 0)
        throw InputError("GIM1_pi_ETAQA: input matrix A has an incorrect number of rows");
    const std::size_t dega = A.rows() / m - 1;

    const Matrix<double> ImR = msub(eye<double>(m), R);
    const Matrix<double> temp0 = inverse(ImR);
    bool all_cols_negative = true;
    for (std::size_t j = 0; j < m && all_cols_negative; ++j) {
        bool any = false;
        for (std::size_t i = 0; i < m; ++i)
            if (temp0(i, j) < -100.0 * 2.220446049250313e-16) any = true;
        all_cols_negative = any;
    }
    if (all_cols_negative)
        throw NumericError(
            "GIM1_pi_ETAQA: the spectral radius of R is not below 1, GIM1 is not positive "
            "recurrent");

    Matrix<double> B0;
    if (!B0in.empty()) {
        if (B0in.rows() != mb)
            throw InputError("GIM1_pi_ETAQA: Boundary has an incorrect number of rows");
        if (B0in.cols() != m)
            throw InputError("GIM1_pi_ETAQA: Boundary has an incorrect number of columns");
        B0 = B0in;
    } else {
        B0 = Matrix<double>(m, A.cols(), 0.0);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < A.cols(); ++j) B0(i, j) = A(i, j);
    }

    // A transition matrix is turned into a generator, as the reference tests it.
    Matrix<double> Btop(mb, B.cols(), 0.0);
    for (std::size_t i = 0; i < mb; ++i)
        for (std::size_t j = 0; j < B.cols(); ++j) Btop(i, j) = B(i, j);
    const Matrix<double> test = madd(Btop, B0);
    const std::vector<double> trs = rowsums(test);
    double tot = 0.0;
    for (double x : trs) tot += x - 1.0;
    if (std::fabs(tot) < 1e-10) {
        for (std::size_t i = 0; i < mb; ++i) B(i, i) -= 1.0;
        for (std::size_t i = 0; i < m; ++i) A(m + i, i) -= 1.0;
    }

    const Blocks Ab = vblocks_of(A, m);

    const Matrix<double> Firstc(mb + 2 * m, 1, 1.0);

    // sum_{i>=2} R^(i-2) (I - R) B(i)
    Matrix<double> temp = msub(eye<double>(m), R);
    Matrix<double> tempsum(m, mb, 0.0);
    for (std::size_t i = 2; i <= degb; ++i) {
        tempsum = madd(tempsum, matmul(temp, vblock_after(B, mb, m, i)));
        temp = matmul(R, temp);
    }
    std::vector<Matrix<double>> secondp;
    secondp.push_back(Btop);
    secondp.push_back(vblock_after(B, mb, m, 1));
    secondp.push_back(tempsum);
    const Matrix<double> Secondc = vcat2(secondp);

    // sum_{i>=2} R^(i-2) (I - R) A(i)
    temp = msub(eye<double>(m), R);
    tempsum = Matrix<double>(m, m, 0.0);
    for (std::size_t i = 2; i <= dega; ++i) {
        tempsum = madd(tempsum, matmul(temp, Ab[i]));
        temp = matmul(R, temp);
    }
    std::vector<Matrix<double>> thirdp;
    thirdp.push_back(B0);
    thirdp.push_back(Ab[1]);
    thirdp.push_back(tempsum);
    const Matrix<double> Thirdc = vcat2(thirdp);

    // sum_{i>=2} R^(i-1) A(i)
    temp = R;
    tempsum = Matrix<double>(m, m, 0.0);
    for (std::size_t i = 2; i <= dega; ++i) {
        tempsum = madd(tempsum, matmul(temp, Ab[i]));
        temp = matmul(R, temp);
    }
    std::vector<Matrix<double>> fourthp;
    fourthp.push_back(Matrix<double>(mb, m, 0.0));
    fourthp.push_back(Ab[0]);
    fourthp.push_back(madd(madd(Ab[0], Ab[1]), tempsum));
    Matrix<double> Fourthc = vcat2(fourthp);
    Fourthc = cols_of(Fourthc, 0, Fourthc.cols() - 1);

    std::vector<Matrix<double>> xparts;
    xparts.push_back(Firstc);
    xparts.push_back(Secondc);
    xparts.push_back(Thirdc);
    xparts.push_back(Fourthc);
    const Matrix<double> X = hcat2(xparts);

    std::vector<double> rside(mb + 2 * m, 0.0);
    rside[0] = 1.0;
    return right_divide(rside, X);
}

// ---------------------------------------------------------------------------
// GIM1_qlen_ETAQA.m
// ---------------------------------------------------------------------------

/**
 * n-th moment of the level of a GI/M/1-type chain from the ETAQA aggregates.
 * Port of `GIM1_qlen_ETAQA.m`, including the scalar `A(3)` of defect 3.
 */
inline double gim1_qlen_etaqa(const Matrix<double>& Bin, const Matrix<double>& Ain,
                              const Matrix<double>& R, const std::vector<double>& pi,
                              std::size_t n, const Matrix<double>& B0in = Matrix<double>()) {
    using namespace etaqa_detail;
    Matrix<double> B = Bin, A = Ain;
    const std::size_t m = R.rows();
    const std::size_t mb = B.cols();
    if ((B.rows() - mb) % m != 0)
        throw InputError("GIM1_qlen_ETAQA: input matrix B has an incorrect number of rows");
    const std::size_t degb = (B.rows() - mb) / m;
    if (A.rows() % m != 0)
        throw InputError("GIM1_qlen_ETAQA: input matrix A has an incorrect number of rows");
    const std::size_t dega = A.rows() / m - 1;

    const Matrix<double> temp0 = inverse(msub(eye<double>(m), R));
    bool all_cols_negative = true;
    for (std::size_t j = 0; j < m && all_cols_negative; ++j) {
        bool any = false;
        for (std::size_t i = 0; i < m; ++i)
            if (temp0(i, j) < -100.0 * 2.220446049250313e-16) any = true;
        all_cols_negative = any;
    }
    if (all_cols_negative)
        throw NumericError(
            "GIM1_qlen_ETAQA: the spectral radius of R is not below 1, GIM1 is not positive "
            "recurrent");

    Matrix<double> B0;
    if (!B0in.empty()) {
        if (B0in.rows() != mb || B0in.cols() != m)
            throw InputError("GIM1_qlen_ETAQA: Boundary has an incorrect dimension");
        B0 = B0in;
    } else {
        B0 = Matrix<double>(m, A.cols(), 0.0);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < A.cols(); ++j) B0(i, j) = A(i, j);
    }

    const std::vector<double> pi0(pi.begin(), pi.begin() + mb);
    const std::vector<double> pi1(pi.begin() + mb, pi.begin() + mb + m);
    const std::vector<double> pistar(pi.begin() + mb + m, pi.begin() + mb + 2 * m);

    if (n == 0) return 1.0;

    // The scalar of defect 3, read at column-major linear index 3 of A.
    const double a3 = A(2 % A.rows(), 2 / A.rows());

    Matrix<double> Btop(mb, B.cols(), 0.0);
    for (std::size_t i = 0; i < mb; ++i)
        for (std::size_t j = 0; j < B.cols(); ++j) Btop(i, j) = B(i, j);
    const Matrix<double> test = madd(Btop, B0);
    const std::vector<double> trs = rowsums(test);
    double tot = 0.0;
    for (double x : trs) tot += x - 1.0;
    if (std::fabs(tot) < 1e-10) {
        for (std::size_t i = 0; i < mb; ++i) B(i, i) -= 1.0;
        for (std::size_t i = 0; i < m; ++i) A(m + i, i) -= 1.0;
    }

    const Blocks Ab = vblocks_of(A, m);

    Matrix<double> Rpower = eye<double>(m);
    Matrix<double> lsum = madd(Ab[0], Ab[1]);
    for (std::size_t i = 2; i <= dega; ++i) {
        lsum = madd(lsum, matmul(Rpower, Ab[i]));
        Rpower = matmul(R, Rpower);
    }

    std::vector<double> leftr(m, 0.0);
    if (degb >= 2 && dega >= 2) {
        // leftr_part1 starts as the SCALAR a3 and is broadcast over the m x m
        // accumulator only if the loop below runs at all.
        const bool loop_runs = dega >= 3;
        Matrix<double> part1(m, m, a3), part2(m, m, 0.0);
        Rpower = R;
        for (std::size_t i = 1; i + 2 <= dega; ++i) {
            part1 = madd(part1, matmul(Rpower, Ab[i + 2]));
            part2 = madd(part2, mscale(matmul(Rpower, Ab[i + 2]), static_cast<double>(i)));
            Rpower = matmul(R, Rpower);
        }
        for (std::size_t i = 0; i < m; ++i) {
            double s = 0.0;
            if (loop_runs) {
                for (std::size_t j = 0; j < m; ++j) s += part1(i, j) + part2(i, j);
            } else {
                s = a3;  // scalar * ones(m,1), never expanded to a matrix
            }
            for (std::size_t j = 0; j < m; ++j) s -= Ab[0](i, j);
            leftr[i] = s;
        }
    } else if (degb == 1 && dega != 1) {
        Matrix<double> acc(m, m, 0.0);
        Rpower = eye<double>(m);
        for (std::size_t i = 2; i <= dega; ++i) {
            acc = madd(acc, mscale(matmul(Rpower, Ab[i]), static_cast<double>(i - 1)));
            Rpower = matmul(R, Rpower);
        }
        for (std::size_t i = 0; i < m; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < m; ++j) s += acc(i, j) - Ab[0](i, j);
            leftr[i] = s;
        }
    } else {
        throw InputError(
            "GIM1_qlen_ETAQA: the number of A blocks is not enough, this is a reducible Markov "
            "chain");
    }

    Matrix<double> lsleft(m, m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j + 1 < m; ++j) lsleft(i, j) = lsum(i, j);
        lsleft(i, m - 1) = leftr[i];
    }

    std::vector<std::vector<double>> r;
    r.push_back(pistar);

    std::vector<std::vector<double>> reuse;  // rsidepart2reuse, one column per k

    for (std::size_t k = 1; k <= n; ++k) {
        const double dk = static_cast<double>(k);
        std::vector<double> bk(m, 0.0);
        {
            const std::vector<double> t0 = vmul(pi0, B0);
            const Matrix<double> M =
                madd(mscale(Ab[1], std::pow(2.0, dk)), mscale(Ab[0], std::pow(3.0, dk)));
            const std::vector<double> t1 = vmul(pi1, M);
            for (std::size_t j = 0; j < m; ++j) bk[j] = -(std::pow(2.0, dk) * t0[j] + t1[j]);
            for (std::size_t l = 1; l <= k; ++l) {
                const Matrix<double> Ml =
                    madd(Ab[1], mscale(Ab[0], std::pow(2.0, static_cast<double>(l))));
                const std::vector<double> t = vmul(r[k - l], Ml);
                for (std::size_t j = 0; j < m; ++j) bk[j] -= bino(k, l) * t[j];
            }
        }

        Matrix<double> tempsum(m, m, 0.0);
        Rpower = R;
        for (std::size_t i = 1; i + 2 <= dega; ++i) {
            double t = 0.0;
            for (std::size_t z = 1; z <= i; ++z) t += std::pow(static_cast<double>(z), dk);
            tempsum = madd(tempsum, mscale(matmul(Rpower, Ab[i + 2]), t));
            Rpower = matmul(Rpower, R);
        }
        std::vector<double> col(m, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < m; ++j) s += Ab[0](i, j) - tempsum(i, j);
            col[i] = s;
        }
        reuse.push_back(col);

        double ck = std::pow(2.0, dk) * vmul_e(pi1, Ab[0]);
        Matrix<double> tempsum2(m, mb, 0.0);
        Rpower = R;
        for (std::size_t i = 2; i <= degb; ++i) {
            double t = 0.0;
            for (std::size_t z = 2; z <= i; ++z) t += std::pow(static_cast<double>(z), dk);
            tempsum2 = madd(tempsum2, mscale(matmul(Rpower, vblock_after(B, mb, m, i)), t));
            Rpower = matmul(R, Rpower);
        }
        ck -= vmul_e(pi1, tempsum2);
        for (std::size_t l = 1; l <= k; ++l) {
            double s = 0.0;
            for (std::size_t j = 0; j < m; ++j) s += r[k - l][j] * reuse[l - 1][j];
            ck += bino(k, l) * s;
        }

        std::vector<double> rside(m, 0.0);
        for (std::size_t j = 0; j + 1 < m; ++j) rside[j] = bk[j];
        rside[m - 1] = ck;
        r.push_back(right_divide(rside, lsleft));
    }

    return vsum(r.back()) + vsum(pi1);
}

}  // namespace smc
}  // namespace line

#endif  // LINE_LIB_SMC_ETAQA_H
