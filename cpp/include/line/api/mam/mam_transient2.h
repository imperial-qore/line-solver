/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAM_TRANSIENT2_H
#define LINE_API_MAM_MAM_TRANSIENT2_H

/**
 * Laplace-domain transient level-to-level transform V(s,n,m) of a piecewise
 * level-dependent QBD, the port of `matlab/src/solvers/MAM/mam_transient2.m`
 * (finite, closed at a top level) and `mam_transient2_open.m` (the last regime
 * repeating to infinity).
 *
 * V(s,n,m) is the Laplace transform, at complex argument s, of the transient
 * transition-probability matrix from level n to level m. The construction
 * follows the Horvath et al. formulation the reference is ported from: per
 * regime, the fundamental matrices G and R of the shifted blocks (L - sI) are
 * computed in both directions (up and down), the boundary levels are linked by
 * the "SY / SYh" backward and forward recursions, the threshold-to-threshold
 * transforms SV are assembled, and a level strictly inside a regime is reached
 * by interpolating between its two enclosing thresholds with the R matrices.
 *
 * INDEXING. The block cell arrays are indexed by regime and are kept 1-based
 * here, exactly as in the reference, because the recursions mix k, k+1 and k-1
 * with the threshold vector T and every off-by-one is silent -- it returns a
 * transform of the wrong level rather than failing. The accessors `B1`, `L1`,
 * `F1`, `Lv1` and `T1` take the reference's own index.
 *
 * ARITHMETIC. Complex double only. The blocks are shifted by -sI at a complex
 * quadrature node, so G, R and every intermediate are complex; see
 * `line/num/complex_number.h` for why the generic linear algebra is available
 * at that element type.
 */

#include <algorithm>
#include <complex>
#include <cstddef>
#include <vector>

#include "line/num/complex_number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

using CMat = Matrix<Complex>;

namespace transient_detail {

/** A \ B: LU of A once, then one back substitution per column of B. */
inline CMat mldivide(const CMat& A, const CMat& B) {
    if (A.rows() != A.cols()) throw InputError("mam_transient2: left division needs a square matrix");
    if (A.rows() != B.rows()) throw InputError("mam_transient2: left division dimension mismatch");
    CMat LU = A;
    const std::vector<std::size_t> piv = lu_factor(LU);
    CMat X(B.rows(), B.cols());
    std::vector<Complex> col(B.rows());
    for (std::size_t j = 0; j < B.cols(); ++j) {
        for (std::size_t i = 0; i < B.rows(); ++i) col[i] = B(i, j);
        lu_solve(LU, piv, col);
        for (std::size_t i = 0; i < B.rows(); ++i) X(i, j) = col[i];
    }
    return X;
}

/** A / B, as the transposed left division so it is one factorization too. */
inline CMat mrdivide(const CMat& A, const CMat& B) {
    return mldivide(B.transpose(), A.transpose()).transpose();
}

/** MATLAB's mxpow: A^k, and the identity when k is zero. */
inline CMat mxpow(const CMat& A, long k) {
    if (k < 0) throw InputError("mam_transient2: mxpow requires a non-negative exponent");
    if (k == 0) return eye<Complex>(A.rows());
    return matpow(A, static_cast<unsigned>(k));
}

inline CMat madd(const CMat& A, const CMat& B) {
    if (A.rows() != B.rows() || A.cols() != B.cols())
        throw InputError("mam_transient2: addition shape mismatch");
    CMat C(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = A(i, j) + B(i, j);
    return C;
}

inline CMat msub(const CMat& A, const CMat& B) {
    if (A.rows() != B.rows() || A.cols() != B.cols())
        throw InputError("mam_transient2: subtraction shape mismatch");
    CMat C(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) C(i, j) = A(i, j) - B(i, j);
    return C;
}

/** s*I(n) of complex order n. */
inline CMat sI(const Complex& s, std::size_t n) {
    CMat M(n, n, Complex(0.0, 0.0));
    for (std::size_t i = 0; i < n; ++i) M(i, i) = s;
    return M;
}

/** [A, B; C, D] with all four blocks of equal order. */
inline CMat block2(const CMat& A, const CMat& B, const CMat& C, const CMat& D) {
    const std::size_t n = A.rows();
    CMat M(2 * n, 2 * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            M(i, j) = A(i, j);
            M(i, n + j) = B(i, j);
            M(n + i, j) = C(i, j);
            M(n + i, n + j) = D(i, j);
        }
    return M;
}

/** [A; B] stacked vertically. */
inline CMat vcat(const CMat& A, const CMat& B) {
    if (A.cols() != B.cols()) throw InputError("mam_transient2: vcat width mismatch");
    CMat M(A.rows() + B.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) M(i, j) = A(i, j);
    for (std::size_t i = 0; i < B.rows(); ++i)
        for (std::size_t j = 0; j < B.cols(); ++j) M(A.rows() + i, j) = B(i, j);
    return M;
}

/** Rows [r0, r0+nr) and columns [c0, c0+nc). */
inline CMat sub(const CMat& A, std::size_t r0, std::size_t nr, std::size_t c0, std::size_t nc) {
    CMat M(nr, nc);
    for (std::size_t i = 0; i < nr; ++i)
        for (std::size_t j = 0; j < nc; ++j) M(i, j) = A(r0 + i, c0 + j);
    return M;
}

/** Maximum absolute row sum, MATLAB's norm(A, inf) at complex argument. */
inline double norminf(const CMat& A) {
    double best = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < A.cols(); ++j) s += std::abs(A(i, j));
        if (s > best) best = s;
    }
    return best;
}

}  // namespace transient_detail

/** G and R of a QBD whose local block is complex. */
struct CQbdFundMat {
    CMat G;
    CMat R;
    unsigned iterations = 0;
};

/**
 * Port of `qbd_fundmat.m` at complex argument: cyclic reduction (Bini-Meini
 * logarithmic reduction) on the raw level blocks, with R recovered from G.
 *
 * The reference uniformizes by `lamb = max(-real(diag(L)))` -- the REAL part,
 * because the shift -sI makes the diagonal complex while the uniformization
 * constant must stay a positive real, and stops on `min(norm(BB,inf),
 * norm(BF,inf))`, a real quantity too. Both are reproduced.
 */
inline CQbdFundMat qbd_fundmat_laplace(const CMat& B, const CMat& L, const CMat& F,
                                       double precision = 1e-14, unsigned maxNumIt = 50) {
    using namespace transient_detail;
    const std::size_t m = L.rows();
    if (L.cols() != m || B.rows() != m || B.cols() != m || F.rows() != m || F.cols() != m)
        throw InputError("qbd_fundmat_laplace: B, L and F must be square and of equal order");
    const CMat II = eye<Complex>(m);

    double lamb = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < m; ++i) lamb = std::max(lamb, -L(i, i).real());
    if (!(lamb > 0.0))
        throw NumericError("qbd_fundmat_laplace: the local block has no negative real diagonal");

    CMat Bm(m, m), Lm(m, m), Fm(m, m);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) {
            Bm(i, j) = B(i, j) / lamb;
            Lm(i, j) = L(i, j) / lamb + (i == j ? Complex(1.0, 0.0) : Complex(0.0, 0.0));
            Fm(i, j) = F(i, j) / lamb;
        }

    CMat BF = mldivide(msub(II, Lm), II);
    CMat BB = matmul(BF, Fm);
    BF = matmul(BF, Bm);
    CMat G = BF;
    CMat PI = BB;
    double check = 1.0;
    unsigned numit = 0;
    while (check > precision && numit < maxNumIt) {
        const CMat Lstar = madd(matmul(BF, BB), matmul(BB, BF));
        const CMat Bstar = matmul(BB, BB);
        const CMat Fstar = matmul(BF, BF);
        // The reference reuses BB for the inverse before overwriting it.
        BB = mldivide(msub(II, Lstar), II);
        BF = matmul(BB, Fstar);
        BB = matmul(BB, Bstar);
        G = madd(G, matmul(PI, BF));
        PI = matmul(PI, BB);
        check = std::min(norminf(BB), norminf(BF));
        ++numit;
    }

    CQbdFundMat out;
    out.G = G;
    out.R = matmul(Fm, mldivide(msub(II, madd(Lm, matmul(Fm, G))), II));
    out.iterations = numit;
    return out;
}

/**
 * The piecewise QBD blocks, 1-based exactly as the reference's cell arrays.
 *
 * `T` holds the regime thresholds. For the OPEN form there are K = |T| regimes
 * and the last repeats to infinity; for the FINITE form there are K = |T| - 1
 * regimes and `Lv[K+1]` is the top boundary level.
 */
struct TransientQbd {
    std::vector<CMat> B, L, F, Lv;  ///< 1-based: entry 0 is unused padding
    std::vector<long> T;            ///< 1-based: entry 0 is unused padding

    const CMat& B1(std::size_t k) const { return B.at(k); }
    const CMat& L1(std::size_t k) const { return L.at(k); }
    const CMat& F1(std::size_t k) const { return F.at(k); }
    const CMat& Lv1(std::size_t k) const { return Lv.at(k); }
    long T1(std::size_t k) const { return T.at(k); }
    std::size_t nT() const { return T.size() - 1; }
};

/** Build the 1-based padded form from plain 0-based vectors. */
inline TransientQbd make_transient_qbd(const std::vector<CMat>& B, const std::vector<CMat>& L,
                                       const std::vector<CMat>& F, const std::vector<CMat>& Lv,
                                       const std::vector<long>& T) {
    TransientQbd q;
    q.B.push_back(CMat());
    q.L.push_back(CMat());
    q.F.push_back(CMat());
    q.Lv.push_back(CMat());
    q.T.push_back(0);
    q.B.insert(q.B.end(), B.begin(), B.end());
    q.L.insert(q.L.end(), L.begin(), L.end());
    q.F.insert(q.F.end(), F.begin(), F.end());
    q.Lv.insert(q.Lv.end(), Lv.begin(), Lv.end());
    q.T.insert(q.T.end(), T.begin(), T.end());
    return q;
}

namespace transient_detail {

/** The regime containing level n: the reference's `find(T>n,1)-1`, floored. */
inline std::size_t regime_of(const TransientQbd& q, long n, std::size_t fallback) {
    for (std::size_t i = 1; i <= q.nT(); ++i)
        if (q.T1(i) > n) return i - 1;
    return fallback;
}

/**
 * The four H blocks of a d-step level range, shared by both forms.
 *
 * `num / den` with num = [Gh^(d-1), G; Gh, G^(d-1)] and den = [I, G^d; Gh^d, I],
 * sliced into the (n->n, n->0, hat n->n, hat n->0) quadrants.
 */
struct HBlocks {
    CMat nn, n0, hnn, hn0;
};

inline HBlocks h_blocks(const CMat& G, const CMat& Gh, long d, std::size_t NN) {
    const CMat II = eye<Complex>(NN);
    const CMat num = block2(mxpow(Gh, d - 1), G, Gh, mxpow(G, d - 1));
    const CMat den = block2(II, mxpow(G, d), mxpow(Gh, d), II);
    const CMat Tmp = mrdivide(num, den);
    HBlocks h;
    h.nn = sub(Tmp, 0, NN, 0, NN);
    h.n0 = sub(Tmp, 0, NN, NN, NN);
    h.hnn = sub(Tmp, NN, NN, 0, NN);
    h.hn0 = sub(Tmp, NN, NN, NN, NN);
    return h;
}

/**
 * The R-interpolation that closes both forms: a level m strictly between the
 * thresholds Ll and Lu is a convex combination, in the transform domain, of the
 * transforms at those two thresholds.
 */
inline CMat interpolate(const CMat& Rk, const CMat& Rhk, const CMat& Vl, const CMat& Vu, long m,
                        long Ll, long Lu) {
    if (Rk.rows() == 0)
        throw NumericError(
            "mam_transient2: the interpolation needs the repeating R of a regime wider than one "
            "level, but that regime has none");
    if (Vl.rows() == 0 || Vu.rows() == 0)
        throw NumericError("mam_transient2: the interpolation needs both enclosing transforms");
    const std::size_t NN = Rk.rows();
    const CMat II = eye<Complex>(NN);
    const CMat Zden = block2(II, mxpow(Rk, Lu - Ll), mxpow(Rhk, Lu - Ll), II);
    const CMat Znum = vcat(mxpow(Rk, m - Ll), mxpow(Rhk, Lu - m));
    const CMat Z = mldivide(Zden, Znum);
    return madd(matmul(Vl, sub(Z, 0, NN, 0, Z.cols())),
                matmul(Vu, sub(Z, NN, NN, 0, Z.cols())));
}

}  // namespace transient_detail

/**
 * V(s,n,m) for an OPEN piecewise QBD, the port of `mam_transient2_open.m`.
 *
 * K = |T| regimes; regime K repeats to infinity, so a target level above the
 * last threshold is reached by post-multiplying with a power of R{K}.
 */
inline CMat mam_transient2_open(const TransientQbd& q, long n, long m, const Complex& s) {
    using namespace transient_detail;
    const std::size_t K = q.nT();
    if (K == 0) throw InputError("mam_transient2_open: at least one regime is required");
    if (n < 0 || m < 0) throw InputError("mam_transient2_open: levels must be non-negative");

    std::vector<CMat> Gs(K + 1), Rs(K + 1), Ghs(K + 1), Rhs(K + 1);
    for (std::size_t k = 1; k <= K; ++k) {
        if (k < K && q.T1(k + 1) - q.T1(k) == 1) continue;  // no repeating level in the regime
        const CMat Lk = msub(q.L1(k), sI(s, q.L1(k).rows()));
        const CQbdFundMat gr = qbd_fundmat_laplace(q.B1(k), Lk, q.F1(k));
        Gs[k] = gr.G;
        Rs[k] = gr.R;
        const CQbdFundMat grh = qbd_fundmat_laplace(q.F1(k), Lk, q.B1(k));
        Ghs[k] = grh.G;
        Rhs[k] = grh.R;
    }

    std::vector<CMat> SvHn(K + 1), SvH0(K + 1), SvHhn(K + 1), SvHh0(K + 1);
    for (std::size_t k = 1; k + 1 <= K; ++k) {
        const std::size_t NN = q.Lv1(k).rows();
        if (q.T1(k + 1) - q.T1(k) > 1) {
            const HBlocks h = h_blocks(Gs[k], Ghs[k], q.T1(k + 1) - q.T1(k), NN);
            SvHn[k] = h.nn;
            SvH0[k] = h.n0;
            SvHhn[k] = h.hnn;
            SvHh0[k] = h.hn0;
        } else {
            const std::size_t NN1 = q.Lv1(k + 1).rows();
            SvH0[k] = CMat(NN1, NN, Complex(0.0, 0.0));
            SvHh0[k] = eye<Complex>(NN);
            SvHhn[k] = CMat(NN, NN1, Complex(0.0, 0.0));
            SvHn[k] = eye<Complex>(NN1);
        }
    }

    std::vector<CMat> SY(K + 2);
    SY[K] = Gs[K];
    for (std::size_t k = K - 1; k >= 1; --k) {
        const std::size_t NNk1 = q.Lv1(k + 1).rows();
        const CMat M = msub(msub(msub(sI(s, NNk1), q.Lv1(k + 1)), matmul(q.F1(k + 1), SY[k + 1])),
                            matmul(q.B1(k), SvHhn[k]));
        SY[k] = madd(SvH0[k], matmul(SvHn[k], mldivide(M, matmul(q.B1(k), SvHh0[k]))));
        if (k == 1) break;
    }

    const std::size_t NN1 = q.Lv1(1).rows();
    std::vector<CMat> SYh(K + 1);
    if (K >= 2) {
        const CMat M = msub(msub(sI(s, NN1), q.Lv1(1)), matmul(q.F1(1), SvH0[1]));
        SYh[1] = madd(SvHhn[1], matmul(SvHh0[1], mldivide(M, matmul(q.F1(1), SvHn[1]))));
    }
    for (std::size_t k = 2; k + 1 <= K; ++k) {
        const std::size_t NNk = q.Lv1(k).rows();
        const CMat M = msub(msub(msub(sI(s, NNk), q.Lv1(k)), matmul(q.B1(k - 1), SYh[k - 1])),
                            matmul(q.F1(k), SvH0[k]));
        SYh[k] = madd(SvHhn[k], matmul(SvHh0[k], mldivide(M, matmul(q.F1(k), SvHn[k]))));
    }

    // SV[k][l] is the threshold-to-threshold transform, 1-based in both indices.
    std::vector<std::vector<CMat>> SV(K + 2, std::vector<CMat>(K + 2));
    for (std::size_t l = 0; l + 1 <= K; ++l) {
        if (l == 0) {
            SV[1][1] = mldivide(msub(msub(sI(s, NN1), q.Lv1(1)), matmul(q.F1(1), SY[1])),
                                eye<Complex>(NN1));
        } else {
            const std::size_t NNl1 = q.Lv1(l + 1).rows();
            const CMat M = msub(msub(msub(sI(s, NNl1), q.Lv1(l + 1)), matmul(q.F1(l + 1), SY[l + 1])),
                                matmul(q.B1(l), SYh[l]));
            SV[l + 1][l + 1] = mldivide(M, eye<Complex>(NNl1));
        }
        for (std::size_t k = l + 1; k + 1 <= K; ++k) {
            const std::size_t NNk1 = q.Lv1(k + 1).rows();
            const CMat M = msub(msub(msub(sI(s, NNk1), q.Lv1(k + 1)), matmul(q.F1(k + 1), SY[k + 1])),
                                matmul(q.B1(k), SvHhn[k]));
            SV[k + 1][l + 1] = mldivide(M, matmul(matmul(q.B1(k), SvHh0[k]), SV[k][l + 1]));
        }
        for (std::size_t k = l; k-- > 1;) {
            const std::size_t NNk1 = q.Lv1(k + 1).rows();
            const CMat M = msub(msub(msub(sI(s, NNk1), q.Lv1(k + 1)), matmul(q.F1(k + 1), SvH0[k + 1])),
                                matmul(q.B1(k), SYh[k]));
            SV[k + 1][l + 1] = mldivide(M, matmul(matmul(q.F1(k + 1), SvHn[k + 1]), SV[k + 2][l + 1]));
        }
        if (l > 0) {
            const CMat M = msub(msub(sI(s, NN1), q.Lv1(1)), matmul(q.F1(1), SvH0[1]));
            SV[1][l + 1] = mldivide(M, matmul(matmul(q.F1(1), SvHn[1]), SV[2][l + 1]));
        }
    }

    const std::size_t kn = regime_of(q, n, K);
    const std::size_t km = regime_of(q, m, K);

    const std::size_t NNkn = q.Lv1(kn).rows();
    const CMat IIkn = eye<Complex>(NNkn);

    CMat Vu, Vl;
    long Lu = 0, Ll = 0;

    if (q.T1(kn) == n) {
        if (q.T1(km) == m) return SV[kn][km];
        Vl = SV[kn][km];
        Ll = q.T1(km);
        if (km < K) {
            Vu = SV[kn][km + 1];
            Lu = q.T1(km + 1);
        }
    } else {
        HBlocks hT, hn;
        CMat Yn;
        if (kn < K) {
            hT = h_blocks(Gs[kn], Ghs[kn], q.T1(kn + 1) - n, NNkn);
            hn = h_blocks(Gs[kn], Ghs[kn], n - q.T1(kn), NNkn);
            const std::size_t NNkn1 = q.Lv1(kn + 1).rows();
            const CMat M = msub(msub(msub(sI(s, NNkn1), q.Lv1(kn + 1)), matmul(q.F1(kn + 1), SY[kn + 1])),
                                matmul(q.B1(kn), hT.hnn));
            Yn = madd(hT.n0, matmul(hT.nn, mldivide(M, matmul(q.B1(kn), hT.hn0))));
        } else {
            hn = h_blocks(Gs[kn], Ghs[kn], n - q.T1(kn), NNkn);
            Yn = Gs[kn];
        }

        CMat Yhn;
        if (kn == 1) {
            const CMat M = msub(msub(sI(s, NN1), q.Lv1(1)), matmul(q.F1(1), hn.n0));
            Yhn = madd(hn.hnn, matmul(hn.hn0, mldivide(M, matmul(q.F1(1), hn.nn))));
        } else {
            const std::size_t NNk = q.Lv1(kn).rows();
            const CMat M = msub(msub(msub(sI(s, NNk), q.Lv1(kn)), matmul(q.B1(kn - 1), SYh[kn - 1])),
                                matmul(q.F1(kn), hn.n0));
            Yhn = madd(hn.hnn, matmul(hn.hn0, mldivide(M, matmul(q.F1(kn), hn.nn))));
        }

        const CMat Mkn = msub(sI(s, q.L1(kn).rows()), q.L1(kn));
        CMat Vnl;
        if (q.T1(km) < n) {
            const CMat M = msub(msub(Mkn, matmul(q.B1(kn), hn.hnn)), matmul(q.F1(kn), Yn));
            Vnl = mldivide(M, matmul(matmul(q.B1(kn), hn.hn0), SV[kn][km]));
        } else {
            const CMat M = msub(msub(Mkn, matmul(q.F1(kn), hT.n0)), matmul(q.B1(kn), Yhn));
            Vnl = mldivide(M, matmul(matmul(q.F1(kn), hT.nn), SV[kn + 1][km]));
        }
        if (m == q.T1(km)) return Vnl;

        CMat Vnu;
        if (km < K) {
            if (q.T1(km + 1) < n) {
                const CMat M = msub(msub(Mkn, matmul(q.B1(kn), hn.hnn)), matmul(q.F1(kn), Yn));
                Vnu = mldivide(M, matmul(matmul(q.B1(kn), hn.hn0), SV[kn][km + 1]));
            } else {
                const CMat M = msub(msub(Mkn, matmul(q.F1(kn), hT.n0)), matmul(q.B1(kn), Yhn));
                Vnu = mldivide(M, matmul(matmul(q.F1(kn), hT.nn), SV[kn + 1][km + 1]));
            }
        }
        const CMat Vnn = mldivide(
            msub(msub(msub(sI(s, NNkn), q.L1(kn)), matmul(q.B1(kn), Yhn)), matmul(q.F1(kn), Yn)),
            IIkn);
        if (n == m) return Vnn;
        if (km == K && n < m) {
            if (kn < K) return matmul(Vnl, mxpow(Rs[km], m - q.T1(km)));
            return matmul(Vnn, mxpow(Rs[km], m - n));
        }
        if (kn != km) {
            Vu = Vnu; Vl = Vnl; Lu = q.T1(km + 1); Ll = q.T1(km);
        } else if (n <= m) {
            Vu = Vnu; Vl = Vnn; Lu = q.T1(km + 1); Ll = n;
        } else {
            Vu = Vnn; Vl = Vnl; Lu = n; Ll = q.T1(km);
        }
    }

    if (km == K && n < m) {
        if (kn < K) return matmul(Vl, mxpow(Rs[km], m - q.T1(km)));
        return matmul(SV[kn][kn], mxpow(Rs[km], m - n));
    }

    return interpolate(Rs[km], Rhs[km], Vl, Vu, m, Ll, Lu);
}

/**
 * V(s,n,m) for a FINITE piecewise QBD, the port of `mam_transient2.m`.
 *
 * K = |T| - 1 regimes and `Lv[K+1]` is the top boundary level, so the chain is
 * closed above and there is no R tail.
 */
inline CMat mam_transient2(const TransientQbd& q, long n, long m, const Complex& s) {
    using namespace transient_detail;
    if (q.nT() < 2) throw InputError("mam_transient2: the finite form needs at least two thresholds");
    const std::size_t K = q.nT() - 1;
    if (n < 0 || m < 0) throw InputError("mam_transient2: levels must be non-negative");

    std::vector<CMat> Gs(K + 1), Rs(K + 1), Ghs(K + 1), Rhs(K + 1);
    for (std::size_t k = 1; k <= K; ++k) {
        if (q.T1(k + 1) - q.T1(k) == 1) continue;
        const CMat Lk = msub(q.L1(k), sI(s, q.L1(k).rows()));
        const CQbdFundMat gr = qbd_fundmat_laplace(q.B1(k), Lk, q.F1(k));
        Gs[k] = gr.G;
        Rs[k] = gr.R;
        const CQbdFundMat grh = qbd_fundmat_laplace(q.F1(k), Lk, q.B1(k));
        Ghs[k] = grh.G;
        Rhs[k] = grh.R;
    }

    std::vector<CMat> SvHn(K + 1), SvH0(K + 1), SvHhn(K + 1), SvHh0(K + 1);
    for (std::size_t k = 1; k <= K; ++k) {
        const std::size_t NN = q.Lv1(k).rows();
        if (q.T1(k + 1) - q.T1(k) > 1) {
            const HBlocks h = h_blocks(Gs[k], Ghs[k], q.T1(k + 1) - q.T1(k), NN);
            SvHn[k] = h.nn;
            SvH0[k] = h.n0;
            SvHhn[k] = h.hnn;
            SvHh0[k] = h.hn0;
        } else {
            const std::size_t NN1 = q.Lv1(k + 1).rows();
            SvH0[k] = CMat(NN1, NN, Complex(0.0, 0.0));
            SvHh0[k] = eye<Complex>(NN);
            SvHhn[k] = CMat(NN, NN1, Complex(0.0, 0.0));
            SvHn[k] = eye<Complex>(NN1);
        }
    }

    const std::size_t NNK = q.Lv1(K + 1).rows();
    std::vector<CMat> SY(K + 2);
    {
        const CMat M = msub(msub(sI(s, NNK), q.Lv1(K + 1)), matmul(q.B1(K), SvHhn[K]));
        SY[K] = madd(SvH0[K], matmul(SvHn[K], mldivide(M, matmul(q.B1(K), SvHh0[K]))));
    }
    for (std::size_t k = K; k-- > 1;) {
        const std::size_t NNk1 = q.Lv1(k + 1).rows();
        const CMat M = msub(msub(msub(sI(s, NNk1), q.Lv1(k + 1)), matmul(q.F1(k + 1), SY[k + 1])),
                            matmul(q.B1(k), SvHhn[k]));
        SY[k] = madd(SvH0[k], matmul(SvHn[k], mldivide(M, matmul(q.B1(k), SvHh0[k]))));
    }

    const std::size_t NN1 = q.Lv1(1).rows();
    std::vector<CMat> SYh(K + 1);
    {
        const CMat M = msub(msub(sI(s, NN1), q.Lv1(1)), matmul(q.F1(1), SvH0[1]));
        SYh[1] = madd(SvHhn[1], matmul(SvHh0[1], mldivide(M, matmul(q.F1(1), SvHn[1]))));
    }
    for (std::size_t k = 2; k <= K; ++k) {
        const std::size_t NNk = q.Lv1(k).rows();
        const CMat M = msub(msub(msub(sI(s, NNk), q.Lv1(k)), matmul(q.B1(k - 1), SYh[k - 1])),
                            matmul(q.F1(k), SvH0[k]));
        SYh[k] = madd(SvHhn[k], matmul(SvHh0[k], mldivide(M, matmul(q.F1(k), SvHn[k]))));
    }

    std::vector<std::vector<CMat>> SV(K + 3, std::vector<CMat>(K + 3));
    for (std::size_t l = 0; l <= K; ++l) {
        if (l == 0) {
            SV[1][1] = mldivide(msub(msub(sI(s, NN1), q.Lv1(1)), matmul(q.F1(1), SY[1])),
                                eye<Complex>(NN1));
        } else if (l == K) {
            SV[K + 1][K + 1] = mldivide(
                msub(msub(sI(s, NNK), q.Lv1(K + 1)), matmul(q.B1(K), SYh[K])), eye<Complex>(NNK));
        } else {
            const std::size_t NNl1 = q.Lv1(l + 1).rows();
            const CMat M = msub(msub(msub(sI(s, NNl1), q.Lv1(l + 1)), matmul(q.F1(l + 1), SY[l + 1])),
                                matmul(q.B1(l), SYh[l]));
            SV[l + 1][l + 1] = mldivide(M, eye<Complex>(NNl1));
        }
        for (std::size_t k = l + 1; k + 1 <= K; ++k) {
            const std::size_t NNk1 = q.Lv1(k + 1).rows();
            const CMat M = msub(msub(msub(sI(s, NNk1), q.Lv1(k + 1)), matmul(q.F1(k + 1), SY[k + 1])),
                                matmul(q.B1(k), SvHhn[k]));
            SV[k + 1][l + 1] = mldivide(M, matmul(matmul(q.B1(k), SvHh0[k]), SV[k][l + 1]));
        }
        if (l < K) {
            const CMat M = msub(msub(sI(s, NNK), q.Lv1(K + 1)), matmul(q.B1(K), SvHhn[K]));
            SV[K + 1][l + 1] = mldivide(M, matmul(matmul(q.B1(K), SvHh0[K]), SV[K][l + 1]));
        }
        for (std::size_t k = l; k-- > 1;) {
            const std::size_t NNk1 = q.Lv1(k + 1).rows();
            const CMat M = msub(msub(msub(sI(s, NNk1), q.Lv1(k + 1)), matmul(q.F1(k + 1), SvH0[k + 1])),
                                matmul(q.B1(k), SYh[k]));
            SV[k + 1][l + 1] = mldivide(M, matmul(matmul(q.F1(k + 1), SvHn[k + 1]), SV[k + 2][l + 1]));
        }
        if (l > 0) {
            const CMat M = msub(msub(sI(s, NN1), q.Lv1(1)), matmul(q.F1(1), SvH0[1]));
            SV[1][l + 1] = mldivide(M, matmul(matmul(q.F1(1), SvHn[1]), SV[2][l + 1]));
        }
    }

    const std::size_t kn = regime_of(q, n, K + 1);
    const std::size_t km = regime_of(q, m, K + 1);

    const std::size_t NNkn = q.Lv1(kn).rows();
    const CMat IIkn = eye<Complex>(NNkn);

    CMat Vu, Vl;
    long Lu = 0, Ll = 0;

    if (q.T1(kn) == n) {
        if (q.T1(km) == m) return SV[kn][km];
        Vu = SV[kn][km + 1];
        Vl = SV[kn][km];
        Lu = q.T1(km + 1);
        Ll = q.T1(km);
    } else {
        const HBlocks hT = h_blocks(Gs[kn], Ghs[kn], q.T1(kn + 1) - n, NNkn);
        const HBlocks hn = h_blocks(Gs[kn], Ghs[kn], n - q.T1(kn), NNkn);

        CMat Yn;
        if (kn == K) {
            const CMat M = msub(msub(sI(s, NNK), q.Lv1(K + 1)), matmul(q.B1(K), hT.hnn));
            Yn = madd(hT.n0, matmul(hT.nn, mldivide(M, matmul(q.B1(K), hT.hn0))));
        } else {
            const std::size_t NNkn1 = q.Lv1(kn + 1).rows();
            const CMat M = msub(msub(msub(sI(s, NNkn1), q.Lv1(kn + 1)), matmul(q.F1(kn + 1), SY[kn + 1])),
                                matmul(q.B1(kn), hT.hnn));
            Yn = madd(hT.n0, matmul(hT.nn, mldivide(M, matmul(q.B1(kn), hT.hn0))));
        }

        CMat Yhn;
        if (kn == 1) {
            const CMat M = msub(msub(sI(s, NN1), q.Lv1(1)), matmul(q.F1(1), hn.n0));
            Yhn = madd(hn.hnn, matmul(hn.hn0, mldivide(M, matmul(q.F1(1), hn.nn))));
        } else {
            const std::size_t NNk = q.Lv1(kn).rows();
            const CMat M = msub(msub(msub(sI(s, NNk), q.Lv1(kn)), matmul(q.B1(kn - 1), SYh[kn - 1])),
                                matmul(q.F1(kn), hn.n0));
            Yhn = madd(hn.hnn, matmul(hn.hn0, mldivide(M, matmul(q.F1(kn), hn.nn))));
        }

        const CMat Mkn = msub(sI(s, q.L1(kn).rows()), q.L1(kn));
        CMat Vnl;
        if (q.T1(km) < n) {
            const CMat M = msub(msub(Mkn, matmul(q.B1(kn), hn.hnn)), matmul(q.F1(kn), Yn));
            Vnl = mldivide(M, matmul(matmul(q.B1(kn), hn.hn0), SV[kn][km]));
        } else {
            const CMat M = msub(msub(Mkn, matmul(q.F1(kn), hT.n0)), matmul(q.B1(kn), Yhn));
            Vnl = mldivide(M, matmul(matmul(q.F1(kn), hT.nn), SV[kn + 1][km]));
        }
        if (m == q.T1(km)) return Vnl;

        CMat Vnu;
        if (q.T1(km + 1) < n) {
            const CMat M = msub(msub(Mkn, matmul(q.B1(kn), hn.hnn)), matmul(q.F1(kn), Yn));
            Vnu = mldivide(M, matmul(matmul(q.B1(kn), hn.hn0), SV[kn][km + 1]));
        } else {
            const CMat M = msub(msub(Mkn, matmul(q.F1(kn), hT.n0)), matmul(q.B1(kn), Yhn));
            Vnu = mldivide(M, matmul(matmul(q.F1(kn), hT.nn), SV[kn + 1][km + 1]));
        }
        const CMat Vnn = mldivide(
            msub(msub(msub(sI(s, NNkn), q.L1(kn)), matmul(q.B1(kn), Yhn)), matmul(q.F1(kn), Yn)),
            IIkn);
        if (n == m) return Vnn;
        if (kn != km) {
            Vu = Vnu; Vl = Vnl; Lu = q.T1(km + 1); Ll = q.T1(km);
        } else if (n <= m) {
            Vu = Vnu; Vl = Vnn; Lu = q.T1(km + 1); Ll = n;
        } else {
            Vu = Vnn; Vl = Vnl; Lu = n; Ll = q.T1(km);
        }
    }

    return interpolate(Rs[km], Rhs[km], Vl, Vu, m, Ll, Lu);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAM_TRANSIENT2_H
