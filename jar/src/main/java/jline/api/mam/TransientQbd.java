/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 *
 * Laplace-domain transient analysis of piecewise level-dependent QBD processes.
 * Port of the MATLAB transient-QBD engine (mam_transient2_open / mam_transient2
 * / qbd_fundmat). Block arrays are 1-based (index 0 unused) to mirror the
 * reference indexing exactly. All arithmetic is complex (argument s).
 *
 * Uses commons-math3 dense FieldMatrix<Complex> (not jline's sparse-backed
 * ComplexMatrix): the transient level blocks are small and dense, and dense
 * arithmetic avoids the large sparse-CSC overhead per operation.
 */
package jline.api.mam;

import jline.util.matrix.Matrix;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexField;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

public final class TransientQbd {

    private TransientQbd() {}

    private static final ComplexField CF = ComplexField.getInstance();

    // ---------------------------------------------------------------------
    // Dense complex-matrix helpers
    // ---------------------------------------------------------------------

    /** Wrap a real jline Matrix as a dense complex FieldMatrix. */
    public static FieldMatrix<Complex> cm(Matrix real) {
        int r = real.getNumRows();
        int c = real.getNumCols();
        Complex[][] d = new Complex[r][c];
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                d[i][j] = new Complex(real.get(i, j));
            }
        }
        return new Array2DRowFieldMatrix<Complex>(d, false);
    }

    public static FieldMatrix<Complex> eye(int n) {
        return MatrixUtils.createFieldIdentityMatrix(CF, n);
    }

    static FieldMatrix<Complex> zeros(int r, int c) {
        Complex[][] d = new Complex[r][c];
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                d[i][j] = Complex.ZERO;
            }
        }
        return new Array2DRowFieldMatrix<Complex>(d, false);
    }

    public static FieldMatrix<Complex> sEye(int n, Complex s) {
        return eye(n).scalarMultiply(s);
    }

    /** MATLAB A \ B = inv(A) * B (A square). */
    public static FieldMatrix<Complex> ldiv(FieldMatrix<Complex> A, FieldMatrix<Complex> B) {
        return new FieldLUDecomposition<Complex>(A).getSolver().solve(B);
    }

    /** MATLAB X / Y = X * inv(Y), solved as (Y^T \ X^T)^T (plain transpose). */
    static FieldMatrix<Complex> rdiv(FieldMatrix<Complex> X, FieldMatrix<Complex> Y) {
        return ldiv(Y.transpose(), X.transpose()).transpose();
    }

    static FieldMatrix<Complex> mxpow(FieldMatrix<Complex> A, int k) {
        if (k == 0) {
            return eye(A.getRowDimension());
        }
        FieldMatrix<Complex> P = A.copy();
        for (int i = 1; i < k; i++) {
            P = P.multiply(A);
        }
        return P;
    }

    /** Assemble [[A, B], [C, D]]. */
    static FieldMatrix<Complex> blk22(FieldMatrix<Complex> A, FieldMatrix<Complex> B,
                                      FieldMatrix<Complex> C, FieldMatrix<Complex> D) {
        int r1 = A.getRowDimension(), c1 = A.getColumnDimension();
        int r2 = C.getRowDimension(), c2 = B.getColumnDimension();
        FieldMatrix<Complex> out = zeros(r1 + r2, c1 + c2);
        out.setSubMatrix(A.getData(), 0, 0);
        out.setSubMatrix(B.getData(), 0, c1);
        out.setSubMatrix(C.getData(), r1, 0);
        out.setSubMatrix(D.getData(), r1, c1);
        return out;
    }

    /** Submatrix rows [r0,r1) x cols [c0,c1). */
    static FieldMatrix<Complex> sub(FieldMatrix<Complex> M, int r0, int r1, int c0, int c1) {
        return M.getSubMatrix(r0, r1 - 1, c0, c1 - 1);
    }

    private static double infNorm(FieldMatrix<Complex> A) {
        double best = 0.0;
        for (int i = 0; i < A.getRowDimension(); i++) {
            double rs = 0.0;
            for (int j = 0; j < A.getColumnDimension(); j++) {
                rs += A.getEntry(i, j).abs();
            }
            best = Math.max(best, rs);
        }
        return best;
    }

    // ---------------------------------------------------------------------
    // Fundamental matrices (cyclic reduction), complex-capable
    // ---------------------------------------------------------------------

    /** Returns {G, R} of the QBD with raw blocks (B, L, F). */
    public static FieldMatrix<Complex>[] qbdFundMat(FieldMatrix<Complex> B, FieldMatrix<Complex> L,
                                                    FieldMatrix<Complex> F) {
        int m = L.getRowDimension();
        FieldMatrix<Complex> II = eye(m);
        double lamb = 0.0;
        for (int i = 0; i < m; i++) {
            lamb = Math.max(lamb, -L.getEntry(i, i).getReal());
        }
        Complex invLamb = new Complex(1.0 / lamb);
        FieldMatrix<Complex> Bm = B.scalarMultiply(invLamb);
        FieldMatrix<Complex> Lm = L.scalarMultiply(invLamb).add(II);
        FieldMatrix<Complex> Fm = F.scalarMultiply(invLamb);

        FieldMatrix<Complex> BF = ldiv(II.subtract(Lm), II);
        FieldMatrix<Complex> BB = BF.multiply(Fm);
        BF = BF.multiply(Bm);
        FieldMatrix<Complex> G = BF.copy();
        FieldMatrix<Complex> PI = BB.copy();
        double check = 1.0;
        int numit = 0;
        while (check > 1e-14 && numit < 50) {
            FieldMatrix<Complex> Lstar = BF.multiply(BB).add(BB.multiply(BF));
            FieldMatrix<Complex> Bstar = BB.multiply(BB);
            FieldMatrix<Complex> Fstar = BF.multiply(BF);
            BB = ldiv(II.subtract(Lstar), II);
            BF = BB.multiply(Fstar);
            BB = BB.multiply(Bstar);
            G = G.add(PI.multiply(BF));
            PI = PI.multiply(BB);
            check = Math.min(infNorm(BB), infNorm(BF));
            numit++;
        }
        FieldMatrix<Complex> R = Fm.multiply(ldiv(II.subtract(Lm.add(Fm.multiply(G))), II));
        @SuppressWarnings("unchecked")
        FieldMatrix<Complex>[] out = new FieldMatrix[]{G, R};
        return out;
    }

    // ---------------------------------------------------------------------
    // Transient transform: open (infinite) QBD
    // ---------------------------------------------------------------------

    /**
     * Laplace-domain transient V(s,n,m) for an open QBD. Block arrays are
     * 1-based (index 0 unused); K = T.length - 1 regimes, last repeats.
     */
    public static FieldMatrix<Complex> transient2Open(FieldMatrix<Complex>[] B, FieldMatrix<Complex>[] L,
            FieldMatrix<Complex>[] F, FieldMatrix<Complex>[] Lv, int[] T, int n, int m, Complex s) {
        int K = T.length - 1;

        FieldMatrix<Complex>[] Gs = newArr(K + 1);
        FieldMatrix<Complex>[] Rs = newArr(K + 1);
        FieldMatrix<Complex>[] Ghs = newArr(K + 1);
        FieldMatrix<Complex>[] Rhs = newArr(K + 1);
        for (int k = 1; k <= K; k++) {
            if (k < K && T[k + 1] - T[k] == 1) {
                continue;
            }
            FieldMatrix<Complex> Lk = L[k].subtract(sEye(L[k].getRowDimension(), s));
            FieldMatrix<Complex>[] gr = qbdFundMat(B[k], Lk, F[k]);
            Gs[k] = gr[0];
            Rs[k] = gr[1];
            FieldMatrix<Complex>[] grh = qbdFundMat(F[k], Lk, B[k]);
            Ghs[k] = grh[0];
            Rhs[k] = grh[1];
        }

        FieldMatrix<Complex>[] SvHn = newArr(K + 1);
        FieldMatrix<Complex>[] SvH0 = newArr(K + 1);
        FieldMatrix<Complex>[] SvHhn = newArr(K + 1);
        FieldMatrix<Complex>[] SvHh0 = newArr(K + 1);
        for (int k = 1; k < K; k++) {
            int NN = Lv[k].getRowDimension();
            if (T[k + 1] - T[k] > 1) {
                int d = T[k + 1] - T[k];
                FieldMatrix<Complex> num = blk22(mxpow(Ghs[k], d - 1), Gs[k], Ghs[k], mxpow(Gs[k], d - 1));
                FieldMatrix<Complex> den = blk22(eye(NN), mxpow(Gs[k], d), mxpow(Ghs[k], d), eye(NN));
                FieldMatrix<Complex> SH = rdiv(num, den);
                SvHn[k] = sub(SH, 0, NN, 0, NN);
                SvH0[k] = sub(SH, 0, NN, NN, 2 * NN);
                SvHhn[k] = sub(SH, NN, 2 * NN, 0, NN);
                SvHh0[k] = sub(SH, NN, 2 * NN, NN, 2 * NN);
            } else {
                int NN1 = Lv[k + 1].getRowDimension();
                SvH0[k] = zeros(NN1, NN);
                SvHh0[k] = eye(NN);
                SvHhn[k] = zeros(NN, NN1);
                SvHn[k] = eye(NN1);
            }
        }

        FieldMatrix<Complex>[] SY = newArr(K + 1);
        SY[K] = Gs[K];
        for (int k = K - 1; k >= 1; k--) {
            int NNk1 = Lv[k + 1].getRowDimension();
            FieldMatrix<Complex> Aik = sEye(NNk1, s).subtract(Lv[k + 1]).subtract(F[k + 1].multiply(SY[k + 1])).subtract(B[k].multiply(SvHhn[k]));
            SY[k] = SvH0[k].add(SvHn[k].multiply(ldiv(Aik, B[k].multiply(SvHh0[k]))));
        }

        int NN1 = Lv[1].getRowDimension();
        FieldMatrix<Complex>[] SYh = newArr(K + 1);
        FieldMatrix<Complex> A1 = sEye(NN1, s).subtract(Lv[1]).subtract(F[1].multiply(SvH0[1]));
        SYh[1] = SvHhn[1].add(SvHh0[1].multiply(ldiv(A1, F[1].multiply(SvHn[1]))));
        for (int k = 2; k < K; k++) {
            int NNk = Lv[k].getRowDimension();
            FieldMatrix<Complex> Ak = sEye(NNk, s).subtract(Lv[k]).subtract(B[k - 1].multiply(SYh[k - 1])).subtract(F[k].multiply(SvH0[k]));
            SYh[k] = SvHhn[k].add(SvHh0[k].multiply(ldiv(Ak, F[k].multiply(SvHn[k]))));
        }

        @SuppressWarnings("unchecked") // generic array creation is an inherent Java limitation
        FieldMatrix<Complex>[][] SV = new FieldMatrix[K + 1][K + 1];
        for (int l = 0; l <= K - 1; l++) {
            if (l == 0) {
                FieldMatrix<Complex> A = sEye(NN1, s).subtract(Lv[1]).subtract(F[1].multiply(SY[1]));
                SV[1][1] = ldiv(A, eye(NN1));
            } else {
                int NNl1 = Lv[l + 1].getRowDimension();
                FieldMatrix<Complex> A = sEye(NNl1, s).subtract(Lv[l + 1]).subtract(F[l + 1].multiply(SY[l + 1])).subtract(B[l].multiply(SYh[l]));
                SV[l + 1][l + 1] = ldiv(A, eye(NNl1));
            }
            for (int k = l + 1; k <= K - 1; k++) {
                int NNk1 = Lv[k + 1].getRowDimension();
                FieldMatrix<Complex> A = sEye(NNk1, s).subtract(Lv[k + 1]).subtract(F[k + 1].multiply(SY[k + 1])).subtract(B[k].multiply(SvHhn[k]));
                SV[k + 1][l + 1] = ldiv(A, B[k].multiply(SvHh0[k]).multiply(SV[k][l + 1]));
            }
            for (int k = l - 1; k >= 1; k--) {
                int NNk1 = Lv[k + 1].getRowDimension();
                FieldMatrix<Complex> A = sEye(NNk1, s).subtract(Lv[k + 1]).subtract(F[k + 1].multiply(SvH0[k + 1])).subtract(B[k].multiply(SYh[k]));
                SV[k + 1][l + 1] = ldiv(A, F[k + 1].multiply(SvHn[k + 1]).multiply(SV[k + 2][l + 1]));
            }
            if (l > 0) {
                FieldMatrix<Complex> A = sEye(NN1, s).subtract(Lv[1]).subtract(F[1].multiply(SvH0[1]));
                SV[1][l + 1] = ldiv(A, F[1].multiply(SvHn[1]).multiply(SV[2][l + 1]));
            }
        }

        int kn = regimeOf(T, n, K);
        int km = regimeOf(T, m, K);

        int NN = Lv[kn].getRowDimension();
        FieldMatrix<Complex> II = eye(NN);

        FieldMatrix<Complex> Vu = null, Vl;
        int Lu = -1, Ll;
        boolean thresholdN = (T[kn] == n);
        if (thresholdN) {
            if (T[km] == m) {
                return SV[kn][km];
            }
            Vl = SV[kn][km];
            Ll = T[km];
            if (km < K) {
                Vu = SV[kn][km + 1];
                Lu = T[km + 1];
            }
        } else {
            FieldMatrix<Complex> HTnn = null, HTn0 = null, HhTnn = null, HhTn0 = null;
            FieldMatrix<Complex> HnTn, HnT0, HhnTn, HhnT0, Yn;
            if (kn < K) {
                int d1 = T[kn + 1] - n;
                FieldMatrix<Complex> num = blk22(mxpow(Ghs[kn], d1 - 1), Gs[kn], Ghs[kn], mxpow(Gs[kn], d1 - 1));
                FieldMatrix<Complex> den = blk22(II, mxpow(Gs[kn], d1), mxpow(Ghs[kn], d1), II);
                FieldMatrix<Complex> Tmp = rdiv(num, den);
                HTnn = sub(Tmp, 0, NN, 0, NN);
                HTn0 = sub(Tmp, 0, NN, NN, 2 * NN);
                HhTnn = sub(Tmp, NN, 2 * NN, 0, NN);
                HhTn0 = sub(Tmp, NN, 2 * NN, NN, 2 * NN);

                int d2 = n - T[kn];
                num = blk22(mxpow(Ghs[kn], d2 - 1), Gs[kn], Ghs[kn], mxpow(Gs[kn], d2 - 1));
                den = blk22(II, mxpow(Gs[kn], d2), mxpow(Ghs[kn], d2), II);
                Tmp = rdiv(num, den);
                HnTn = sub(Tmp, 0, NN, 0, NN);
                HnT0 = sub(Tmp, 0, NN, NN, 2 * NN);
                HhnTn = sub(Tmp, NN, 2 * NN, 0, NN);
                HhnT0 = sub(Tmp, NN, 2 * NN, NN, 2 * NN);

                int NNkn1 = Lv[kn + 1].getRowDimension();
                FieldMatrix<Complex> A = sEye(NNkn1, s).subtract(Lv[kn + 1]).subtract(F[kn + 1].multiply(SY[kn + 1])).subtract(B[kn].multiply(HhTnn));
                Yn = HTn0.add(HTnn.multiply(ldiv(A, B[kn].multiply(HhTn0))));
            } else {
                int d2 = n - T[kn];
                FieldMatrix<Complex> num = blk22(mxpow(Ghs[kn], d2 - 1), Gs[kn], Ghs[kn], mxpow(Gs[kn], d2 - 1));
                FieldMatrix<Complex> den = blk22(II, mxpow(Gs[kn], d2), mxpow(Ghs[kn], d2), II);
                FieldMatrix<Complex> Tmp = rdiv(num, den);
                HnTn = sub(Tmp, 0, NN, 0, NN);
                HnT0 = sub(Tmp, 0, NN, NN, 2 * NN);
                HhnTn = sub(Tmp, NN, 2 * NN, 0, NN);
                HhnT0 = sub(Tmp, NN, 2 * NN, NN, 2 * NN);
                Yn = Gs[kn];
            }

            FieldMatrix<Complex> Yhn;
            if (kn == 1) {
                FieldMatrix<Complex> A = sEye(NN1, s).subtract(Lv[1]).subtract(F[1].multiply(HnT0));
                Yhn = HhnTn.add(HhnT0.multiply(ldiv(A, F[1].multiply(HnTn))));
            } else {
                int NNkn = Lv[kn].getRowDimension();
                FieldMatrix<Complex> A = sEye(NNkn, s).subtract(Lv[kn]).subtract(B[kn - 1].multiply(SYh[kn - 1])).subtract(F[kn].multiply(HnT0));
                Yhn = HhnTn.add(HhnT0.multiply(ldiv(A, F[kn].multiply(HnTn))));
            }

            FieldMatrix<Complex> Mkn = sEye(L[kn].getRowDimension(), s).subtract(L[kn]);
            FieldMatrix<Complex> Vnl;
            if (T[km] < n) {
                Vnl = ldiv(Mkn.subtract(B[kn].multiply(HhnTn)).subtract(F[kn].multiply(Yn)), B[kn].multiply(HhnT0).multiply(SV[kn][km]));
            } else {
                Vnl = ldiv(Mkn.subtract(F[kn].multiply(HTn0)).subtract(B[kn].multiply(Yhn)), F[kn].multiply(HTnn).multiply(SV[kn + 1][km]));
            }
            if (m == T[km]) {
                return Vnl;
            }

            FieldMatrix<Complex> Vnu = null;
            if (km < K) {
                if (T[km + 1] < n) {
                    Vnu = ldiv(Mkn.subtract(B[kn].multiply(HhnTn)).subtract(F[kn].multiply(Yn)), B[kn].multiply(HhnT0).multiply(SV[kn][km + 1]));
                } else {
                    Vnu = ldiv(Mkn.subtract(F[kn].multiply(HTn0)).subtract(B[kn].multiply(Yhn)), F[kn].multiply(HTnn).multiply(SV[kn + 1][km + 1]));
                }
            }
            FieldMatrix<Complex> Vnn = ldiv(sEye(NN, s).subtract(L[kn]).subtract(B[kn].multiply(Yhn)).subtract(F[kn].multiply(Yn)), II);
            if (n == m) {
                return Vnn;
            }
            if (km == K && n < m) {
                if (kn < K) {
                    return Vnl.multiply(mxpow(Rs[km], m - T[km]));
                }
                return Vnn.multiply(mxpow(Rs[km], m - n));
            }
            if (kn != km) {
                Vu = Vnu; Vl = Vnl; Lu = T[km + 1]; Ll = T[km];
            } else if (n <= m) {
                Vu = Vnu; Vl = Vnn; Lu = T[km + 1]; Ll = n;
            } else {
                Vu = Vnn; Vl = Vnl; Lu = n; Ll = T[km];
            }
        }

        if (km == K && n < m) {
            if (kn < K) {
                return Vl.multiply(mxpow(Rs[km], m - T[km]));
            }
            return SV[kn][kn].multiply(mxpow(Rs[km], m - n));
        }

        NN = Rs[km].getRowDimension();
        II = eye(NN);
        FieldMatrix<Complex> Zden = blk22(II, mxpow(Rs[km], Lu - Ll), mxpow(Rhs[km], Lu - Ll), II);
        FieldMatrix<Complex> Znum = zeros(2 * NN, NN);
        Znum.setSubMatrix(mxpow(Rs[km], m - Ll).getData(), 0, 0);
        Znum.setSubMatrix(mxpow(Rhs[km], Lu - m).getData(), NN, 0);
        FieldMatrix<Complex> Z = ldiv(Zden, Znum);
        return Vl.multiply(sub(Z, 0, NN, 0, NN)).add(Vu.multiply(sub(Z, NN, 2 * NN, 0, NN)));
    }

    // ---------------------------------------------------------------------
    // Transient transform: closed (finite) QBD
    // ---------------------------------------------------------------------

    /**
     * Laplace-domain transient V(s,n,m) for a finite QBD. Block arrays are
     * 1-based; T is 1-based of length K+2 with Lv[K+1] the top boundary level.
     */
    public static FieldMatrix<Complex> transient2(FieldMatrix<Complex>[] B, FieldMatrix<Complex>[] L,
            FieldMatrix<Complex>[] F, FieldMatrix<Complex>[] Lv, int[] T, int n, int m, Complex s) {
        int K = T.length - 2;

        FieldMatrix<Complex>[] Gs = newArr(K + 1);
        FieldMatrix<Complex>[] Rs = newArr(K + 1);
        FieldMatrix<Complex>[] Ghs = newArr(K + 1);
        FieldMatrix<Complex>[] Rhs = newArr(K + 1);
        for (int k = 1; k <= K; k++) {
            if (T[k + 1] - T[k] == 1) {
                continue;
            }
            FieldMatrix<Complex> Lk = L[k].subtract(sEye(L[k].getRowDimension(), s));
            FieldMatrix<Complex>[] gr = qbdFundMat(B[k], Lk, F[k]);
            Gs[k] = gr[0];
            Rs[k] = gr[1];
            FieldMatrix<Complex>[] grh = qbdFundMat(F[k], Lk, B[k]);
            Ghs[k] = grh[0];
            Rhs[k] = grh[1];
        }

        FieldMatrix<Complex>[] SvHn = newArr(K + 2);
        FieldMatrix<Complex>[] SvH0 = newArr(K + 2);
        FieldMatrix<Complex>[] SvHhn = newArr(K + 2);
        FieldMatrix<Complex>[] SvHh0 = newArr(K + 2);
        for (int k = 1; k <= K; k++) {
            int NN = Lv[k].getRowDimension();
            if (T[k + 1] - T[k] > 1) {
                int d = T[k + 1] - T[k];
                FieldMatrix<Complex> num = blk22(mxpow(Ghs[k], d - 1), Gs[k], Ghs[k], mxpow(Gs[k], d - 1));
                FieldMatrix<Complex> den = blk22(eye(NN), mxpow(Gs[k], d), mxpow(Ghs[k], d), eye(NN));
                FieldMatrix<Complex> SH = rdiv(num, den);
                SvHn[k] = sub(SH, 0, NN, 0, NN);
                SvH0[k] = sub(SH, 0, NN, NN, 2 * NN);
                SvHhn[k] = sub(SH, NN, 2 * NN, 0, NN);
                SvHh0[k] = sub(SH, NN, 2 * NN, NN, 2 * NN);
            } else {
                int NN1 = Lv[k + 1].getRowDimension();
                SvH0[k] = zeros(NN1, NN);
                SvHh0[k] = eye(NN);
                SvHhn[k] = zeros(NN, NN1);
                SvHn[k] = eye(NN1);
            }
        }

        int NNK = Lv[K + 1].getRowDimension();
        FieldMatrix<Complex>[] SY = newArr(K + 1);
        FieldMatrix<Complex> AK = sEye(NNK, s).subtract(Lv[K + 1]).subtract(B[K].multiply(SvHhn[K]));
        SY[K] = SvH0[K].add(SvHn[K].multiply(ldiv(AK, B[K].multiply(SvHh0[K]))));
        for (int k = K - 1; k >= 1; k--) {
            int NNk1 = Lv[k + 1].getRowDimension();
            FieldMatrix<Complex> A = sEye(NNk1, s).subtract(Lv[k + 1]).subtract(F[k + 1].multiply(SY[k + 1])).subtract(B[k].multiply(SvHhn[k]));
            SY[k] = SvH0[k].add(SvHn[k].multiply(ldiv(A, B[k].multiply(SvHh0[k]))));
        }

        int NN1 = Lv[1].getRowDimension();
        FieldMatrix<Complex>[] SYh = newArr(K + 1);
        FieldMatrix<Complex> A1 = sEye(NN1, s).subtract(Lv[1]).subtract(F[1].multiply(SvH0[1]));
        SYh[1] = SvHhn[1].add(SvHh0[1].multiply(ldiv(A1, F[1].multiply(SvHn[1]))));
        for (int k = 2; k <= K; k++) {
            int NNk = Lv[k].getRowDimension();
            FieldMatrix<Complex> A = sEye(NNk, s).subtract(Lv[k]).subtract(B[k - 1].multiply(SYh[k - 1])).subtract(F[k].multiply(SvH0[k]));
            SYh[k] = SvHhn[k].add(SvHh0[k].multiply(ldiv(A, F[k].multiply(SvHn[k]))));
        }

        @SuppressWarnings("unchecked") // generic array creation is an inherent Java limitation
        FieldMatrix<Complex>[][] SV = new FieldMatrix[K + 2][K + 2];
        for (int l = 0; l <= K; l++) {
            if (l == 0) {
                FieldMatrix<Complex> A = sEye(NN1, s).subtract(Lv[1]).subtract(F[1].multiply(SY[1]));
                SV[1][1] = ldiv(A, eye(NN1));
            } else if (l == K) {
                FieldMatrix<Complex> A = sEye(NNK, s).subtract(Lv[K + 1]).subtract(B[K].multiply(SYh[K]));
                SV[K + 1][K + 1] = ldiv(A, eye(NNK));
            } else {
                int NNl1 = Lv[l + 1].getRowDimension();
                FieldMatrix<Complex> A = sEye(NNl1, s).subtract(Lv[l + 1]).subtract(F[l + 1].multiply(SY[l + 1])).subtract(B[l].multiply(SYh[l]));
                SV[l + 1][l + 1] = ldiv(A, eye(NNl1));
            }
            for (int k = l + 1; k <= K - 1; k++) {
                int NNk1 = Lv[k + 1].getRowDimension();
                FieldMatrix<Complex> A = sEye(NNk1, s).subtract(Lv[k + 1]).subtract(F[k + 1].multiply(SY[k + 1])).subtract(B[k].multiply(SvHhn[k]));
                SV[k + 1][l + 1] = ldiv(A, B[k].multiply(SvHh0[k]).multiply(SV[k][l + 1]));
            }
            if (l < K) {
                FieldMatrix<Complex> A = sEye(NNK, s).subtract(Lv[K + 1]).subtract(B[K].multiply(SvHhn[K]));
                SV[K + 1][l + 1] = ldiv(A, B[K].multiply(SvHh0[K]).multiply(SV[K][l + 1]));
            }
            for (int k = l - 1; k >= 1; k--) {
                int NNk1 = Lv[k + 1].getRowDimension();
                FieldMatrix<Complex> A = sEye(NNk1, s).subtract(Lv[k + 1]).subtract(F[k + 1].multiply(SvH0[k + 1])).subtract(B[k].multiply(SYh[k]));
                SV[k + 1][l + 1] = ldiv(A, F[k + 1].multiply(SvHn[k + 1]).multiply(SV[k + 2][l + 1]));
            }
            if (l > 0) {
                FieldMatrix<Complex> A = sEye(NN1, s).subtract(Lv[1]).subtract(F[1].multiply(SvH0[1]));
                SV[1][l + 1] = ldiv(A, F[1].multiply(SvHn[1]).multiply(SV[2][l + 1]));
            }
        }

        int kn = regimeOf(T, n, K + 1);
        int km = regimeOf(T, m, K + 1);

        int NN = Lv[kn].getRowDimension();
        FieldMatrix<Complex> II = eye(NN);

        FieldMatrix<Complex> Vu, Vl;
        int Lu, Ll;
        if (T[kn] == n) {
            if (T[km] == m) {
                return SV[kn][km];
            }
            Vu = SV[kn][km + 1];
            Vl = SV[kn][km];
            Lu = T[km + 1];
            Ll = T[km];
        } else {
            int d1 = T[kn + 1] - n;
            FieldMatrix<Complex> num = blk22(mxpow(Ghs[kn], d1 - 1), Gs[kn], Ghs[kn], mxpow(Gs[kn], d1 - 1));
            FieldMatrix<Complex> den = blk22(II, mxpow(Gs[kn], d1), mxpow(Ghs[kn], d1), II);
            FieldMatrix<Complex> Tmp = rdiv(num, den);
            FieldMatrix<Complex> HTnn = sub(Tmp, 0, NN, 0, NN);
            FieldMatrix<Complex> HTn0 = sub(Tmp, 0, NN, NN, 2 * NN);
            FieldMatrix<Complex> HhTnn = sub(Tmp, NN, 2 * NN, 0, NN);
            FieldMatrix<Complex> HhTn0 = sub(Tmp, NN, 2 * NN, NN, 2 * NN);

            int d2 = n - T[kn];
            num = blk22(mxpow(Ghs[kn], d2 - 1), Gs[kn], Ghs[kn], mxpow(Gs[kn], d2 - 1));
            den = blk22(II, mxpow(Gs[kn], d2), mxpow(Ghs[kn], d2), II);
            Tmp = rdiv(num, den);
            FieldMatrix<Complex> HnTn = sub(Tmp, 0, NN, 0, NN);
            FieldMatrix<Complex> HnT0 = sub(Tmp, 0, NN, NN, 2 * NN);
            FieldMatrix<Complex> HhnTn = sub(Tmp, NN, 2 * NN, 0, NN);
            FieldMatrix<Complex> HhnT0 = sub(Tmp, NN, 2 * NN, NN, 2 * NN);

            FieldMatrix<Complex> Yn;
            if (kn == K) {
                FieldMatrix<Complex> A = sEye(NNK, s).subtract(Lv[K + 1]).subtract(B[K].multiply(HhTnn));
                Yn = HTn0.add(HTnn.multiply(ldiv(A, B[K].multiply(HhTn0))));
            } else {
                int NNkn1 = Lv[kn + 1].getRowDimension();
                FieldMatrix<Complex> A = sEye(NNkn1, s).subtract(Lv[kn + 1]).subtract(F[kn + 1].multiply(SY[kn + 1])).subtract(B[kn].multiply(HhTnn));
                Yn = HTn0.add(HTnn.multiply(ldiv(A, B[kn].multiply(HhTn0))));
            }
            FieldMatrix<Complex> Yhn;
            if (kn == 1) {
                FieldMatrix<Complex> A = sEye(NN1, s).subtract(Lv[1]).subtract(F[1].multiply(HnT0));
                Yhn = HhnTn.add(HhnT0.multiply(ldiv(A, F[1].multiply(HnTn))));
            } else {
                int NNkn = Lv[kn].getRowDimension();
                FieldMatrix<Complex> A = sEye(NNkn, s).subtract(Lv[kn]).subtract(B[kn - 1].multiply(SYh[kn - 1])).subtract(F[kn].multiply(HnT0));
                Yhn = HhnTn.add(HhnT0.multiply(ldiv(A, F[kn].multiply(HnTn))));
            }

            FieldMatrix<Complex> Mkn = sEye(L[kn].getRowDimension(), s).subtract(L[kn]);
            FieldMatrix<Complex> Vnl;
            if (T[km] < n) {
                Vnl = ldiv(Mkn.subtract(B[kn].multiply(HhnTn)).subtract(F[kn].multiply(Yn)), B[kn].multiply(HhnT0).multiply(SV[kn][km]));
            } else {
                Vnl = ldiv(Mkn.subtract(F[kn].multiply(HTn0)).subtract(B[kn].multiply(Yhn)), F[kn].multiply(HTnn).multiply(SV[kn + 1][km]));
            }
            if (m == T[km]) {
                return Vnl;
            }
            FieldMatrix<Complex> Vnu;
            if (T[km + 1] < n) {
                Vnu = ldiv(Mkn.subtract(B[kn].multiply(HhnTn)).subtract(F[kn].multiply(Yn)), B[kn].multiply(HhnT0).multiply(SV[kn][km + 1]));
            } else {
                Vnu = ldiv(Mkn.subtract(F[kn].multiply(HTn0)).subtract(B[kn].multiply(Yhn)), F[kn].multiply(HTnn).multiply(SV[kn + 1][km + 1]));
            }
            FieldMatrix<Complex> Vnn = ldiv(sEye(NN, s).subtract(L[kn]).subtract(B[kn].multiply(Yhn)).subtract(F[kn].multiply(Yn)), II);
            if (n == m) {
                return Vnn;
            }
            if (kn != km) {
                Vu = Vnu; Vl = Vnl; Lu = T[km + 1]; Ll = T[km];
            } else if (n <= m) {
                Vu = Vnu; Vl = Vnn; Lu = T[km + 1]; Ll = n;
            } else {
                Vu = Vnn; Vl = Vnl; Lu = n; Ll = T[km];
            }
        }

        NN = Rs[km].getRowDimension();
        II = eye(NN);
        FieldMatrix<Complex> Zden = blk22(II, mxpow(Rs[km], Lu - Ll), mxpow(Rhs[km], Lu - Ll), II);
        FieldMatrix<Complex> Znum = zeros(2 * NN, NN);
        Znum.setSubMatrix(mxpow(Rs[km], m - Ll).getData(), 0, 0);
        Znum.setSubMatrix(mxpow(Rhs[km], Lu - m).getData(), NN, 0);
        FieldMatrix<Complex> Z = ldiv(Zden, Znum);
        return Vl.multiply(sub(Z, 0, NN, 0, NN)).add(Vu.multiply(sub(Z, NN, 2 * NN, 0, NN)));
    }

    @SuppressWarnings("unchecked")
    private static FieldMatrix<Complex>[] newArr(int n) {
        return new FieldMatrix[n];
    }

    /** 1-based regime index k with T[k] <= level < T[k+1], clamped to kmax. */
    static int regimeOf(int[] T, int level, int kmax) {
        for (int k = 1; k < T.length; k++) {
            if (T[k] > level) {
                return k - 1;
            }
        }
        return kmax;
    }
}
