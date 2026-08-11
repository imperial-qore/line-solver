/*
 * Ported from BUTools-family fluid tools (G. Horvath).
 *
 * First- and second-order (Brownian) level-dependent, multi-regime Markovian
 * fluid queues. The generator, drift and (optionally) variance change at
 * threshold fluid levels; setting the variance cells to zero reduces the model
 * to first order.
 */
package jline.lib.butools.mam;

import java.util.ArrayList;
import java.util.List;

import jline.lib.butools.QBDFundamentalMatrices;
import jline.util.matrix.Matrix;

import static jline.lib.butools.mam.FluidTools.*;

public final class SecondOrderLevelDependentFluidSolve {
    private SecondOrderLevelDependentFluidSolve() {}

    public static LevelDependentFluidSolution solve(List<Matrix> Q, List<Matrix> R, List<Matrix> S,
                                                    double[] T, double[] boundaryL, double[] boundaryU,
                                                    List<Matrix> Qt, double prec) {
        int K = T.length;
        int N = Q.get(0).getNumRows();
        if (boundaryL == null) boundaryL = new double[N];
        if (boundaryU == null) boundaryU = boundaryL;

        if (Qt == null || Qt.isEmpty()) {
            Qt = new ArrayList<Matrix>();
            for (int k = 0; k < K; k++) Qt.add(Q.get(k));
            Qt.add(Q.get(K - 1));
        } else if (Qt.size() == 1) {
            Matrix q0 = Qt.get(0);
            Qt = new ArrayList<Matrix>();
            for (int k = 0; k < K + 1; k++) Qt.add(q0);
        }

        double[] Tarr = new double[K + 1];
        for (int k = 0; k < K; k++) Tarr[k + 1] = T[k];

        List<Matrix> KF = nulls(K), KB = nulls(K), cloF = nulls(K), cloB = nulls(K);
        int[] Np = new int[K], Nn = new int[K], Ns = new int[K], NbF = new int[K], NbB = new int[K];
        int[][] vixp = new int[K][], vixn = new int[K][], vix0 = new int[K][], vixs = new int[K][];

        for (int k = 0; k < K; k++) {
            Matrix Sk = S.get(k).scale(0.5);
            Matrix Qk = Q.get(k), Rk = R.get(k);
            double[] dR = diagArray(Rk), dS = diagArray(Sk);
            int[] ix0 = where2(dR, dS, prec, true);       // |R|<=prec && S<=prec
            int[] ixn0 = complement(ix0, N);

            Matrix Q00 = sub(Qk, ix0, ix0), Q0n = sub(Qk, ix0, ixn0);
            Matrix Qn0 = sub(Qk, ixn0, ix0), Qnn = sub(Qk, ixn0, ixn0);
            Matrix Qv = ix0.length > 0 ? Qnn.add(Qn0.mult(Q00.scale(-1.0).pinv()).mult(Q0n)) : Qnn;
            Matrix Rv = sub(Rk, ixn0, ixn0), Sv = sub(Sk, ixn0, ixn0);
            int Nnz = Qv.getNumRows();

            double[] dRv = diagArray(Rv), dSv = diagArray(Sv);
            int[] ixp = whereRS(dRv, dSv, prec, +1);
            int[] ixn = whereRS(dRv, dSv, prec, -1);
            int[] ixs = whereS(dSv, prec);
            Np[k] = ixp.length; Nn[k] = ixn.length; Ns[k] = ixs.length;

            // ----- FORWARD -----
            double c1 = ixp.length > 0 ? maxRatio(neg(diagArray(sub(Qv, ixp, ixp))), diagArray(sub(Rv, ixp, ixp))) : Double.NEGATIVE_INFINITY;
            double c2 = cForward(Rv, Sv, Qv, ixs);
            double c = Math.max(Math.max(c1, c2), 1.0);
            int[] ixbF = concat(ixs, ixp);
            NbF[k] = ixbF.length;
            Matrix Bm = blkdiag(sub(Sv, ixbF, ixbF).scale(c), sub(Rv, ixn, ixn).scale(-1.0));
            Matrix Lm = block(
                    sub(Rv, ixbF, ixbF).scale(-1.0).sub(sub(Sv, ixbF, ixbF).scale(2 * c)), Matrix.zeros(NbF[k], Nn[k]),
                    sub(Qv, ixn, ixbF).scale(1.0 / c), sub(Qv, ixn, ixn).scale(1.0 / c).add(sub(Rv, ixn, ixn)));
            Matrix Fm = block(
                    sub(Qv, ixbF, ixbF).scale(1.0 / c).add(sub(Sv, ixbF, ixbF).scale(c)).add(sub(Rv, ixbF, ixbF)), sub(Qv, ixbF, ixn).scale(1.0 / c),
                    Matrix.zeros(Nn[k], NbF[k]), Matrix.zeros(Nn[k], Nn[k]));
            Matrix QBDR = QBDFundamentalMatrices.QBDFundamentalMatrices(Bm, Lm, Fm, prec, null, null, null).get("R");
            KF.set(k, sub(QBDR, range(NbF[k]), range(NbF[k])).sub(Matrix.eye(NbF[k])).scale(c));
            Matrix PsiF = sub(QBDR, range(NbF[k]), shiftRange(NbF[k], Nn[k]));
            Matrix clovF = Matrix.zeros(NbF[k], NbF[k] + Nn[k]);
            setCols(clovF, ixbF, Matrix.eye(NbF[k]));
            if (ixn.length > 0) setCols(clovF, ixn, PsiF);
            Matrix cloFk = Matrix.zeros(NbF[k], N);
            setCols(cloFk, ixn0, clovF);
            if (ix0.length > 0) setCols(cloFk, ix0, clovF.mult(Qn0).mult(Q00.scale(-1.0).pinv()));
            cloF.set(k, cloFk);

            // ----- BACKWARD -----
            c1 = ixn.length > 0 ? maxRatio(neg(diagArray(sub(Qv, ixn, ixn))), neg(diagArray(sub(Rv, ixn, ixn)))) : Double.NEGATIVE_INFINITY;
            c2 = cBackward(Rv, Sv, Qv, ixs);
            c = Math.max(Math.max(c1, c2), 1.0);
            int[] ixbB = concat(ixs, ixn);
            NbB[k] = ixbB.length;
            Bm = blkdiag(sub(Sv, ixbB, ixbB).scale(c), sub(Rv, ixp, ixp));
            Lm = block(
                    sub(Rv, ixbB, ixbB).sub(sub(Sv, ixbB, ixbB).scale(2 * c)), Matrix.zeros(NbB[k], Np[k]),
                    sub(Qv, ixp, ixbB).scale(1.0 / c), sub(Qv, ixp, ixp).scale(1.0 / c).sub(sub(Rv, ixp, ixp)));
            Fm = block(
                    sub(Qv, ixbB, ixbB).scale(1.0 / c).add(sub(Sv, ixbB, ixbB).scale(c)).sub(sub(Rv, ixbB, ixbB)), sub(Qv, ixbB, ixp).scale(1.0 / c),
                    Matrix.zeros(Np[k], NbB[k]), Matrix.zeros(Np[k], Np[k]));
            QBDR = QBDFundamentalMatrices.QBDFundamentalMatrices(Bm, Lm, Fm, prec, null, null, null).get("R");
            KB.set(k, sub(QBDR, range(NbB[k]), range(NbB[k])).sub(Matrix.eye(NbB[k])).scale(c));
            Matrix PsiB = sub(QBDR, range(NbB[k]), shiftRange(NbB[k], Np[k]));
            Matrix clovB = Matrix.zeros(NbB[k], NbB[k] + Np[k]);
            setCols(clovB, ixbB, Matrix.eye(NbB[k]));
            if (ixp.length > 0) setCols(clovB, ixp, PsiB);
            Matrix cloBk = Matrix.zeros(NbB[k], N);
            setCols(cloBk, ixn0, clovB);
            if (ix0.length > 0) setCols(cloBk, ix0, clovB.mult(Qn0).mult(Q00.scale(-1.0).pinv()));
            cloB.set(k, cloBk);

            vixp[k] = gather(ixn0, ixp);
            vixn[k] = gather(ixn0, ixn);
            vix0[k] = ix0;
            vixs[k] = gather(ixn0, ixs);
            S.set(k, Sk); // keep the halved variance for the boundary equations
        }

        // ---- boundary equations ----
        int Neqns = (K + 1) * N + sum(Np) + sum(Nn) + 2 * sum(Ns);
        Matrix M = Matrix.zeros(Neqns, Neqns);
        int[] pos = new int[1 + 3 * K];
        pos[0] = 1;
        for (int k = 0; k < K; k++) { pos[1 + 3 * k] = N; pos[2 + 3 * k] = NbF[k]; pos[3 + 3 * k] = NbB[k]; }
        int[] pp = cumsum(pos);

        int i = 0;
        for (int k = 0; k <= K; k++) {
            setBlock(M, pp[i] - 1, N, k * N, Qt.get(k).scale(-1.0));
            if (k > 0) {
                Matrix eKF = KF.get(k - 1).scale(Tarr[k] - Tarr[k - 1]).expm();
                Matrix fwd = eKF.mult(cloF.get(k - 1).mult(R.get(k - 1)).scale(-1.0).add(KF.get(k - 1).mult(cloF.get(k - 1)).mult(S.get(k - 1))));
                setBlock(M, pp[i - 2] - 1, NbF[k - 1], k * N, fwd);
                Matrix bwd = cloB.get(k - 1).mult(R.get(k - 1)).scale(-1.0).sub(KB.get(k - 1).mult(cloB.get(k - 1)).mult(S.get(k - 1)));
                setBlock(M, pp[i - 1] - 1, NbB[k - 1], k * N, bwd);
            }
            if (k < K) {
                Matrix fwd = cloF.get(k).mult(R.get(k)).sub(KF.get(k).mult(cloF.get(k)).mult(S.get(k)));
                setBlock(M, pp[i + 1] - 1, NbF[k], k * N, fwd);
                Matrix eKB = KB.get(k).scale(Tarr[k + 1] - Tarr[k]).expm();
                Matrix bwd = eKB.mult(cloB.get(k).mult(R.get(k)).add(KB.get(k).mult(cloB.get(k)).mult(S.get(k))));
                setBlock(M, pp[i + 2] - 1, NbB[k], k * N, bwd);
            }
            i += 3;
        }

        int col = (K + 1) * N + 1;
        int[] ixAll = range(N);
        i = 0;
        for (int k = 0; k <= K; k++) {
            if (k == 0) {
                int[] ixr0 = intersect(idxWhereEq(boundaryL, 0.0), vixs[0]);
                int Nr0 = ixr0.length;
                Matrix ms0 = Matrix.zeros(N, Np[0] + Nr0);
                int[] sel = concat(vixp[0], ixr0);
                for (int r = 0; r < sel.length; r++) ms0.set(sel[r], r, 1.0);
                setBlockCols(M, pp[i] - 1, N, col - 1, ms0);
                col += Np[0] + Nr0;
                int[] ixa0 = intersect(idxWhereEq(boundaryL, 1.0), vixs[0]);
                int Na0 = ixa0.length;
                if (Na0 > 0) {
                    Matrix pdfF = cloF.get(0);
                    Matrix pdfB = KB.get(0).scale(Tarr[1]).expm().mult(cloB.get(0));
                    setBlockCols(M, pp[i + 1] - 1, NbF[0], col - 1, sub(pdfF, range(NbF[0]), ixa0));
                    setBlockCols(M, pp[i + 2] - 1, NbB[0], col - 1, sub(pdfB, range(NbB[0]), ixa0));
                }
                col += Na0;
            } else if (k == K) {
                int[] ixrB = intersect(idxWhereEq(boundaryU, 0.0), vixs[K - 1]);
                int NrB = ixrB.length;
                Matrix ms0 = Matrix.zeros(N, Nn[K - 1] + NrB);
                int[] sel = concat(vixn[K - 1], ixrB);
                for (int r = 0; r < sel.length; r++) ms0.set(sel[r], r, 1.0);
                setBlockCols(M, pp[i] - 1, N, col - 1, ms0);
                col += Nn[K - 1] + NrB;
                int[] ixaB = intersect(idxWhereEq(boundaryU, 1.0), vixs[K - 1]);
                int NaB = ixaB.length;
                if (NaB > 0) {
                    Matrix pdfF = KF.get(K - 1).scale(Tarr[K] - Tarr[K - 1]).expm().mult(cloF.get(K - 1));
                    Matrix pdfB = cloB.get(K - 1);
                    setBlockCols(M, pp[i - 2] - 1, NbF[K - 1], col - 1, sub(pdfF, range(NbF[K - 1]), ixaB));
                    setBlockCols(M, pp[i - 1] - 1, NbB[K - 1], col - 1, sub(pdfB, range(NbB[K - 1]), ixaB));
                }
                col += NaB;
            } else {
                int[] st0 = complement(union(intersect(vixp[k - 1], vixn[k]), union(vix0[k - 1], vix0[k])), N);
                int N0 = st0.length;
                Matrix ms0 = Matrix.zeros(N, N0);
                for (int r = 0; r < st0.length; r++) ms0.set(st0[r], r, 1.0);
                setBlockCols(M, pp[i] - 1, N, col - 1, ms0);
                col += N0;
                int[] sts = complementSet(union(vixn[k], vixp[k - 1]), union(vixs[k - 1], vixs[k]));
                int Nss = sts.length;
                if (Nss > 0) {
                    Matrix sqrtSk = S.get(k - 1).sqrt(), sqrtSk1 = S.get(k).sqrt();
                    Matrix BelowF = KF.get(k - 1).scale(Tarr[k] - Tarr[k - 1]).expm().mult(cloF.get(k - 1).scale(-1.0)).mult(sqrtSk);
                    Matrix BelowB = cloB.get(k - 1).scale(-1.0).mult(sqrtSk);
                    Matrix AboveF = cloF.get(k).mult(sqrtSk1);
                    Matrix AboveB = KB.get(k).scale(Tarr[k + 1] - Tarr[k]).expm().mult(cloB.get(k)).mult(sqrtSk1);
                    setBlockCols(M, pp[i - 2] - 1, NbF[k - 1], col - 1, sub(BelowF, range(NbF[k - 1]), sts));
                    setBlockCols(M, pp[i - 1] - 1, NbB[k - 1], col - 1, sub(BelowB, range(NbB[k - 1]), sts));
                    setBlockCols(M, pp[i + 1] - 1, NbF[k], col - 1, sub(AboveF, range(NbF[k]), sts));
                    setBlockCols(M, pp[i + 2] - 1, NbB[k], col - 1, sub(AboveB, range(NbB[k]), sts));
                }
                col += Nss;
            }
            i += 3;
        }

        // normalizing condition
        Matrix h = Matrix.ones(N, 1);
        for (int k = 0; k < K; k++) {
            Matrix[] ie = integExp2(KF.get(k), KB.get(k), Tarr[k + 1] - Tarr[k]);
            h = vstack(h, rowSums(ie[0].mult(cloF.get(k))));
            h = vstack(h, rowSums(ie[1].mult(cloB.get(k))));
            h = vstack(h, Matrix.ones(N, 1));
        }
        for (int r = 0; r < Neqns; r++) M.set(r, 0, h.get(r, 0));

        Matrix e1 = Matrix.zeros(Neqns, 1);
        e1.set(0, 0, 1.0);
        Matrix b = M.transpose().leftMatrixDivide(e1); // (Neqns,1)

        List<Matrix> masses = nulls(K + 1), iniF = nulls(K), iniB = nulls(K);
        masses.set(0, sub(b, range(N), new int[]{0}).transpose());
        int ii = 1;
        for (int k = 0; k < K; k++) {
            iniF.set(k, sub(b, shiftRange(pp[ii] - 1, NbF[k]), new int[]{0}).transpose());
            iniB.set(k, sub(b, shiftRange(pp[ii + 1] - 1, NbB[k]), new int[]{0}).transpose());
            masses.set(k + 1, sub(b, shiftRange(pp[ii + 2] - 1, N), new int[]{0}).transpose());
            ii += 3;
        }
        return new LevelDependentFluidSolution(masses, iniF, KF, cloF, iniB, KB, cloB);
    }

    // ---------- small helpers ----------
    private static List<Matrix> nulls(int n) {
        List<Matrix> l = new ArrayList<Matrix>();
        for (int i = 0; i < n; i++) l.add(null);
        return l;
    }
    private static int sum(int[] a) { int s = 0; for (int x : a) s += x; return s; }
    private static int[] cumsum(int[] a) { int[] r = new int[a.length]; int s = 0; for (int i = 0; i < a.length; i++) { s += a[i]; r[i] = s; } return r; }
    private static double[] neg(double[] a) { double[] r = new double[a.length]; for (int i = 0; i < a.length; i++) r[i] = -a[i]; return r; }
    private static double maxRatio(double[] num, double[] den) {
        double m = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < num.length; i++) m = Math.max(m, num[i] / den[i]);
        return m;
    }
    private static double cForward(Matrix Rv, Matrix Sv, Matrix Qv, int[] ixs) {
        double c2 = Double.NEGATIVE_INFINITY;
        for (int idx : ixs) {
            double rv = Rv.get(idx, idx), sv = Sv.get(idx, idx), qv = Qv.get(idx, idx);
            double discr = rv * rv - 2 * (2 * sv) * qv;
            if (discr > 0) c2 = Math.max(c2, (-rv + Math.sqrt(discr)) / (2 * sv));
        }
        return c2;
    }
    private static double cBackward(Matrix Rv, Matrix Sv, Matrix Qv, int[] ixs) {
        double c2 = Double.NEGATIVE_INFINITY;
        for (int idx : ixs) {
            double rv = Rv.get(idx, idx), sv = Sv.get(idx, idx), qv = Qv.get(idx, idx);
            double discr = rv * rv - 2 * (2 * sv) * qv;
            if (discr > 0) c2 = Math.max(c2, (rv + Math.sqrt(discr)) / (2 * sv));
        }
        return c2;
    }
    private static int[] where2(double[] dR, double[] dS, double prec, boolean zero) {
        List<Integer> l = new ArrayList<Integer>();
        for (int i = 0; i < dR.length; i++) if (Math.abs(dR[i]) <= prec && dS[i] <= prec) l.add(i);
        return toArr(l);
    }
    private static int[] whereRS(double[] dRv, double[] dSv, double prec, int sign) {
        List<Integer> l = new ArrayList<Integer>();
        for (int i = 0; i < dRv.length; i++)
            if (dSv[i] <= prec && ((sign > 0 && dRv[i] > prec) || (sign < 0 && dRv[i] < -prec))) l.add(i);
        return toArr(l);
    }
    private static int[] whereS(double[] dSv, double prec) {
        List<Integer> l = new ArrayList<Integer>();
        for (int i = 0; i < dSv.length; i++) if (dSv[i] > prec) l.add(i);
        return toArr(l);
    }
    private static int[] complement(int[] subset, int N) {
        boolean[] in = new boolean[N];
        for (int x : subset) in[x] = true;
        List<Integer> l = new ArrayList<Integer>();
        for (int i = 0; i < N; i++) if (!in[i]) l.add(i);
        return toArr(l);
    }
    private static int[] complementSet(int[] subset, int[] universe) {
        boolean[] in = new boolean[0];
        java.util.Set<Integer> rm = new java.util.HashSet<Integer>();
        for (int x : subset) rm.add(x);
        List<Integer> l = new ArrayList<Integer>();
        for (int x : universe) if (!rm.contains(x)) l.add(x);
        return toArr(l);
    }
    private static int[] gather(int[] base, int[] idx) {
        int[] r = new int[idx.length];
        for (int i = 0; i < idx.length; i++) r[i] = base[idx[i]];
        return r;
    }
    private static int[] intersect(int[] a, int[] b) {
        java.util.Set<Integer> sb = new java.util.TreeSet<Integer>();
        for (int x : b) sb.add(x);
        java.util.Set<Integer> sa = new java.util.TreeSet<Integer>();
        for (int x : a) sa.add(x);
        List<Integer> l = new ArrayList<Integer>();
        for (int x : sa) if (sb.contains(x)) l.add(x);
        return toArr(l);
    }
    private static int[] union(int[] a, int[] b) {
        java.util.TreeSet<Integer> s = new java.util.TreeSet<Integer>();
        for (int x : a) s.add(x);
        for (int x : b) s.add(x);
        List<Integer> l = new ArrayList<Integer>(s);
        return toArr(l);
    }
    private static int[] idxWhereEq(double[] a, double v) {
        List<Integer> l = new ArrayList<Integer>();
        for (int i = 0; i < a.length; i++) if (a[i] == v) l.add(i);
        return toArr(l);
    }
    private static int[] toArr(List<Integer> l) {
        int[] r = new int[l.size()];
        for (int i = 0; i < r.length; i++) r[i] = l.get(i);
        return r;
    }
    private static void setBlock(Matrix M, int row0, int nrows, int col0, Matrix src) {
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < src.getNumCols(); c++)
                M.set(row0 + r, col0 + c, src.get(r, c));
    }
    private static void setBlockCols(Matrix M, int row0, int nrows, int col0, Matrix src) {
        for (int r = 0; r < nrows; r++)
            for (int c = 0; c < src.getNumCols(); c++)
                M.set(row0 + r, col0 + c, src.get(r, c));
    }
}
