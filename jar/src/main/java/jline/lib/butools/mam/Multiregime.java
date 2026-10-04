/*
 * Ported from BUTools-family fluid tools (Kankaya & Akar).
 *
 * Multi-regime feedback fluid queue. The generator and drift rates are regime
 * dependent, with separate boundary (feedback) generators/rates at each
 * threshold.
 */
package jline.lib.butools.mam;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

import static jline.lib.butools.mam.FluidTools.*;

public final class Multiregime {
    private Multiregime() {}

    /**
     * Solves a multi-regime feedback Markovian fluid queue. Q, Qt are lists of
     * generators (per regime k=1..K and per boundary k=0..K); R, Rt are lists of
     * diagonal drift-rate VECTORS; T is the vector of regime thresholds.
     *
     * @return {pdf, pdfd, cdf, cdfm}, each (points x N).
     */
    public static Matrix[] multiregime(List<Matrix> Q, List<double[]> R, List<Matrix> Qt, List<double[]> Rt,
                                       double[] T, double[] pdfpoints, double[] cdfpoints) {
        int K = R.size();
        int N = R.get(0).length;
        double[] Tarr = new double[K + 1];
        for (int k = 0; k < K; k++) Tarr[k + 1] = T[k];

        if (Q.size() == 1) {
            Matrix q0 = Q.get(0);
            Q = new ArrayList<Matrix>();
            for (int k = 0; k < K; k++) Q.add(q0);
        }
        if (Qt == null || Qt.isEmpty()) {
            Qt = new ArrayList<Matrix>();
            Qt.add(Q.get(0));
            for (int k = 0; k < K; k++) Qt.add(Q.get(k));
        } else if (Qt.size() == 1) {
            Matrix q0 = Qt.get(0);
            Qt = new ArrayList<Matrix>();
            for (int k = 0; k < K + 1; k++) Qt.add(q0);
        }
        if (Rt == null || Rt.isEmpty()) {
            Rt = new ArrayList<double[]>();
            Rt.add(R.get(0));
            for (int k = 0; k < K; k++) Rt.add(R.get(k));
        }

        Matrix[] An = new Matrix[K], Ap = new Matrix[K];
        Matrix[] L0 = new Matrix[K], Ln = new Matrix[K], Lp = new Matrix[K];
        Matrix[] M0 = new Matrix[K], MT = new Matrix[K], Mi = new Matrix[K];
        int[] zeig = new int[K], neig = new int[K], peig = new int[K], Nnz = new int[K];

        for (int k = 0; k < K; k++) {
            int[] zix = idxEq(R.get(k), 0.0);
            int[] nzix = idxNe(R.get(k), 0.0);
            int Nn = nzix.length;
            Matrix Qk = Q.get(k);
            Matrix Qnn = sub(Qk, nzix, nzix);
            Matrix Qnk;
            if (zix.length > 0) {
                Qnk = Qnn.add(sub(Qk, nzix, zix).mult(sub(Qk, zix, zix).scale(-1.0).pinv()).mult(sub(Qk, zix, nzix)));
            } else {
                Qnk = Qnn;
            }
            double[] invR = new double[Nn];
            for (int i = 0; i < Nn; i++) invR[i] = 1.0 / R.get(k)[nzix[i]];
            Matrix A = Qnk.mult(diagFrom(invR));

            Object[] os = orderedSchur(A, 1e-10);
            Matrix Z = (Matrix) os[0], D = (Matrix) os[1];
            int[] cnt = (int[]) os[2];
            int ze = cnt[0], ne = cnt[1], pe = cnt[2];
            zeig[k] = ze; neig[k] = ne; peig[k] = pe;

            // decoupling Sylvester solves (D is quasi-triangular from ordered Schur):
            // X1 solves 0*X - X*D22 = D12 ; X2 solves Dnn*X - X*Dpp = Dnp
            Matrix X1 = ze > 0
                    ? triSylvester(Matrix.zeros(ze, ze), sub(D, shiftRange(ze, Nn - ze), shiftRange(ze, Nn - ze)),
                    sub(D, range(ze), shiftRange(ze, Nn - ze)))
                    : Matrix.zeros(0, Nn - ze);
            Matrix X2 = (ne > 0 && pe > 0)
                    ? triSylvester(sub(D, shiftRange(ze, ne), shiftRange(ze, ne)), sub(D, shiftRange(ze + ne, pe), shiftRange(ze + ne, pe)),
                    sub(D, shiftRange(ze, ne), shiftRange(ze + ne, pe)))
                    : Matrix.zeros(ne, pe);

            Matrix blk1 = Matrix.eye(Nn);
            if (ze > 0) setSub(blk1, range(ze), shiftRange(ze, Nn - ze), X1.scale(-1.0));
            Matrix blk2 = Matrix.eye(Nn);
            if (ne > 0 && pe > 0) setSub(blk2, shiftRange(ze, ne), shiftRange(ze + ne, pe), X2.scale(-1.0));
            Matrix Y = Z.mult(blk1).mult(blk2);
            Matrix iY = Y.inv();
            Matrix iY0 = sub(iY, range(ze), range(N));
            Matrix iYn = sub(iY, shiftRange(ze, ne), range(N));
            Matrix iYp = sub(iY, shiftRange(ze + ne, pe), range(N));
            Matrix At = iY.mult(A).mult(Y);
            An[k] = sub(At, shiftRange(ze, ne), shiftRange(ze, ne));
            Ap[k] = sub(At, shiftRange(ze + ne, pe), shiftRange(ze + ne, pe));

            Matrix corr = zix.length > 0 ? sub(Qk, nzix, zix).mult(sub(Qk, zix, zix).scale(-1.0).pinv()) : null;
            L0[k] = buildL(iY0, nzix, zix, corr, N);
            Ln[k] = buildL(iYn, nzix, zix, corr, N);
            Lp[k] = buildL(iYp, nzix, zix, corr, N);

            double Tk = Tarr[k + 1] - Tarr[k];
            Matrix eAnT = ne > 0 ? An[k].scale(Tk).expm() : Matrix.zeros(0, 0);
            Matrix eApT = pe > 0 ? Ap[k].scale(-Tk).expm() : Matrix.zeros(0, 0);
            M0[k] = vstack(vstack(L0[k], Ln[k]), pe > 0 ? eApT.mult(Lp[k]) : Matrix.zeros(0, N));
            MT[k] = vstack(vstack(L0[k], ne > 0 ? eAnT.mult(Ln[k]) : Matrix.zeros(0, N)), Lp[k]);
            Matrix MiN = ne > 0 ? An[k].scale(-1.0).pinv().mult(Matrix.eye(ne).sub(eAnT)).mult(Ln[k]) : Matrix.zeros(0, N);
            Matrix MiP = pe > 0 ? Ap[k].pinv().mult(Matrix.eye(pe).sub(eApT)).mult(Lp[k]) : Matrix.zeros(0, N);
            Mi[k] = vstack(vstack(L0[k].scale(Tk), MiN), MiP);
            Nnz[k] = Nn;
        }

        int d = (K + 1) * N;
        int Neqns = d;
        for (int k = 0; k < K; k++) Neqns += Nnz[k];
        Matrix M = Matrix.zeros(Neqns, Neqns);
        int[] dstart = new int[K + 1];
        dstart[0] = d;
        for (int k = 1; k <= K; k++) dstart[k] = dstart[k - 1] + Nnz[k - 1];

        // eq. (12)
        int p = 0;
        setBlk(M, 0, N, p, Qt.get(0).scale(-1.0));
        setBlk(M, dstart[0], Nnz[0], p, M0[0].mult(diagFrom(R.get(0))));
        p += N;
        // eq. (13)
        for (int k = 1; k < K; k++) {
            setBlk(M, k * N, N, p, Qt.get(k).scale(-1.0));
            setBlk(M, dstart[k], Nnz[k], p, M0[k].mult(diagFrom(R.get(k))));
            setBlk(M, dstart[k - 1], Nnz[k - 1], p, MT[k - 1].mult(diagFrom(R.get(k - 1))).scale(-1.0));
            p += N;
        }
        // eq. (16)
        setBlk(M, K * N, N, p, Qt.get(K).scale(-1.0));
        setBlk(M, dstart[K - 1], Nnz[K - 1], p, MT[K - 1].mult(diagFrom(R.get(K - 1))).scale(-1.0));
        p += N;
        // eq. (8)
        for (int m = 0; m < N; m++) if (R.get(0)[m] > 0) { zeroCol(M, p); M.set(m, p, 1.0); p++; }
        // eq. (11)
        for (int m = 0; m < N; m++) if (R.get(K - 1)[m] < 0) { zeroCol(M, p); M.set(K * N + m, p, 1.0); p++; }
        // eq. (9)
        for (int k = 1; k < K; k++) for (int m = 0; m < N; m++)
            if ((R.get(k - 1)[m] > 0 && R.get(k)[m] > 0) || (R.get(k - 1)[m] < 0 && R.get(k)[m] < 0)) { zeroCol(M, p); M.set(k * N + m, p, 1.0); p++; }
        // eq. (10)
        for (int k = 1; k < K; k++) for (int m = 0; m < N; m++)
            if (R.get(k - 1)[m] < 0 && R.get(k)[m] > 0 && Rt.get(k)[m] != 0) { zeroCol(M, p); M.set(k * N + m, p, 1.0); p++; }
        // eq. (14)
        for (int k = 1; k < K; k++) for (int m = 0; m < N; m++)
            if (R.get(k - 1)[m] < 0 && Rt.get(k)[m] >= 0) { zeroCol(M, p); setCol(M, dstart[k - 1], Nnz[k - 1], p, colOf(MT[k - 1], m)); p++; }
        // eq. (15)
        for (int k = 1; k < K; k++) for (int m = 0; m < N; m++)
            if (R.get(k)[m] > 0 && Rt.get(k)[m] <= 0) { zeroCol(M, p); setCol(M, dstart[k], Nnz[k], p, colOf(M0[k], m)); p++; }

        // normalization (column 0)
        for (int r = 0; r < d; r++) M.set(r, 0, 1.0);
        for (int k = 0; k < K; k++) {
            Matrix rs = rowSums(Mi[k]);
            for (int r = 0; r < Nnz[k]; r++) M.set(dstart[k] + r, 0, rs.get(r, 0));
        }
        Matrix rhs = Matrix.zeros(Neqns, 1);
        rhs.set(0, 0, 1.0);
        Matrix sol = linsolve(M.transpose(), rhs); // partial-pivot LU (Neqns,1)

        Matrix[] masses = new Matrix[K + 1];
        Matrix[] a0 = new Matrix[K], an = new Matrix[K], ap = new Matrix[K];
        masses[0] = rowFrom(sol, 0, N);
        for (int k = 0; k < K; k++) {
            masses[k + 1] = rowFrom(sol, (k + 1) * N, N);
            int off = dstart[k];
            a0[k] = rowFrom(sol, off, zeig[k]);
            an[k] = rowFrom(sol, off + zeig[k], neig[k]);
            ap[k] = rowFrom(sol, off + zeig[k] + neig[k], peig[k]);
        }

        Matrix pdf = evalPdf(pdfpoints, K, N, Tarr, a0, an, ap, L0, Ln, Lp, An, Ap, neig, peig, false);
        Matrix pdfd = evalPdf(pdfpoints, K, N, Tarr, a0, an, ap, L0, Ln, Lp, An, Ap, neig, peig, true);
        Matrix[] cc = evalCdf(cdfpoints, K, N, Tarr, a0, an, ap, L0, Ln, Lp, An, Ap, Mi, masses, neig, peig);
        return new Matrix[]{pdf, pdfd, cc[0], cc[1]};
    }

    // ---------------- helpers ----------------

    private static Matrix buildL(Matrix iYx, int[] nzix, int[] zix, Matrix corr, int N) {
        Matrix L = Matrix.zeros(iYx.getNumRows(), N);
        setCols(L, nzix, iYx);
        if (zix.length > 0) setCols(L, zix, iYx.mult(corr));
        return L;
    }

    private static int regime(double p, int K, double[] Tarr) {
        int k = 0;
        while (k < K - 1 && p >= Tarr[k + 1]) k++;
        return k;
    }

    private static Matrix evalPdf(double[] pts, int K, int N, double[] Tarr, Matrix[] a0, Matrix[] an, Matrix[] ap,
                                  Matrix[] L0, Matrix[] Ln, Matrix[] Lp, Matrix[] An, Matrix[] Ap,
                                  int[] neig, int[] peig, boolean deriv) {
        Matrix res = Matrix.zeros(pts.length, N);
        for (int pi = 0; pi < pts.length; pi++) {
            double p = pts[pi];
            int k = regime(p, K, Tarr);
            Matrix pr = Matrix.zeros(1, N);
            if (!deriv) pr = pr.add(a0[k].mult(L0[k]));
            if (neig[k] > 0) {
                Matrix e = An[k].scale(p - Tarr[k]).expm();
                pr = pr.add((deriv ? an[k].mult(An[k]) : an[k]).mult(e).mult(Ln[k]));
            }
            if (peig[k] > 0) {
                Matrix e = Ap[k].scale(-(Tarr[k + 1] - p)).expm();
                pr = pr.add((deriv ? ap[k].mult(Ap[k]) : ap[k]).mult(e).mult(Lp[k]));
            }
            for (int j = 0; j < N; j++) res.set(pi, j, pr.get(0, j));
        }
        return res;
    }

    private static Matrix[] evalCdf(double[] pts, int K, int N, double[] Tarr, Matrix[] a0, Matrix[] an, Matrix[] ap,
                                    Matrix[] L0, Matrix[] Ln, Matrix[] Lp, Matrix[] An, Matrix[] Ap,
                                    Matrix[] Mi, Matrix[] masses, int[] neig, int[] peig) {
        Matrix cdf = Matrix.zeros(pts.length, N), cdfm = Matrix.zeros(pts.length, N);
        for (int pi = 0; pi < pts.length; pi++) {
            double c = pts[pi];
            Matrix cres = Matrix.zeros(1, N), cresm = Matrix.zeros(1, N);
            int k = 0;
            while (k < K && c >= Tarr[k]) {
                if (k > 0) {
                    Matrix avec = hstack(hstack(a0[k - 1], an[k - 1]), ap[k - 1]);
                    cres = cres.add(avec.mult(Mi[k - 1]));
                    cresm = cresm.add(avec.mult(Mi[k - 1]));
                }
                cresm = cresm.add(masses[k]);
                if (c > Tarr[k]) cres = cres.add(masses[k]);
                k++;
            }
            if (k == K && c == Tarr[k]) cresm = cresm.add(masses[k]);
            int kk = k - 1;
            double crem = c - Tarr[kk];
            double Tk = Tarr[kk + 1] - Tarr[kk];
            Matrix val = a0[kk].mult(L0[kk]).scale(crem);
            if (neig[kk] > 0)
                val = val.add(an[kk].mult(An[kk].scale(-1.0).pinv()).mult(Matrix.eye(neig[kk]).sub(An[kk].scale(crem).expm())).mult(Ln[kk]));
            if (peig[kk] > 0)
                val = val.add(ap[kk].mult(Ap[kk].scale(-1.0).pinv()).mult(Ap[kk].scale(-Tk).expm().sub(Ap[kk].scale(-(Tk - crem)).expm())).mult(Lp[kk]));
            cres = cres.add(val);
            cresm = cresm.add(val);
            for (int j = 0; j < N; j++) { cdf.set(pi, j, cres.get(0, j)); cdfm.set(pi, j, cresm.get(0, j)); }
        }
        return new Matrix[]{cdf, cdfm};
    }

    private static Matrix rowFrom(Matrix col, int start, int len) {
        Matrix r = Matrix.zeros(1, len);
        for (int i = 0; i < len; i++) r.set(0, i, col.get(start + i, 0));
        return r;
    }
    private static Matrix colOf(Matrix M, int j) {
        Matrix c = Matrix.zeros(M.getNumRows(), 1);
        for (int i = 0; i < M.getNumRows(); i++) c.set(i, 0, M.get(i, j));
        return c;
    }
    private static void setBlk(Matrix M, int row0, int nrows, int col0, Matrix src) {
        for (int r = 0; r < nrows; r++) for (int c = 0; c < src.getNumCols(); c++) M.set(row0 + r, col0 + c, src.get(r, c));
    }
    private static void setCol(Matrix M, int row0, int nrows, int col, Matrix src) {
        for (int r = 0; r < nrows; r++) M.set(row0 + r, col, src.get(r, 0));
    }
    private static void zeroCol(Matrix M, int col) {
        for (int r = 0; r < M.getNumRows(); r++) M.set(r, col, 0.0);
    }
    private static int[] idxEq(double[] a, double v) {
        List<Integer> l = new ArrayList<Integer>();
        for (int i = 0; i < a.length; i++) if (a[i] == v) l.add(i);
        int[] r = new int[l.size()]; for (int i = 0; i < r.length; i++) r[i] = l.get(i); return r;
    }
    private static int[] idxNe(double[] a, double v) {
        List<Integer> l = new ArrayList<Integer>();
        for (int i = 0; i < a.length; i++) if (a[i] != v) l.add(i);
        int[] r = new int[l.size()]; for (int i = 0; i < r.length; i++) r[i] = l.get(i); return r;
    }
}
