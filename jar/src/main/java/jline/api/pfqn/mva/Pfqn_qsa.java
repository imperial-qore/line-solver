/**
 * @file Queue-Shift Approximation (QSA) for closed product-form networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

/**
 * Schweitzer, Serazzi and Broglia, "A Queue-Shift Approximation Technique for
 * Product-Form Queueing Networks", Tools'98, LNCS 1469, pp. 267-279.
 *
 * QSA approximates the arrival-instant queue lengths through the absolute shift
 * of the AGGREGATE queue length, Y_ri(K) = 1 + Q_i(K - e_r) - Q_i(K), in place
 * of the fractional deviations of Linearizer, so the unknowns are one per
 * station rather than one per station-class. The core equation (13a) is imposed
 * at K, at every K - e_s and, in the three-level variant of eq. (16), at every
 * K - e_s - e_t through the affine extrapolation of eq. (15).
 *
 * Port of matlab/src/api/pfqn/pfqn_qsa.m. The quintuple (16) is solved as ONE
 * system by damped Newton, as Sect. 4 of the paper prescribes: the decomposed
 * successive substitution that works for Linearizer drifts to the degenerate
 * root in which the bottleneck absorbs the whole population.
 */
public final class Pfqn_qsa {
    private Pfqn_qsa() {}

    public static Ret.pfqnAMVA pfqn_qsa(Matrix L, Matrix N) {
        return pfqn_qsa(L, N, new Matrix(1, L.getNumCols()), null, 1e-10, 100, 3, null);
    }

    public static Ret.pfqnAMVA pfqn_qsa(Matrix L, Matrix N, Matrix Z) {
        return pfqn_qsa(L, N, Z, null, 1e-10, 100, 3, null);
    }

    public static Ret.pfqnAMVA pfqn_qsa(Matrix L, Matrix N, Matrix Z, SchedStrategy[] type,
                                        double tol, int maxiter) {
        return pfqn_qsa(L, N, Z, type, tol, maxiter, 3, null);
    }

    /**
     * @param L       service demand matrix (M x R)
     * @param N       population vector (1 x R)
     * @param Z       think time vector (1 x R)
     * @param type    scheduling strategy per station; SchedStrategy.INF marks a delay centre
     * @param tol     residual tolerance of the Newton iteration
     * @param maxiter maximum Newton iterations
     * @param levels  2 for the two-level QSA of eq. (14), 3 for eq. (16)
     * @param QN0     warm start for the Bard-Schweitzer initialization (M x R), may be null
     */
    public static Ret.pfqnAMVA pfqn_qsa(Matrix L, Matrix N, Matrix Z, SchedStrategy[] type,
                                        double tol, int maxiter, int levels, Matrix QN0) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (Z == null || Z.isEmpty()) {
            Z = new Matrix(1, R);
        }
        boolean[] isQC = new boolean[M];
        for (int i = 0; i < M; i++) {
            isQC[i] = type == null || i >= type.length || type[i] != SchedStrategy.INF;
        }

        Matrix Q = new Matrix(M, R);
        Matrix U = new Matrix(M, R);
        Matrix W = new Matrix(M, R);
        Matrix C = new Matrix(1, R);
        Matrix X = new Matrix(1, R);
        Matrix T = new Matrix(M, R);

        boolean emptyDemands = true;
        for (int j = 0; j < R && emptyDemands; j++) {
            for (int i = 0; i < M; i++) {
                if (L.get(i, j) != 0) {
                    emptyDemands = false;
                    break;
                }
            }
        }
        boolean emptyPop = true;
        for (int r = 0; r < R; r++) {
            if (N.get(0, r) > 0) {
                emptyPop = false;
            }
        }
        if (M == 0 || emptyDemands || emptyPop) {
            for (int r = 0; r < R; r++) {
                if (N.get(0, r) > 0 && Z.get(0, r) > 0) {
                    X.set(0, r, N.get(0, r) / Z.get(0, r));
                }
                for (int i = 0; i < M; i++) {
                    U.set(i, r, X.get(0, r) * L.get(i, r));
                }
            }
            return new Ret.pfqnAMVA(Q, U, W, T, C, X, 0);
        }

        // Populations touched by (16): K, every K - e_s, every K - e_s - e_t.
        List<double[]> popList = new ArrayList<double[]>();
        double[] base = new double[R];
        for (int r = 0; r < R; r++) {
            base[r] = N.get(0, r);
        }
        popList.add(base);
        int[] sIdx = new int[R];
        int[][] pIdx = new int[R][R];
        for (int s = 0; s < R; s++) {
            sIdx[s] = -1;
            for (int t = 0; t < R; t++) {
                pIdx[s][t] = -1;
            }
        }
        for (int s = 0; s < R; s++) {
            double[] n = base.clone();
            n[s] -= 1.0;
            if (nonNegative(n)) {
                popList.add(n);
                sIdx[s] = popList.size() - 1;
            }
        }
        if (levels >= 3) {
            for (int s = 0; s < R; s++) {
                for (int t = s; t < R; t++) {
                    double[] n = base.clone();
                    n[s] -= 1.0;
                    n[t] -= 1.0;
                    if (nonNegative(n)) {
                        popList.add(n);
                        pIdx[s][t] = popList.size() - 1;
                        pIdx[t][s] = pIdx[s][t];
                    }
                }
            }
        }
        double[][] pops = popList.toArray(new double[popList.size()][]);
        int nP = pops.length;

        int mq = 0;
        for (int i = 0; i < M; i++) {
            if (isQC[i]) {
                mq++;
            }
        }
        int[] qc = new int[mq];
        int k = 0;
        for (int i = 0; i < M; i++) {
            if (isQC[i]) {
                qc[k++] = i;
            }
        }
        double[] Ldc = new double[R];
        for (int r = 0; r < R; r++) {
            for (int i = 0; i < M; i++) {
                if (!isQC[i]) {
                    Ldc[r] += L.get(i, r);
                }
            }
        }

        // Bard-Schweitzer at every population supplies the Newton starting point.
        double[][] q = new double[M][nP];
        for (int p = 0; p < nP; p++) {
            double[] qp = aggbs(L, pops[p], Z, isQC, QN0);
            for (int i = 0; i < M; i++) {
                q[i][p] = qp[i];
            }
        }

        int nUnk = mq * nP;
        double[] x = new double[nUnk];
        for (int a = 0; a < mq; a++) {
            for (int p = 0; p < nP; p++) {
                x[a * nP + p] = q[qc[a]][p];
            }
        }
        boolean[] admFlag = new boolean[1];
        double[] F = resid(x, L, Z, pops, sIdx, pIdx, qc, Ldc, levels, admFlag);
        double fnrm = norm(F);
        int totiter = 0;
        Matrix Jm = new Matrix(nUnk, nUnk);
        Matrix rhs = new Matrix(nUnk, 1);
        Matrix step = new Matrix(nUnk, 1);
        for (int it = 1; it <= maxiter; it++) {
            if (fnrm < tol) {
                break;
            }
            totiter = it;
            for (int col = 0; col < nUnk; col++) {
                double h = 1e-7 * Math.max(1.0, Math.abs(x[col]));
                double[] xp = x.clone();
                xp[col] += h;
                double[] Fp = resid(xp, L, Z, pops, sIdx, pIdx, qc, Ldc, levels, new boolean[1]);
                for (int row = 0; row < nUnk; row++) {
                    Jm.set(row, col, (Fp[row] - F[row]) / h);
                }
            }
            for (int row = 0; row < nUnk; row++) {
                rhs.set(row, 0, -F[row]);
            }
            if (!Matrix.solveSafe(Jm, rhs, step)) {
                break;
            }
            boolean accepted = false;
            double lambda = 1.0;
            for (int ls = 0; ls < 40; ls++) {
                double[] xn = new double[nUnk];
                boolean finite = true;
                for (int row = 0; row < nUnk; row++) {
                    xn[row] = x[row] + lambda * step.get(row, 0);
                    if (Double.isNaN(xn[row]) || Double.isInfinite(xn[row])) {
                        finite = false;
                    }
                }
                if (finite) {
                    boolean[] admn = new boolean[1];
                    double[] Fn = resid(xn, L, Z, pops, sIdx, pIdx, qc, Ldc, levels, admn);
                    double nn = norm(Fn);
                    if (admn[0] && nn < fnrm) {
                        x = xn;
                        F = Fn;
                        fnrm = nn;
                        accepted = true;
                        break;
                    }
                }
                lambda /= 2.0;
            }
            if (!accepted) {
                break;
            }
        }

        // Disaggregate (13) at K into the per-class measures
        for (int a = 0; a < mq; a++) {
            for (int p = 0; p < nP; p++) {
                q[qc[a]][p] = x[a * nP + p];
            }
        }
        double[][] Y0 = shift(q, M, R, pops, sIdx, pIdx, -1, -1, levels);
        for (int r = 0; r < R; r++) {
            if (N.get(0, r) < 1) {
                continue;
            }
            double sumW = 0.0;
            for (int i = 0; i < M; i++) {
                double w = isQC[i] ? L.get(i, r) * (q[i][0] + Y0[i][r]) : L.get(i, r);
                W.set(i, r, w);
                sumW += w;
            }
            double xr = N.get(0, r) / (Z.get(0, r) + sumW);
            X.set(0, r, xr);
            for (int i = 0; i < M; i++) {
                Q.set(i, r, xr * W.get(i, r));
                U.set(i, r, xr * L.get(i, r));
                T.set(i, r, xr);
            }
            C.set(0, r, N.get(0, r) / xr - Z.get(0, r));
        }
        return new Ret.pfqnAMVA(Q, U, W, T, C, X, totiter);
    }

    /**
     * Residual of (13) imposed simultaneously at every population of (16). The
     * adm flag carries the side conditions of Remark 2 (non-negative queue
     * lengths, positive cycle times).
     */
    private static double[] resid(double[] x, Matrix L, Matrix Z, double[][] pops,
                                  int[] sIdx, int[][] pIdx, int[] qc, double[] Ldc,
                                  int levels, boolean[] adm) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        int nP = pops.length;
        int mq = qc.length;
        adm[0] = true;
        double[][] q = new double[M][nP];
        for (int a = 0; a < mq; a++) {
            for (int p = 0; p < nP; p++) {
                q[qc[a]][p] = x[a * nP + p];
                if (x[a * nP + p] < 0) {
                    adm[0] = false;
                }
            }
        }
        double[] F = new double[mq * nP];
        for (int p = 0; p < nP; p++) {
            double[] np = pops[p];
            int[] st = which(p, sIdx, pIdx);
            double[][] Y = shift(q, M, R, pops, sIdx, pIdx, st[0], st[1], levels);
            double[] acc = new double[mq];
            for (int r = 0; r < R; r++) {
                if (np[r] < 1) {
                    continue;
                }
                double c = Z.get(0, r) + Ldc[r];
                for (int a = 0; a < mq; a++) {
                    c += L.get(qc[a], r) * (q[qc[a]][p] + Y[qc[a]][r]);
                }
                if (!(c > 0) || Double.isNaN(c) || Double.isInfinite(c)) {
                    adm[0] = false;
                    c = Math.ulp(1.0);
                }
                double xr = np[r] / c;
                for (int a = 0; a < mq; a++) {
                    acc[a] += xr * L.get(qc[a], r) * (q[qc[a]][p] + Y[qc[a]][r]);
                }
            }
            for (int a = 0; a < mq; a++) {
                F[a * nP + p] = q[qc[a]][p] - acc[a];
            }
        }
        return F;
    }

    /** Shift matrix of (16d)-(16e), or (15) when both s and t are set. */
    private static double[][] shift(double[][] q, int M, int R, double[][] pops,
                                    int[] sIdx, int[][] pIdx, int s, int t, int levels) {
        double[][] Y = new double[M][R];
        if (s < 0) {
            for (int r = 0; r < R; r++) {
                if (sIdx[r] >= 0) {
                    for (int i = 0; i < M; i++) {
                        Y[i][r] = 1.0 + q[i][sIdx[r]] - q[i][0];
                    }
                }
            }
        } else if (t < 0) {
            if (levels < 3) {
                return shift(q, M, R, pops, sIdx, pIdx, -1, -1, levels);   // (14)
            }
            for (int r = 0; r < R; r++) {
                if (pIdx[s][r] >= 0 && pops[sIdx[s]][r] >= 1) {
                    for (int i = 0; i < M; i++) {
                        Y[i][r] = 1.0 + q[i][pIdx[s][r]] - q[i][sIdx[s]];
                    }
                }
            }
        } else {
            double[][] Ys = shift(q, M, R, pops, sIdx, pIdx, s, -1, levels);
            double[][] Yt = shift(q, M, R, pops, sIdx, pIdx, t, -1, levels);
            double[][] Y0 = shift(q, M, R, pops, sIdx, pIdx, -1, -1, levels);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    Y[i][r] = Ys[i][r] + Yt[i][r] - Y0[i][r];
                }
            }
        }
        return Y;
    }

    /** Decode a population index into the removed classes. */
    private static int[] which(int p, int[] sIdx, int[][] pIdx) {
        if (p == 0) {
            return new int[]{-1, -1};
        }
        for (int s = 0; s < sIdx.length; s++) {
            if (sIdx[s] == p) {
                return new int[]{s, -1};
            }
        }
        for (int s = 0; s < pIdx.length; s++) {
            for (int t = 0; t < pIdx.length; t++) {
                if (pIdx[s][t] == p) {
                    return new int[]{s, t};
                }
            }
        }
        return new int[]{-1, -1};
    }

    /**
     * Aggregate Bard-Schweitzer queue lengths at population n, with the
     * delay-centre demands folded into the think time.
     */
    private static double[] aggbs(Matrix L, double[] nIn, Matrix Z, boolean[] isQC, Matrix QN0) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double[] q = new double[M];
        Matrix n = new Matrix(1, R);
        boolean empty = true;
        for (int r = 0; r < R; r++) {
            double v = Math.max(nIn[r], 0.0);
            n.set(0, r, v);
            if (v > 0) {
                empty = false;
            }
        }
        if (empty) {
            return q;
        }
        Matrix Zeff = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double zd = Z.get(0, r);
            for (int i = 0; i < M; i++) {
                if (!isQC[i]) {
                    zd += L.get(i, r);
                }
            }
            Zeff.set(0, r, zd);
        }
        int mq = 0;
        for (int i = 0; i < M; i++) {
            if (isQC[i]) {
                mq++;
            }
        }
        double[] X = new double[R];
        if (mq > 0) {
            Matrix Lq = new Matrix(mq, R);
            Matrix Q0q = QN0 == null || QN0.isEmpty() ? null : new Matrix(mq, R);
            int a = 0;
            for (int i = 0; i < M; i++) {
                if (isQC[i]) {
                    for (int r = 0; r < R; r++) {
                        Lq.set(a, r, L.get(i, r));
                        if (Q0q != null) {
                            Q0q.set(a, r, QN0.get(i, r));
                        }
                    }
                    a++;
                }
            }
            Ret.pfqnAMVA bs = Pfqn_bs.pfqn_bs(Lq, n, Zeff, 1e-6, 1000, Q0q);
            a = 0;
            for (int i = 0; i < M; i++) {
                if (isQC[i]) {
                    double s = 0.0;
                    for (int r = 0; r < R; r++) {
                        s += bs.Q.get(a, r);
                    }
                    q[i] = s;
                    a++;
                }
            }
            for (int r = 0; r < R; r++) {
                X[r] = bs.X.get(0, r);
            }
        } else {
            for (int r = 0; r < R; r++) {
                if (n.get(0, r) >= 1 && Zeff.get(0, r) > 0) {
                    X[r] = n.get(0, r) / Zeff.get(0, r);
                }
            }
        }
        for (int i = 0; i < M; i++) {
            if (!isQC[i]) {
                double s = 0.0;
                for (int r = 0; r < R; r++) {
                    s += X[r] * L.get(i, r);
                }
                q[i] = s;
            }
        }
        return q;
    }

    private static boolean nonNegative(double[] n) {
        for (int i = 0; i < n.length; i++) {
            if (n[i] < 0) {
                return false;
            }
        }
        return true;
    }

    private static double norm(double[] v) {
        double s = 0.0;
        for (int i = 0; i < v.length; i++) {
            s += v[i] * v[i];
        }
        return Math.sqrt(s);
    }
}
