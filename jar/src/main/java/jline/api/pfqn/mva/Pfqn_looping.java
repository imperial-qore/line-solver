/**
 * @file Eager Looping bounds for closed multiclass product-form networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

/**
 * Eager Looping approximate MVA bounds.
 *
 * <p>D. L. Eager, "Bounding Algorithms for Queueing Network Models of Computer
 * Systems", Ph.D. thesis, Tech. Rept. CSRG-156, University of Toronto, 1984.
 * Looping supplies the initial pessimistic and optimistic estimates that the
 * multiple-class performance bound hierarchy starts from, so it carries a pair
 * of bounds rather than a single fixed point.</p>
 *
 * <p>It is built on the convolution identity of Zahorjan (1980)</p>
 *
 * <pre>Q_jk(N - 1_c) = [X_j^{+k}(N - 1_c) / X_j(N)] Q_jk(N),</pre>
 *
 * <p>with X_j^{+k}(N - 1_c) estimated from the level-0 multiple-class PBH upper
 * bound B_j and X_j(N) from the optimistic response time R_j^(opt). A HEAP H_j
 * is the class-j congestion that the current queue-length lower bounds have not
 * yet accounted for; it is charged back at the pessimistic inflation factor
 * V_c = max_k D_ck or the optimistic one L_c = min_k D_ck, which are the
 * largest and smallest delays one customer can inflict. The level-0
 * multiple-class PBH bounds on the mean response time are</p>
 *
 * <pre>J_j(n) = sum_k D_jk,   B_j(n) = sum_k D_jk + (sum(n) - 1) max_k D_jk,</pre>
 *
 * <p>i.e. an arriving customer queues behind nobody, respectively behind every
 * other customer in the network at its own worst centre.</p>
 */
public final class Pfqn_looping {
    private Pfqn_looping() {}

    /** Result of the Looping iteration. */
    public static final class Result {
        /** Pessimistic (lower) throughput bound, one entry per class. */
        public Matrix Xlo;
        /** Optimistic (upper) throughput bound, one entry per class. */
        public Matrix Xup;
        /** Mean queue lengths on the pessimistic side. */
        public Matrix Q;
        /** Residence times. */
        public Matrix R;
        /** Iterations performed. */
        public int iter;
    }

    public static Result pfqn_looping(Matrix L, Matrix N, Matrix Z) {
        return pfqn_looping(L, N, Z, 1e-6, 1000);
    }

    public static Result pfqn_looping(Matrix L, Matrix N, Matrix Zin, double tol, int maxiter) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix Z = (Zin == null || Zin.isEmpty()) ? new Matrix(1, R) : Zin;

        double[] Dtot = new double[R];
        double[] Vpess = new double[R];
        double[] Lopt = new double[R];
        double Ntot = 0.0;
        for (int c = 0; c < R; c++) {
            double s = 0.0;
            double mx = Double.NEGATIVE_INFINITY;
            double mn = Double.POSITIVE_INFINITY;
            for (int i = 0; i < M; i++) {
                double d = L.get(i, c);
                s += d;
                if (d > mx) mx = d;
                if (d < mn) mn = d;
            }
            Dtot[c] = s;
            Vpess[c] = (M > 0) ? mx : 0.0;
            Lopt[c] = (M > 0) ? mn : 0.0;
            Ntot += N.get(c);
        }
        // level-0 multiple-class PBH bounds on R_j, at N - 1_c
        double[] Jbnd = new double[R];
        double[] Bm = new double[R];
        for (int c = 0; c < R; c++) {
            Jbnd[c] = Dtot[c];
            Bm[c] = Dtot[c] + Math.max(Ntot - 2, 0) * Vpess[c];
        }

        // Qm[c][j][k] = Q_jk(N - 1_c)
        double[][][] Qm = new double[R][R][M];
        for (int c = 0; c < R; c++) {
            for (int j = 0; j < R; j++) {
                double nj = Math.max(N.get(j) - (c == j ? 1.0 : 0.0), 0.0);
                for (int k = 0; k < M; k++) {
                    Qm[c][j][k] = (M > 0) ? nj / M : 0.0;
                }
            }
        }
        double[][] Hopt = new double[R][R];    // Hopt[j][c]
        double[][] Hpess = new double[R][R];

        Matrix QN = new Matrix(M, R);
        Matrix RN = new Matrix(M, R);
        Matrix XN = new Matrix(1, R);
        double[] Rc = new double[R];
        double[] Rpess = new double[R];
        double[] Ropt = new double[R];

        int it = 1;
        while (it <= maxiter) {
            Matrix QN_1 = QN.copy();
            for (int c = 0; c < R; c++) {
                if (N.get(c) == 0.0) {
                    for (int k = 0; k < M; k++) {
                        RN.set(k, c, 0.0);
                    }
                    XN.set(c, 0.0);
                    Rc[c] = 0.0;
                    Rpess[c] = 0.0;
                    continue;
                }
                double rsum = 0.0;
                for (int k = 0; k < M; k++) {
                    double qk = 0.0;
                    for (int j = 0; j < R; j++) {
                        qk += Qm[c][j][k];
                    }
                    RN.set(k, c, L.get(k, c) * (1 + qk));
                    rsum += RN.get(k, c);
                }
                Rc[c] = rsum;
                double hp = 0.0;
                for (int j = 0; j < R; j++) {
                    hp += Hpess[j][c];
                }
                Rpess[c] = Rc[c] + Vpess[c] * hp;
                XN.set(c, N.get(c) / (Z.get(c) + Rpess[c]));
            }
            for (int c = 0; c < R; c++) {
                if (N.get(c) == 0.0) {
                    Ropt[c] = 0.0;
                    continue;
                }
                double sat = Double.NEGATIVE_INFINITY;
                for (int k = 0; k < M; k++) {
                    double used = 0.0;
                    for (int j = 0; j < R; j++) {
                        if (j != c) {
                            used += XN.get(j) * L.get(k, j);
                        }
                    }
                    double den = 1 - used;
                    if (den > 0) {
                        double v = L.get(k, c) * N.get(c) / den - Z.get(c);
                        if (v > sat) sat = v;
                    }
                }
                double ho = 0.0;
                for (int j = 0; j < R; j++) {
                    ho += Hopt[j][c];
                }
                double heaped = Rc[c] + Lopt[c] * ho;
                double v = Math.max(Math.max(sat, heaped), Dtot[c]);
                Ropt[c] = Math.min(v, Rpess[c]);
            }
            for (int c = 0; c < R; c++) {
                for (int k = 0; k < M; k++) {
                    QN.set(k, c, XN.get(c) * RN.get(k, c));
                }
            }
            for (int c = 0; c < R; c++) {
                for (int j = 0; j < R; j++) {
                    double nj = N.get(j) - (c == j ? 1.0 : 0.0);
                    double qsum = 0.0;
                    if (N.get(j) <= 0 || nj <= 0) {
                        for (int k = 0; k < M; k++) {
                            Qm[c][j][k] = 0.0;
                        }
                    } else {
                        double f = (nj / N.get(j)) * ((Z.get(j) + Ropt[j]) / (Z.get(j) + Bm[j]));
                        for (int k = 0; k < M; k++) {
                            Qm[c][j][k] = f * QN.get(k, j);
                            qsum += Qm[c][j][k];
                        }
                    }
                    if (nj > 0) {
                        Hopt[j][c] = Math.max(0, Jbnd[j] / (Z.get(j) + Jbnd[j]) * nj - qsum);
                        Hpess[j][c] = Math.max(0, Bm[j] / (Z.get(j) + Bm[j]) * nj - qsum);
                    } else {
                        Hopt[j][c] = 0.0;
                        Hpess[j][c] = 0.0;
                    }
                }
            }
            double maxdiff = 0.0;
            boolean anyNonEmpty = false;
            for (int c = 0; c < R; c++) {
                if (N.get(c) <= 0) continue;
                anyNonEmpty = true;
                for (int k = 0; k < M; k++) {
                    maxdiff = Math.max(maxdiff, Math.abs(QN.get(k, c) - QN_1.get(k, c)));
                }
            }
            if (!anyNonEmpty || (it > 1 && maxdiff < tol)) {
                break;
            }
            it++;
        }

        Result out = new Result();
        out.Xlo = XN;
        out.Xup = new Matrix(1, R);
        for (int c = 0; c < R; c++) {
            if (N.get(c) > 0) {
                out.Xup.set(c, N.get(c) / (Z.get(c) + Ropt[c]));
            }
        }
        out.Q = QN;
        out.R = RN;
        out.iter = it;
        return out;
    }
}
