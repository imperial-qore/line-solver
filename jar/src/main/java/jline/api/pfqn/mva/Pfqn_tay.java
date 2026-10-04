/**
 * @file Tay's arrival-instant approximate mean value analysis
 *
 * Approximate MVA for closed multiclass product-form networks in which the arrival-instant
 * queue lengths are estimated from the THROUGHPUT ELASTICITIES rather than from a
 * population-shift heuristic (Tay 1987; eqs. 4.8.2-1..3 of the Schweitzer-Serazzi-Broglia
 * survey). Ported at parity from MATLAB pfqn_tay.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_tay {
    private Pfqn_tay() {}

    /** Mean performance measures plus the arrival-instant queue lengths. */
    public static final class Result {
        /** Per-class throughputs (1 x R). */
        public final Matrix X;
        /** Mean queue lengths (M x R). */
        public final Matrix Q;
        /** Utilizations (M x R). */
        public final Matrix U;
        /** Residence times (M x R). */
        public final Matrix R;
        /** Number of iterations performed. */
        public final int totiter;
        /**
         * Arrival-instant queue lengths, indexed [r][m][k]: the class-k queue length at
         * station m as seen by an arriving class-r job. These are the auxiliary
         * quantities the method is tabulated on, and they are NOT the queue lengths of
         * the model re-solved at N - e_r, which is the same object only for an exact
         * solution.
         */
        public final Matrix[] Qarr;

        public Result(Matrix X, Matrix Q, Matrix U, Matrix R, int totiter, Matrix[] Qarr) {
            this.X = X;
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.totiter = totiter;
            this.Qarr = Qarr;
        }
    }

    /**
     * Tay's arrival-instant AMVA.
     *
     * <p>Let E_mkc = (D_mk/X_c) dX_c/dD_mk be the elasticity of the class-c throughput
     * with respect to the class-k demand at station m. Tay shows that the elasticities
     * satisfy the R linear equations</p>
     *
     * <pre>
     * E_mkj sum_t B_tj Q_jt (1+Q_jt) =
     *     -[(delta_jk + Q_jm) B_mk Q_km + sum_{c!=j} E_mkc sum_t B_tc Q_jt Q_ct]
     * </pre>
     *
     * <p>with B_ir = 1/(1 + D_ir X_r/N_r), and that the arrival-instant queue length is
     * then simply Q_km^(r) = Q_km + E_mkr, which closes the MVA recursion
     * R_rm = D_rm (1 + sum_k Q_km^(r)). One R x R solve per (station, class) pair per
     * iteration.</p>
     *
     * <p>Delay stations enter through Z only. They are "AS" servers in the survey's
     * notation (d_t = 0), so they contribute Z_j X_j to the denominator of the elasticity
     * equations but nothing to its numerator.</p>
     *
     * @param L       service demand matrix (M x R)
     * @param N       population vector (1 x R)
     * @param Z       think time vector (1 x R)
     * @param tol     convergence tolerance on the queue lengths
     * @param maxiter maximum number of iterations
     * @param QN0     initial guess for the queue lengths (M x R), or null
     * @return throughputs, queue lengths, utilizations, residence times, iteration count
     *         and the arrival-instant queue lengths
     */
    public static Result pfqn_tay(Matrix L, Matrix N, Matrix Z, double tol, int maxiter, Matrix QN0) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        double[] Nv = new double[R];
        double[] Zv = new double[R];
        for (int r = 0; r < R; r++) {
            Nv[r] = N.get(r);
            double zr = 0.0;
            if (Z != null) {
                if (Z.getNumRows() > 1) {
                    for (int i = 0; i < Z.getNumRows(); i++) {
                        zr += Z.get(i, r);   // several delay stations aggregate
                    }
                } else if (Z.length() > r) {
                    zr = Z.get(r);
                }
            }
            Zv[r] = zr;
        }

        Matrix X = new Matrix(1, R);
        Matrix Q = new Matrix(M, R);
        Matrix U = new Matrix(M, R);
        Matrix Res = new Matrix(M, R);
        Matrix[] Qarr = new Matrix[R];
        for (int r = 0; r < R; r++) {
            Qarr[r] = new Matrix(M, R);
        }

        // Empty classes contribute no jobs anywhere and make the elasticity system
        // singular (their denominator is identically zero); solve without them and
        // re-expand, as in Pfqn_bs.
        int nact = 0;
        for (int r = 0; r < R; r++) {
            if (Nv[r] > 0) {
                nact++;
            }
        }
        if (nact == 0) {
            return new Result(X, Q, U, Res, 0, Qarr);
        }
        if (nact < R) {
            int[] act = new int[nact];
            int c = 0;
            for (int r = 0; r < R; r++) {
                if (Nv[r] > 0) {
                    act[c++] = r;
                }
            }
            Matrix La = new Matrix(M, nact);
            Matrix Na = new Matrix(1, nact);
            Matrix Za = new Matrix(1, nact);
            for (int q = 0; q < nact; q++) {
                Na.set(0, q, Nv[act[q]]);
                Za.set(0, q, Zv[act[q]]);
                for (int i = 0; i < M; i++) {
                    La.set(i, q, L.get(i, act[q]));
                }
            }
            Result sub = pfqn_tay(La, Na, Za, tol, maxiter, null);
            for (int q = 0; q < nact; q++) {
                X.set(0, act[q], sub.X.get(q));
                for (int i = 0; i < M; i++) {
                    Q.set(i, act[q], sub.Q.get(i, q));
                    U.set(i, act[q], sub.U.get(i, q));
                    Res.set(i, act[q], sub.R.get(i, q));
                    for (int k = 0; k < nact; k++) {
                        Qarr[act[q]].set(i, act[k], sub.Qarr[q].get(i, k));
                    }
                }
            }
            return new Result(X, Q, U, Res, sub.totiter, Qarr);
        }

        if (QN0 == null) {
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    Q.set(i, r, Nv[r] / M);
                }
            }
        } else {
            Q = QN0.copy();
        }
        for (int r = 0; r < R; r++) {
            double sumL = 0.0;
            double sumQ = 0.0;
            for (int i = 0; i < M; i++) {
                sumL += L.get(i, r);
                sumQ += Q.get(i, r);
            }
            X.set(0, r, Nv[r] / (Zv[r] + sumL * (1 + sumQ)));
        }

        double[][] B = new double[M][R];
        double[] den = new double[R];
        double[][] C = new double[R][R];
        Matrix A = new Matrix(R, R);
        Matrix b = new Matrix(R, 1);
        Matrix E = new Matrix(R, 1);
        int it = 0;

        for (it = 1; it <= maxiter; it++) {
            Matrix Qprev = Q.copy();

            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    B[i][r] = 1.0 / (1.0 + L.get(i, r) * X.get(r) / Nv[r]);
                }
            }

            // Denominators of the elasticity system, one per class; the delay term
            // Z_j X_j is the AS-server contribution (d_t = 0 leaves B = 1).
            for (int j = 0; j < R; j++) {
                double acc = 0.0;
                for (int i = 0; i < M; i++) {
                    acc += B[i][j] * Q.get(i, j) * (1 + Q.get(i, j));
                }
                den[j] = acc + Zv[j] * X.get(j);
            }

            for (int j = 0; j < R; j++) {
                for (int cc = 0; cc < R; cc++) {
                    double acc = 0.0;
                    for (int i = 0; i < M; i++) {
                        acc += B[i][cc] * Q.get(i, j) * Q.get(i, cc);
                    }
                    C[j][cc] = acc;
                }
            }

            for (int m = 0; m < M; m++) {
                for (int k = 0; k < R; k++) {
                    for (int j = 0; j < R; j++) {
                        for (int cc = 0; cc < R; cc++) {
                            A.set(j, cc, cc == j ? 1.0 : C[j][cc] / den[j]);
                        }
                        double djk = (j == k) ? 1.0 : 0.0;
                        b.set(j, 0, -(djk + Q.get(m, j)) * B[m][k] * Q.get(m, k) / den[j]);
                    }
                    Matrix.solve(A, b, E);
                    for (int r = 0; r < R; r++) {
                        Qarr[r].set(m, k, Q.get(m, k) + E.get(r, 0));
                    }
                }
            }

            for (int r = 0; r < R; r++) {
                double sumR = 0.0;
                for (int i = 0; i < M; i++) {
                    double seen = 0.0;
                    for (int k = 0; k < R; k++) {
                        seen += Qarr[r].get(i, k);
                    }
                    double ri = L.get(i, r) * (1 + seen);
                    Res.set(i, r, ri);
                    sumR += ri;
                }
                X.set(0, r, Nv[r] / (Zv[r] + sumR));
            }
            double delta = 0.0;
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    Q.set(i, r, X.get(r) * Res.get(i, r));
                    delta = Math.max(delta, Math.abs(Q.get(i, r) - Qprev.get(i, r)));
                }
            }
            if (delta < tol) {
                break;
            }
        }
        if (it > maxiter) {
            it = maxiter;
        }
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                U.set(i, r, L.get(i, r) * X.get(r));
            }
        }
        return new Result(X, Q, U, Res, it, Qarr);
    }

    /** Convenience overload with the default tolerance and iteration budget. */
    public static Result pfqn_tay(Matrix L, Matrix N, Matrix Z) {
        return pfqn_tay(L, N, Z, 1e-6, 1000, null);
    }
}
