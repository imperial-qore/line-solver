/**
 * Moment linearizer: approximate mean queue lengths and their second moments
 * (variance / covariance) for large closed product-form queueing networks.
 *
 * <p>Means come from the Schweitzer-Bard AMVA fixed point. Second moments use
 * the exact product-form identity Cov[n_{i,r},n_{j,s}] = D_{j,s} dQ_{i,r}/dD_{j,s},
 * with the demand derivatives obtained by analytically linearizing the AMVA
 * fixed point (a "moment linearizer" in the sense of Strelen and Akyildiz,
 * developed per class). Both moments carry the AMVA approximation error; for
 * exact moments on tractable models use {@link jline.api.pfqn.sens.Pfqn_sens}.</p>
 *
 * <p>Mirrors the MATLAB reference {@code pfqn_momlin.m} and native-Python
 * {@code pfqn_momlin}.</p>
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_momlin {
    private Pfqn_momlin() {}

    /** Result of {@link #pfqn_momlin}. */
    public static final class MomlinResult {
        public Matrix Q;          // M x R mean queue length
        public Matrix X;          // 1 x R throughput per class
        public Matrix U;          // M x R utilization
        public Matrix R;          // M x R residence time
        public Matrix QVar;       // M x R queue-length variance
        public Matrix[][] QCov;   // QCov[i][r] is M x R : Cov[n_{i,r},n_{j,s}]
        public Matrix[][] dQ;     // dQ[j][s] is M x R : dQ_{i,r}/dD_{j,s}
    }

    public static MomlinResult pfqn_momlin(Matrix L, Matrix N, Matrix Z) {
        return pfqn_momlin(L, N, Z, 1e-8, 1000);
    }

    /**
     * Approximate first and second queue-length moments of a closed
     * product-form network.
     *
     * @param L       service demand matrix (M x R)
     * @param N       closed population vector (1 x R)
     * @param Z       think time vector (1 x R), may be null
     * @param tol     convergence tolerance on the queue-length fixed point
     * @param maxiter maximum iterations
     */
    public static MomlinResult pfqn_momlin(Matrix L, Matrix N, Matrix Z, double tol, int maxiter) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double[] Nn = new double[R];
        for (int r = 0; r < R; r++) {
            Nn[r] = Math.ceil(N.get(N.getNumRows() > 1 ? r : 0, N.getNumRows() > 1 ? 0 : r));
        }
        double[] Zr = new double[R];
        if (Z != null && !Z.isEmpty()) {
            for (int r = 0; r < R; r++) {
                Zr[r] = Z.get(Z.getNumRows() > 1 ? r : 0, Z.getNumRows() > 1 ? 0 : r);
            }
        }
        double[][] Ld = new double[M][R];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Ld[i][r] = L.get(i, r);
            }
        }

        // see _kb/03-api-layer.md for rationale
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(Nn[r]))
                throw new IllegalArgumentException(
                        "pfqn_momlin supports closed classes only");
        }

        // Schweitzer coefficients c[r][s] = (N_s - delta_{rs})/N_s
        double[][] c = new double[R][R];
        for (int r = 0; r < R; r++) {
            for (int s = 0; s < R; s++) {
                if (Nn[s] > 0) {
                    c[r][s] = (Nn[s] - (r == s ? 1.0 : 0.0)) / Nn[s];
                } else {
                    c[r][s] = 0.0;
                }
            }
        }

        // ---- Schweitzer-Bard AMVA fixed point for the means -------------
        double[][] Q = new double[M][R];
        for (int r = 0; r < R; r++) {
            if (Nn[r] > 0) {
                for (int i = 0; i < M; i++) {
                    Q[i][r] = Nn[r] / M;
                }
            }
        }
        double[] X = new double[R];
        double[][] Rr = new double[M][R];
        for (int it = 0; it < maxiter; it++) {
            double diff = 0.0;
            for (int r = 0; r < R; r++) {
                if (Nn[r] == 0) {
                    X[r] = 0.0;
                    for (int i = 0; i < M; i++) {
                        Rr[i][r] = 0.0;
                    }
                    continue;
                }
                double sumR = 0.0;
                for (int i = 0; i < M; i++) {
                    double a = 0.0;
                    for (int s = 0; s < R; s++) {
                        a += c[r][s] * Q[i][s];
                    }
                    Rr[i][r] = Ld[i][r] * (1.0 + a);
                    sumR += Rr[i][r];
                }
                X[r] = Nn[r] / (Zr[r] + sumR);
                for (int i = 0; i < M; i++) {
                    double qnew = X[r] * Rr[i][r];
                    diff = Math.max(diff, Math.abs(qnew - Q[i][r]));
                    Q[i][r] = qnew;
                }
            }
            if (diff < tol) {
                break;
            }
        }

        // ---- analytic linearization of the fixed point ------------------
        // dQ[j][s0] holds dQ_{i,r}/dD_{j,s0} over (i,r).
        Matrix[][] dQ = new Matrix[M][R];
        for (int j = 0; j < M; j++) {
            for (int s0 = 0; s0 < R; s0++) {
                dQ[j][s0] = new Matrix(M, R);
            }
        }
        for (int j = 0; j < M; j++) {
            for (int s0 = 0; s0 < R; s0++) {
                if (Nn[s0] == 0) {
                    continue;
                }
                double[][] dq = new double[M][R];
                for (int it = 0; it < maxiter; it++) {
                    double diff = 0.0;
                    for (int r = 0; r < R; r++) {
                        if (Nn[r] == 0) {
                            continue;
                        }
                        double[] dR = new double[M];
                        double sumdR = 0.0;
                        for (int i = 0; i < M; i++) {
                            double dDir = (i == j && r == s0) ? 1.0 : 0.0;
                            double aQ = 0.0;
                            double adq = 0.0;
                            for (int s = 0; s < R; s++) {
                                aQ += c[r][s] * Q[i][s];
                                adq += c[r][s] * dq[i][s];
                            }
                            dR[i] = dDir * (1.0 + aQ) + Ld[i][r] * adq;
                            sumdR += dR[i];
                        }
                        double dXr = -(X[r] * X[r] / Nn[r]) * sumdR;
                        for (int i = 0; i < M; i++) {
                            double qnew = dXr * Rr[i][r] + X[r] * dR[i];
                            diff = Math.max(diff, Math.abs(qnew - dq[i][r]));
                            dq[i][r] = qnew;
                        }
                    }
                    if (diff < tol) {
                        break;
                    }
                }
                for (int i = 0; i < M; i++) {
                    for (int r = 0; r < R; r++) {
                        dQ[j][s0].set(i, r, dq[i][r]);
                    }
                }
            }
        }

        // ---- second moments via the covariance identity -----------------
        Matrix[][] QCov = new Matrix[M][R];
        Matrix QVar = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Matrix cov = new Matrix(M, R);
                for (int j = 0; j < M; j++) {
                    for (int s = 0; s < R; s++) {
                        cov.set(j, s, Ld[j][s] * dQ[j][s].get(i, r));
                    }
                }
                QCov[i][r] = cov;
                QVar.set(i, r, cov.get(i, r));
            }
        }

        MomlinResult res = new MomlinResult();
        res.Q = new Matrix(M, R);
        res.U = new Matrix(M, R);
        res.R = new Matrix(M, R);
        res.X = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            res.X.set(0, r, X[r]);
            for (int i = 0; i < M; i++) {
                res.Q.set(i, r, Q[i][r]);
                res.R.set(i, r, Rr[i][r]);
                res.U.set(i, r, X[r] * Ld[i][r]);
            }
        }
        res.QVar = QVar;
        res.QCov = QCov;
        res.dQ = dQ;
        return res;
    }
}
