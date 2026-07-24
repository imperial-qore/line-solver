/**
 * Exact per-station queue-length variances and covariances for closed
 * product-form queueing networks, computed by an MVA-like moment recursion that
 * does not require the full sensitivity Jacobian.
 *
 * <p>Mirrors the MATLAB reference {@code pfqn_sens_mva.m}.</p>
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.sens;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_sens_mva {
    private Pfqn_sens_mva() {}

    /**
     * Exact second moments (variances and per-station covariances) of the queue
     * lengths of a closed product-form (BCMP) queueing network. The moments are
     * obtained by an MVA-type recursion evaluated on the same population lattice
     * as {@code pfqn_mva}, so no derivative of the model is ever formed and the
     * cost is O(M*R^2) per lattice point rather than the O(M^2*R^2) of the
     * differentiated-MVA algorithm used by {@link Pfqn_sens}.
     *
     * <p>The recursion is obtained by differentiating the Reiser-Lavenberg MVA
     * equation Q(j,v|N) = X(v|N) * L(j,v) * (mi(j) + Qtot(j|N-e_v)) with respect
     * to the visit ratio theta(i,k) of class k at station i and rescaling.
     * Writing W(k,i;v,j|N) = Cov[n(i,k),n(j,v)] at population N,</p>
     *
     * <pre>
     *   W(k,i;v,j|N) = Q(j,v|N) * ( Q(i,k|N-e_v) - Q(i,k|N) )
     *                + [i==j &amp; k==v] * Q(j,v|N)
     *                + X(v|N) * L(j,v) * sum_t W(k,i;t,j|N-e_v)
     * </pre>
     *
     * <p>with W(.|0) = 0. This routine evaluates the same-station case i==j,
     * which is self-contained: the inner sum then only involves same-station
     * terms, so a single scalar Ssum(j,k|N) = sum_t W(k,j;t,j|N) carried along
     * the lattice closes the recursion. The cross-station case i!=j is not
     * self-contained (it couples every station pair) and costs as much as the
     * full Jacobian, so it is left to {@link Pfqn_sens}.</p>
     *
     * <p>Setting i==j and mi==1 reproduces Corollary 1 of the reference below,
     * i.e. its equations (2.9a) for the variance and (2.10) for the covariance;
     * the station multiplicity mi cancels identically because
     * X(v|N)*L(j,v)*(mi(j)+Qtot(j|N-e_v)) = Q(j,v|N) is the MVA equation for any
     * mi. The equivalent statement for an infinite-server station, equation
     * (2.9b) of the reference, is recovered automatically because LINE folds the
     * delay into the think time Z, which enters only through X(v|N) and carries
     * no queue-length moment of its own.</p>
     *
     * <p>Reference: E. de Souza e Silva and R. R. Muntz, "Simple Relationships
     * Among Moments of Queue Lengths in Product Form Queueing Networks", IEEE
     * Trans. Computers 37(9):1125-1129, 1988 (Theorems 1-2 and Corollary 1). The
     * underlying identity Cov[n(i,k),n(j,v)] = theta(i,k) * dQ(j,v)/dtheta(i,k)
     * is the k=2 case of Theorem 1 of I. F. Akyildiz and J. C. Strelen, "Moment
     * Analysis for Load-Dependent Mixed Product Form Queueing Networks", IEEE
     * Trans. Communications 39(6):828-832, 1991.</p>
     *
     * <p>Restricted to closed populations. Mixed and load-dependent models are
     * handled by {@link Pfqn_sens_mvaldmx}. For a station of multiplicity mi(i)&gt;1,
     * which LINE treats as mi(i) identical replicas sharing the demand row
     * L(i,:), the moments returned are those of the aggregate queue length over
     * the replicas.</p>
     *
     * @param L  service demand matrix (M x R), L(i,r) = visits_ir / rate_ir
     * @param N  population vector (1 x R)
     * @param Z  think time vector (1 x R), null or empty for zeros
     * @param mi station server multiplicity (1 x M), null for ones
     * @return the base measures and their exact second moments
     */
    public static Ret.pfqnSensMva pfqn_sens_mva(Matrix L, Matrix N, Matrix Z, Matrix mi) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        N = N.copy().ceil();
        if (N.getNumRows() > 1) {
            N = N.transpose();
        }
        if (Z == null || Z.isEmpty()) {
            Z = new Matrix(1, R);
        } else if (Z.getNumRows() > 1) {
            Z = Z.transpose();
        }
        if (mi == null) {
            mi = Matrix.ones(1, M);
        } else if (mi.getNumRows() > 1) {
            mi = mi.transpose();
        }
        if (N.length() != R) {
            throw new RuntimeException("pfqn_sens_mva: demand matrix and population vector have "
                    + "different number of classes");
        }
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(N.get(0, r))) {
                throw new RuntimeException("pfqn_sens_mva: requires a closed population; use "
                        + "pfqn_sens_mvaldmx for mixed models");
            }
        }

        Matrix X = new Matrix(1, R);
        Matrix Q = new Matrix(M, R);
        Matrix U = new Matrix(M, R);
        Matrix C = new Matrix(M, R);
        double[][][] QCov = new double[M][R][R];

        if (!N.any()) {
            return pack(X, Q, U, C, QCov);
        }

        // see _kb/03-api-layer.md for rationale
        Matrix prods = new Matrix(1, Math.max(0, R - 1));
        for (int w = 0; w < R - 1; w++) {
            double acc = 1.0;
            for (int i = 0; i < R - (w + 2) + 1; i++) {
                acc *= (1.0 + N.get(0, w + 1 + i));
            }
            prods.set(0, w, acc);
        }
        int firstNonEmpty = R - 1;
        while (N.get(0, firstNonEmpty) == 0.0) {
            firstNonEmpty--;
        }
        double totpop = 1.0;
        for (int r = 0; r < R; r++) {
            totpop *= (N.get(0, r) + 1.0);
        }
        int TP = (int) totpop;
        double ctr = totpop;
        double[][] Qtot = new double[TP][M];        // Qtot[m][i]    = sum_r Q(i,r) at population m
        double[][][] Qcls = new double[TP][M][R];   // Qcls[m][i][r] = Q(i,r) at population m
        double[][] Xall = new double[TP][R];        // Xall[m][r]    = X(r) at population m
        double[][][] Ssum = new double[TP][M][R];   // Ssum[m][j][k] = sum_t Cov[n(j,k),n(j,t)]
        int currentpop = 1;

        Matrix n = new Matrix(1, R);
        n.set(0, firstNonEmpty, 1);
        int[] rows = new int[R];                    // rows[s] = lattice index of n - e_s

        while (ctr > 0) {
            // ---- mean value analysis step at population n -------------------
            int s = 0;
            while (s < R) {
                int pos = 0;
                if (n.get(0, s) > 0) {
                    n.set(0, s, n.get(0, s) - 1);
                    pos = (int) n.get(0, R - 1);
                    int w = 0;
                    while (w < R - 1) {
                        pos = (int) (pos + n.get(0, w) * prods.get(0, w));
                        w++;
                    }
                    n.set(0, s, n.get(0, s) + 1);
                }
                // see _kb/03-api-layer.md for rationale
                rows[s] = pos;
                double CNtot = 0.0;
                for (int i = 0; i < M; i++) {
                    C.set(i, s, L.get(i, s) * (mi.get(0, i) + Qtot[pos][i]));
                    CNtot += C.get(i, s);
                }
                double den = Z.get(0, s) + CNtot;
                double xs = n.get(0, s) / den;
                X.set(0, s, xs);
                Xall[currentpop][s] = xs;
                for (int i = 0; i < M; i++) {
                    double q = xs * C.get(i, s);
                    Q.set(i, s, q);
                    Qcls[currentpop][i][s] = q;
                    Qtot[currentpop][i] += q;
                }
                s++;
            }

            // see _kb/03-api-layer.md for rationale
            for (int j = 0; j < M; j++) {
                for (int k = 0; k < R; k++) {
                    double Qjk = Qcls[currentpop][j][k];
                    double sk = 0.0;
                    for (int t = 0; t < R; t++) {
                        double Qjt = Qcls[currentpop][j][t];
                        double wkt = Qjt * (Qcls[rows[t]][j][k] - Qjk);
                        if (k == t) {
                            wkt += Qjt;
                        }
                        wkt += Xall[currentpop][t] * L.get(j, t) * Ssum[rows[t]][j][k];
                        QCov[j][k][t] = wkt;
                        sk += wkt;
                    }
                    Ssum[currentpop][j][k] = sk;
                }
            }

            // ---- odometer advance -------------------------------------------
            s = R - 1;
            while ((s >= 0 && (n.get(0, s) == N.get(0, s))) || s > firstNonEmpty) {
                s--;
            }
            if (s == -1) {
                break;
            }
            n.set(0, s, n.get(0, s) + 1);
            s++;
            while (s < R) {
                n.set(0, s, 0);
                s++;
            }
            ctr--;
            currentpop++;
        }

        // utilization
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                U.set(i, r, X.get(0, r) * L.get(i, r));
            }
        }

        return pack(X, Q, U, C, QCov);
    }

    public static Ret.pfqnSensMva pfqn_sens_mva(Matrix L, Matrix N, Matrix Z) {
        return pfqn_sens_mva(L, N, Z, null);
    }

    public static Ret.pfqnSensMva pfqn_sens_mva(Matrix L, Matrix N) {
        return pfqn_sens_mva(L, N, null, null);
    }

    // =====================================================================
    private static Ret.pfqnSensMva pack(Matrix X, Matrix Q, Matrix U, Matrix C, double[][][] raw) {
        int M = Q.getNumRows();
        int R = Q.getNumCols();
        // see _kb/03-api-layer.md for rationale
        double asym = 0.0;
        Matrix[] QCov = new Matrix[M];
        for (int i = 0; i < M; i++) {
            Matrix cov = new Matrix(R, R);
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    asym = Math.max(asym, Math.abs(raw[i][r][s] - raw[i][s][r]));
                    cov.set(r, s, (raw[i][r][s] + raw[i][s][r]) / 2.0);
                }
            }
            QCov[i] = cov;
        }
        Matrix QVar = new Matrix(M, R);
        Matrix QTotVar = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            double tot = 0.0;
            for (int r = 0; r < R; r++) {
                QVar.set(i, r, QCov[i].get(r, r));
                for (int s = 0; s < R; s++) {
                    tot += QCov[i].get(r, s);
                }
            }
            QTotVar.set(i, 0, tot);
        }
        return new Ret.pfqnSensMva(X, Q, U, C, QCov, QVar, QTotVar, asym);
    }
}
