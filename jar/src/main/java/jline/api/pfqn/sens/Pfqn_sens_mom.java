/**
 * Exact higher moments (up to order three) of the per-station total queue lengths
 * of a closed product-form queueing network, by second-order differentiation of
 * the MVA recursion.
 *
 * <p>Mirrors the MATLAB reference {@code pfqn_sens_mom.m}.</p>
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.sens;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_sens_mom {
    private Pfqn_sens_mom() {}

    /**
     * Exact moments E[Q_i], E[Q_i^2], E[Q_i^3] and the covariances Cov[Q_i,Q_j] of
     * the TOTAL queue lengths Q_i = sum_r n(i,r) of a closed product-form (BCMP)
     * queueing network.
     *
     * <p>The method is the moment analysis of Strelen. Its Theorem 3.1 states that
     * one further factor Q_i in a moment costs one differentiation with respect to
     * x_i, the reciprocal of the capacity of station i, i.e. a parameter that
     * scales the service times of ALL classes at station i:</p>
     *
     * <pre>
     *   E[Q_i^j] = m_i E[Q_i^(j-1)] + x_i d/dx_i E[Q_i^(j-1)],  Q_i^0 = 1
     * </pre>
     *
     * <p>Iterating from E[Q_i^0] = 1 gives, with m_i = E[Q_i] (equation (3.2)):</p>
     *
     * <pre>
     *   Var[Q_i]     = x_i dm_i/dx_i
     *   Cov[Q_i,Q_j] = x_j dm_i/dx_j = x_i dm_j/dx_i
     *   E[Q_i^2]     = x_i dm_i/dx_i + m_i^2
     *   E[Q_i^3]     = x_i^2 d^2m_i/dx_i^2 + (x_i + 3 x_i m_i) dm_i/dx_i + m_i^3
     * </pre>
     *
     * <p>so the third moment requires the SECOND derivative of the MVA recursion,
     * which is what this routine adds over {@link Pfqn_sens_mva} and
     * {@link Pfqn_sens} (both first order only). The derivatives are obtained by
     * second-order forward-mode differentiation of the Reiser-Lavenberg recursion,
     * i.e. by carrying, for each parameter, the value together with its first and
     * second derivative along the population lattice (Theorem 3.2 for one class,
     * Theorem 3.5 for several).</p>
     *
     * <p>The parameter need not scale a whole column. Theorem 1 of Akyildiz and
     * Strelen states the same recursion for a parameter that scales the service
     * times of an arbitrary class subset T at station i, and the moments it
     * generates are then those of Q_(i,T) = sum_(r in T) n(i,r). {@code groups}
     * supplies that subset structure: it partitions the classes, and the routine
     * reports the moments of each group's queue length at each station. The three
     * useful settings are</p>
     *
     * <ul>
     *   <li>{@code groups = ones(1,R)} -- the whole column: per-station TOTALS
     *       (the default, and Strelen's x_i);</li>
     *   <li>{@code groups = 1..R} -- one class per group: PER-CLASS moments, so
     *       that even the third moment is per class (see {@link #perClassGroups});</li>
     *   <li>{@code groups = chain(r)} -- one group per chain: PER-CHAIN moments.</li>
     * </ul>
     *
     * <p>Strelen states only the first; the generalization is Akyildiz and
     * Strelen's T. The per-class setting reproduces {@link Pfqn_sens_mva}'s
     * {@code QVar} exactly, and brute-force enumeration of the per-class moments
     * including the third.</p>
     *
     * <p>Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks and
     * its Linearizer", Performance Evaluation 11:127-142, 1990, Theorems 2.1, 3.1,
     * 3.2, 3.5 and equation (3.2); I. F. Akyildiz and J. C. Strelen, "Moment
     * Analysis for Load-Dependent Mixed Product Form Queueing Networks", IEEE Trans.
     * Communications 39(6):828-832, 1991, Theorem 1.</p>
     *
     * <p>Restricted to closed populations, as is the moment analysis of the
     * reference. Mixed and load-dependent second moments are in
     * {@link Pfqn_sens_mvaldmx}. Moments of the sojourn times at FCFS centers are
     * built on top of these queue-length moments by {@link Pfqn_sens_respt},
     * following Theorem 4.1. The exact recursion costs O(prod(N+1)) lattice
     * points; {@link Pfqn_sens_linearizer} approximates the same quantities in
     * polynomial time.</p>
     *
     * @param L  service demand matrix (M x R), L(i,r) = visits_ir / rate_ir
     * @param N  population vector (1 x R)
     * @param Z  think time vector (1 x R), null or empty for zeros
     * @param mi station server multiplicity (1 x M), null for ones
     * @param groups class-to-group map (1 x R), a partition of the classes into
     *               G = max(groups) groups labelled consecutively 1..G with no empty
     *               group. The moments returned are those of each group's queue
     *               length at each station. Null selects {@code ones(1,R)}, i.e. one
     *               group holding every class, which is the per-station total.
     * @return the base measures together with the exact moments of each group's
     *         queue length at each station
     */
    public static Ret.pfqnSensMom pfqn_sens_mom(Matrix L, Matrix N, Matrix Z, Matrix mi,
                                                Matrix groups) {
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
        if (groups == null || groups.isEmpty()) {
            groups = Matrix.ones(1, R);
        } else if (groups.getNumRows() > 1) {
            groups = groups.transpose();
        }
        if (N.length() != R) {
            throw new RuntimeException("pfqn_sens_mom: demand matrix and population vector have "
                    + "different number of classes");
        }
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(N.get(0, r))) {
                throw new RuntimeException("pfqn_sens_mom: requires a closed population");
            }
        }
        // see _kb/03-api-layer.md for rationale
        int[] grp = new int[R];
        int G = 0;
        for (int r = 0; r < R; r++) {
            grp[r] = (int) Math.round(groups.get(0, r));
            if (grp[r] < 1) {
                throw new RuntimeException("pfqn_sens_mom: groups must be a (1 x R) vector of "
                        + "group labels starting at 1");
            }
            if (grp[r] > G) {
                G = grp[r];
            }
        }
        if (groups.length() != R) {
            throw new RuntimeException("pfqn_sens_mom: groups must be a (1 x R) vector of group "
                    + "labels starting at 1");
        }
        boolean[] seen = new boolean[G + 1];
        for (int r = 0; r < R; r++) {
            seen[grp[r]] = true;
        }
        for (int g = 1; g <= G; g++) {
            if (!seen[g]) {
                throw new RuntimeException("pfqn_sens_mom: groups must label the classes "
                        + "consecutively from 1 to max(groups), with no empty group");
            }
        }

        Matrix X = new Matrix(1, R);
        Matrix Q = new Matrix(M, R);
        Matrix U = new Matrix(M, R);
        Matrix C = new Matrix(M, R);
        Matrix m = new Matrix(M, G);
        // dm[i][g] is an M x G matrix holding d m_(i,g) / d y_(j,g2) over (j,g2)
        Matrix[][] dm = new Matrix[M][G];
        for (int i = 0; i < M; i++) {
            for (int g = 0; g < G; g++) {
                dm[i][g] = new Matrix(M, G);
            }
        }
        Matrix d2m = new Matrix(M, G);

        if (!N.any()) {
            return pack(X, Q, U, C, m, dm, d2m, G);
        }

        // population-lattice odometer, identical to pfqn_mva
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

        // see _kb/03-api-layer.md for rationale
        int P = M * G;
        int[][] pidx = new int[M][G];
        for (int h = 0; h < M; h++) {
            for (int g = 0; g < G; g++) {
                pidx[h][g] = h * G + g;
            }
        }
        double[][] Qtot = new double[TP][M];
        double[][][] D1 = new double[TP][M][P];
        double[][][] D2 = new double[TP][M][P];
        double[][] Qg = new double[M][G];
        double[][][] D1Qg = new double[M][G][P];
        double[][][] D2Qg = new double[M][G][P];
        int currentpop = 1;

        Matrix n = new Matrix(1, R);
        n.set(0, firstNonEmpty, 1);

        double[] Cs = new double[M];
        double[][] dCs = new double[M][P];
        double[][] d2Cs = new double[M][P];
        double[] dCNtot = new double[P];
        double[] d2CNtot = new double[P];
        double[] dX = new double[P];
        double[] d2X = new double[P];

        while (ctr > 0) {
            // the group accumulators describe one population only
            for (int i = 0; i < M; i++) {
                for (int g = 0; g < G; g++) {
                    Qg[i][g] = 0.0;
                    for (int p = 0; p < P; p++) {
                        D1Qg[i][g][p] = 0.0;
                        D2Qg[i][g][p] = 0.0;
                    }
                }
            }
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
                int row = pos;

                // see _kb/03-api-layer.md for rationale
                int gs = grp[s] - 1;
                double CNtot = 0.0;
                for (int p = 0; p < P; p++) {
                    dCNtot[p] = 0.0;
                    d2CNtot[p] = 0.0;
                }
                for (int i = 0; i < M; i++) {
                    double A = mi.get(0, i) + Qtot[row][i];
                    Cs[i] = L.get(i, s) * A;
                    C.set(i, s, Cs[i]);
                    CNtot += Cs[i];
                    for (int p = 0; p < P; p++) {
                        double dA = D1[row][i][p];
                        double d2A = D2[row][i][p];
                        if (p == pidx[i][gs]) {
                            dCs[i][p] = L.get(i, s) * (A + dA);
                            d2Cs[i][p] = L.get(i, s) * (2 * dA + d2A);
                        } else {
                            dCs[i][p] = L.get(i, s) * dA;
                            d2Cs[i][p] = L.get(i, s) * d2A;
                        }
                        dCNtot[p] += dCs[i][p];
                        d2CNtot[p] += d2Cs[i][p];
                    }
                }

                // see _kb/03-api-layer.md for rationale
                double den = Z.get(0, s) + CNtot;
                double ns = n.get(0, s);
                double xs = ns / den;
                X.set(0, s, xs);
                for (int p = 0; p < P; p++) {
                    dX[p] = -ns * dCNtot[p] / (den * den);
                    d2X[p] = -ns * d2CNtot[p] / (den * den)
                            + 2 * ns * dCNtot[p] * dCNtot[p] / (den * den * den);
                }

                // ---- queue lengths ---------------------------------------------
                // Q = X*C,  dQ = dX*C + X*dC,  d2Q = d2X*C + 2*dX*dC + X*d2C
                for (int i = 0; i < M; i++) {
                    double q = xs * Cs[i];
                    Q.set(i, s, q);
                    Qtot[currentpop][i] += q;
                    Qg[i][gs] += q;
                    for (int p = 0; p < P; p++) {
                        double dQ = dX[p] * Cs[i] + xs * dCs[i][p];
                        double d2Q = d2X[p] * Cs[i] + 2 * dX[p] * dCs[i][p] + xs * d2Cs[i][p];
                        D1[currentpop][i][p] += dQ;
                        D2[currentpop][i][p] += d2Q;
                        D1Qg[i][gs][p] += dQ;
                        D2Qg[i][gs][p] += d2Q;
                    }
                }
                s++;
            }

            // ---- odometer advance ---------------------------------------------
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

        // moments at the full population: Qg, D1Qg and D2Qg were overwritten on every
        // population sweep, so they now hold the values at N.
        for (int i = 0; i < M; i++) {
            for (int g = 0; g < G; g++) {
                m.set(i, g, Qg[i][g]);
                for (int j = 0; j < M; j++) {
                    for (int g2 = 0; g2 < G; g2++) {
                        dm[i][g].set(j, g2, D1Qg[i][g][pidx[j][g2]]);
                    }
                }
                d2m.set(i, g, D2Qg[i][g][pidx[i][g]]);
            }
        }

        return pack(X, Q, U, C, m, dm, d2m, G);
    }

    /**
     * The per-station totals, i.e. the default single group holding every class.
     *
     * @param L  service demand matrix (M x R)
     * @param N  population vector (1 x R)
     * @param Z  think time vector (1 x R), null or empty for zeros
     * @param mi station server multiplicity (1 x M), null for ones
     * @return the base measures together with the exact moments of the total queue
     *         lengths
     */
    public static Ret.pfqnSensMom pfqn_sens_mom(Matrix L, Matrix N, Matrix Z, Matrix mi) {
        return pfqn_sens_mom(L, N, Z, mi, null);
    }

    public static Ret.pfqnSensMom pfqn_sens_mom(Matrix L, Matrix N, Matrix Z) {
        return pfqn_sens_mom(L, N, Z, null, null);
    }

    public static Ret.pfqnSensMom pfqn_sens_mom(Matrix L, Matrix N) {
        return pfqn_sens_mom(L, N, null, null, null);
    }

    /**
     * The class-to-group map that puts one class per group, {@code 1..R}, which
     * yields genuinely per-class moments (Akyildiz-Strelen Theorem 1 with T = {r}).
     *
     * @param R the number of classes
     * @return a (1 x R) group map
     */
    public static Matrix perClassGroups(int R) {
        Matrix g = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            g.set(0, r, r + 1);
        }
        return g;
    }

    // =====================================================================
    private static Ret.pfqnSensMom pack(Matrix X, Matrix Q, Matrix U, Matrix C,
                                        Matrix m, Matrix[][] dm, Matrix d2m, int G) {
        int M = Q.getNumRows();

        // see _kb/03-api-layer.md for rationale
        double covAsym = 0.0;
        Matrix[][] CovG = new Matrix[M][G];
        for (int i = 0; i < M; i++) {
            for (int g = 0; g < G; g++) {
                CovG[i][g] = new Matrix(M, G);
            }
        }
        for (int i = 0; i < M; i++) {
            for (int g = 0; g < G; g++) {
                for (int j = 0; j < M; j++) {
                    for (int g2 = 0; g2 < G; g2++) {
                        double a = dm[i][g].get(j, g2);
                        double b = dm[j][g2].get(i, g);
                        covAsym = Math.max(covAsym, Math.abs(a - b));
                        CovG[i][g].set(j, g2, (a + b) / 2.0);
                    }
                }
            }
        }

        // the collapsed single-group views; the group index carries no information
        // when G == 1, and an (M x 1 x M x 1) array would only be awkward to index
        Matrix Cov = null;
        Matrix dmFlat = null;
        if (G == 1) {
            Cov = new Matrix(M, M);
            dmFlat = new Matrix(M, M);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    Cov.set(i, j, CovG[i][0].get(j, 0));
                    dmFlat.set(i, j, dm[i][0].get(j, 0));
                }
            }
        }

        Matrix Var = new Matrix(M, G);
        Matrix M2 = new Matrix(M, G);
        Matrix M3 = new Matrix(M, G);
        Matrix Skew = new Matrix(M, G);
        for (int i = 0; i < M; i++) {
            for (int g = 0; g < G; g++) {
                double mig = m.get(i, g);
                double d1 = dm[i][g].get(i, g);
                Var.set(i, g, d1);                                                // (3.2)
                M2.set(i, g, d1 + mig * mig);                                     // (3.2)
                M3.set(i, g, d2m.get(i, g) + (1 + 3 * mig) * d1 + mig * mig * mig); // (3.2)
                // third central moment mu3 = E[Q^3] - 3 m E[Q^2] + 2 m^3
                double mu3 = M3.get(i, g) - 3 * mig * M2.get(i, g) + 2 * mig * mig * mig;
                if (Var.get(i, g) > 0) {
                    Skew.set(i, g, mu3 / Math.pow(Var.get(i, g), 1.5));
                } else {
                    Skew.set(i, g, Double.NaN);
                }
            }
        }
        return new Ret.pfqnSensMom(X, Q, U, C, G, m, dmFlat, dm, d2m, Var, Cov, CovG,
                M2, M3, Skew, covAsym);
    }
}
