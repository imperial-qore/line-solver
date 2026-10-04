/**
 * Approximate higher moments of the queue lengths of a closed product-form
 * queueing network, by differentiating the Linearizer fixed point. Polynomial in
 * the population, unlike the exact {@link jline.api.pfqn.sens.Pfqn_sens_mom}.
 *
 * <p>Mirrors the MATLAB reference {@code pfqn_sens_linearizer.m}.</p>
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.sens;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_sens_linearizer {
    private Pfqn_sens_linearizer() {}

    /**
     * Approximate moments E[Q_i], Var[Q_i], Cov[Q_i,Q_j], E[Q_i^2] and E[Q_i^3] of
     * the per-station total queue lengths of a closed product-form queueing
     * network, by the LINEARIZER-2 / LINEARIZER-3 algorithms of the reference
     * (Section 5).
     *
     * <p>Motivation. The exact moment analysis of {@link Pfqn_sens_mom} evaluates
     * the MVA recursion on the whole population lattice, so it costs O(prod(N+1))
     * and is unusable once the populations are large. The Linearizer replaces that
     * lattice by a fixed point over a handful of populations, and the reference
     * observes that the same trick applies to the derivatives: differentiate the
     * Linearizer equations, append the differentiated equations to the originals,
     * and iterate all of them together. This routine does that, carrying both the
     * first and the second derivative, so it returns everything (3.2) needs,
     * including the third moment. Carrying only the first derivative is the
     * reference's LINEARIZER-2; carrying the second as well is its LINEARIZER-3.</p>
     *
     * <p>The approximation. CORE (equations (5.1)-(5.2)) estimates the queue lengths
     * at population n - 1_l from those at n by</p>
     *
     * <pre>
     *   v_i(l)          = m_i^(n)(l) / n(l)
     *   m_i^(n-1_l')(l) = (n - 1_l')_l * ( v_i(l) + delta_i(l',l) )
     * </pre>
     *
     * <p>and substitutes them into the exact MVA equations. Setting the delta terms
     * to zero gives Bard-Schweitzer; Linearizer instead estimates them from (5.3),
     * delta_i^(N)(l',l) = v_i^(N-1_l')(l) - v_i^(N)(l), by running CORE at each of
     * the N - 1_l populations, and holds them fixed across populations (the
     * heuristic (5.4)). Differentiating (5.1)-(5.3) gives (5.5)-(5.8), which are
     * carried alongside.</p>
     *
     * <p>Accuracy. The reference reports, over 51 networks including 34 stress
     * cases, relative errors below 2.1% on E[Q], 4.1% on E[Q^2] and 6.2% on E[Q^3].
     * {@code PfqnSensLinearizerTest} measures the error against the exact
     * {@link Pfqn_sens_mom} on models small enough for both, and asserts bands of
     * that order rather than machine precision: this routine is an approximation and
     * is expected to disagree with the exact answer.</p>
     *
     * <p>Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks and
     * its Linearizer", Performance Evaluation 11:127-142, 1990, Section 5, equations
     * (5.1)-(5.8), the CORE-2 and LINEARIZER-2 algorithms. The delta definition
     * follows the standard Chandy-Neuse Linearizer, v_i^(N-1_l')(l) - v_i^(N)(l),
     * which is what LINE's {@code pfqn_linearizer} implements.</p>
     *
     * <p>Single-server stations plus an optional delay Z, matching
     * {@code pfqn_linearizer}. The exact counterpart is {@link Pfqn_sens_mom}; the
     * per-class exact second moments are in {@link Pfqn_sens_mva}.</p>
     *
     * @param L       service demand matrix (M x R), L(i,r) = visits_ir / rate_ir
     * @param N       population vector (1 x R)
     * @param Z       think time vector (1 x R), null or empty for zeros
     * @param tol     convergence tolerance of the CORE fixed point on the mean queue
     *                lengths; pass a non-positive value for the test of the
     *                reference, 1/(4000+16*sum(n)), which is also applied to the
     *                variances at 1e-3
     * @param maxiter maximum CORE iterations
     * @return the approximate base measures together with the approximate moments
     */
    public static Ret.pfqnSensLinearizer pfqn_sens_linearizer(Matrix L, Matrix N, Matrix Z,
                                                              double tol, int maxiter) {
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
        if (N.length() != R) {
            throw new RuntimeException("pfqn_sens_linearizer: demand matrix and population vector "
                    + "have different number of classes");
        }
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(N.get(0, r))) {
                throw new RuntimeException("pfqn_sens_linearizer: requires a closed population");
            }
        }

        if (!N.any()) {
            return pack(new Matrix(1, R), new Matrix(M, R), new Matrix(M, R), new Matrix(M, R),
                    new Matrix(M, 1), new Matrix(M, M), new Matrix(M, 1), 0);
        }

        // ---- Linearizer state ----------------------------------------------
        // pops: index 0 = N, index 1+l = N - e_l
        int npops = 1 + R;
        double[][] pv = new double[npops][R];
        for (int r = 0; r < R; r++) {
            pv[0][r] = N.get(0, r);
        }
        for (int l = 0; l < R; l++) {
            for (int r = 0; r < R; r++) {
                pv[1 + l][r] = N.get(0, r);
            }
            if (N.get(0, l) > 0) {
                pv[1 + l][l] = N.get(0, l) - 1;
            }
        }

        // mE[p](i,l) = estimate of m_i^{(pv[p])}(l), and its derivatives w.r.t. y_h
        double[][][] mE = new double[npops][M][R];
        double[][][][] dmE = new double[npops][M][R][M];
        double[][][][] d2mE = new double[npops][M][R][M];
        for (int p = 0; p < npops; p++) {
            for (int i = 0; i < M; i++) {
                for (int l = 0; l < R; l++) {
                    mE[p][i][l] = pv[p][l] / M;     // initialization of the reference
                }
            }
        }
        // delta[i][lp][l], indexed by the removed class lp and the class l
        double[][][] delta = new double[M][R][R];
        double[][][][] ddelta = new double[M][R][R][M];
        double[][][][] d2delta = new double[M][R][R][M];

        int totiter = 0;
        Matrix Xf = new Matrix(1, R);
        Matrix Wf = new Matrix(M, R);

        double[] nfull = new double[R];
        for (int r = 0; r < R; r++) {
            nfull[r] = N.get(0, r);
        }

        for (int outer = 1; outer <= 3; outer++) {
            // ---- Step 1: CORE at the full population --------------------------
            Core c0 = core2(L, Z, nfull, delta, ddelta, d2delta, mE[0], dmE[0], d2mE[0],
                    tol, maxiter);
            mE[0] = c0.m;
            dmE[0] = c0.dm;
            d2mE[0] = c0.d2m;
            totiter += c0.it;
            for (int r = 0; r < R; r++) {
                Xf.set(0, r, c0.lam[r]);
            }
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    Wf.set(i, r, c0.w[i][r]);
                }
            }

            if (outer < 3) {
                // ---- Step 2: CORE at each of the N - e_l populations ------------
                for (int l = 0; l < R; l++) {
                    if (N.get(0, l) == 0) {
                        continue;
                    }
                    Core cl = core2(L, Z, pv[1 + l], delta, ddelta, d2delta, mE[1 + l],
                            dmE[1 + l], d2mE[1 + l], tol, maxiter);
                    mE[1 + l] = cl.m;
                    dmE[1 + l] = cl.dm;
                    d2mE[1 + l] = cl.d2m;
                    totiter += cl.it;
                }

                // ---- Step 3: refresh delta from (5.1) and (5.3) ----------------
                Frac fN = fractions(mE[0], dmE[0], d2mE[0], pv[0], M, R);
                for (int lp = 0; lp < R; lp++) {
                    if (N.get(0, lp) == 0) {
                        continue;
                    }
                    Frac fL = fractions(mE[1 + lp], dmE[1 + lp], d2mE[1 + lp], pv[1 + lp], M, R);
                    for (int i = 0; i < M; i++) {
                        for (int l = 0; l < R; l++) {
                            delta[i][lp][l] = fL.v[i][l] - fN.v[i][l];
                            for (int h = 0; h < M; h++) {
                                ddelta[i][lp][l][h] = fL.dv[i][l][h] - fN.dv[i][l][h];
                                d2delta[i][lp][l][h] = fL.d2v[i][l][h] - fN.d2v[i][l][h];
                            }
                        }
                    }
                }
            }
        }

        // ---- final measures ------------------------------------------------------
        Matrix Qf = new Matrix(M, R);
        Matrix m = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            double tot = 0.0;
            for (int r = 0; r < R; r++) {
                Qf.set(i, r, mE[0][i][r]);
                tot += mE[0][i][r];
            }
            m.set(i, 0, tot);
        }
        Matrix dmf = new Matrix(M, M);
        Matrix d2mf = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            for (int h = 0; h < M; h++) {
                double acc = 0.0;
                for (int r = 0; r < R; r++) {
                    acc += dmE[0][i][r][h];
                }
                dmf.set(i, h, acc);
            }
            double acc2 = 0.0;
            for (int r = 0; r < R; r++) {
                acc2 += d2mE[0][i][r][i];
            }
            d2mf.set(i, 0, acc2);
        }
        Matrix U = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                U.set(i, r, Xf.get(0, r) * L.get(i, r));
            }
        }

        return pack(Xf, Qf, U, Wf, m, dmf, d2mf, totiter);
    }

    /** Uses the termination test of the reference and 200 CORE iterations. */
    public static Ret.pfqnSensLinearizer pfqn_sens_linearizer(Matrix L, Matrix N, Matrix Z) {
        return pfqn_sens_linearizer(L, N, Z, 0.0, 200);
    }

    public static Ret.pfqnSensLinearizer pfqn_sens_linearizer(Matrix L, Matrix N) {
        return pfqn_sens_linearizer(L, N, null, 0.0, 200);
    }

    // =========================================================================
    private static final class Frac {
        double[][] v;
        double[][][] dv;
        double[][][] d2v;
    }

    /** (5.1) and (5.5): v_i(l) = m_i(l)/n(l), and likewise for the derivatives. */
    private static Frac fractions(double[][] m, double[][][] dm, double[][][] d2m, double[] n,
                                  int M, int R) {
        Frac f = new Frac();
        f.v = new double[M][R];
        f.dv = new double[M][R][M];
        f.d2v = new double[M][R][M];
        for (int l = 0; l < R; l++) {
            if (n[l] <= 0) {
                continue;
            }
            for (int i = 0; i < M; i++) {
                f.v[i][l] = m[i][l] / n[l];
                for (int h = 0; h < M; h++) {
                    f.dv[i][l][h] = dm[i][l][h] / n[l];
                    f.d2v[i][l][h] = d2m[i][l][h] / n[l];
                }
            }
        }
        return f;
    }

    // =========================================================================
    private static final class Core {
        double[][] m;
        double[][][] dm;
        double[][][] d2m;
        double[] lam;
        double[][] w;
        int it;
    }

    /**
     * CORE-2 of the reference, extended to second derivatives. Iterates (5.1),
     * (5.2) and the MVA equations (1.4)-(1.5) together with their first and second
     * derivatives (5.5)-(5.6) and (3.4), until the mean queue lengths and the
     * variances both stop moving.
     */
    private static Core core2(Matrix L, Matrix Z, double[] n, double[][][] delta,
                              double[][][][] ddelta, double[][][][] d2delta,
                              double[][] m, double[][][] dm, double[][][] d2m,
                              double tol, int maxiter) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double nc = 0.0;
        for (int r = 0; r < R; r++) {
            nc += n[r];
        }
        double tolm;
        if (tol <= 0) {
            tolm = 1.0 / (4000 + 16 * nc);     // termination test of the reference
        } else {
            tolm = tol;
        }
        double tolv = 1e-3;

        double[] lam = new double[R];
        double[][] w = new double[M][R];
        double[] varprev = new double[M];
        int it = 0;
        for (int iter = 1; iter <= maxiter; iter++) {
            it = iter;
            double[][] mprev = new double[M][R];
            for (int i = 0; i < M; i++) {
                System.arraycopy(m[i], 0, mprev[i], 0, R);
            }

            // ---- (5.1)-(5.2) and (5.5)-(5.6): queue lengths one job down --------
            // mtot(i,l) = sum_l2 m_i^{(n-e_l)}(l2)
            Frac f = fractions(m, dm, d2m, n, M, R);
            double[][] mtot = new double[M][R];
            double[][][] dmtot = new double[M][R][M];
            double[][][] d2mtot = new double[M][R][M];
            for (int l = 0; l < R; l++) {
                if (n[l] <= 0) {
                    continue;
                }
                for (int i = 0; i < M; i++) {
                    double acc = 0.0;
                    double[] dacc = new double[M];
                    double[] d2acc = new double[M];
                    for (int l2 = 0; l2 < R; l2++) {
                        double cnt = n[l2] - (l2 == l ? 1 : 0);       // (n - 1_l)_{l2}
                        if (cnt <= 0) {
                            continue;
                        }
                        acc += cnt * (f.v[i][l2] + delta[i][l][l2]);
                        for (int h = 0; h < M; h++) {
                            dacc[h] += cnt * (f.dv[i][l2][h] + ddelta[i][l][l2][h]);
                            d2acc[h] += cnt * (f.d2v[i][l2][h] + d2delta[i][l][l2][h]);
                        }
                    }
                    mtot[i][l] = acc;
                    for (int h = 0; h < M; h++) {
                        dmtot[i][l][h] = dacc[h];
                        d2mtot[i][l][h] = d2acc[h];
                    }
                }
            }

            // ---- MVA (1.4)-(1.5) and its derivatives (3.4) ----------------------
            // w_i(l) = y_i * L(i,l) * (1 + mtot(i,l))
            w = new double[M][R];
            double[][][] dw = new double[M][R][M];
            double[][][] d2w = new double[M][R][M];
            for (int l = 0; l < R; l++) {
                if (n[l] <= 0) {
                    continue;
                }
                for (int i = 0; i < M; i++) {
                    double A = 1 + mtot[i][l];
                    w[i][l] = L.get(i, l) * A;
                    for (int h = 0; h < M; h++) {
                        double dA = dmtot[i][l][h];
                        double d2A = d2mtot[i][l][h];
                        if (i == h) {
                            dw[i][l][h] = L.get(i, l) * (A + dA);
                            d2w[i][l][h] = L.get(i, l) * (2 * dA + d2A);
                        } else {
                            dw[i][l][h] = L.get(i, l) * dA;
                            d2w[i][l][h] = L.get(i, l) * d2A;
                        }
                    }
                }
            }
            lam = new double[R];
            double[][] dlam = new double[R][M];
            double[][] d2lam = new double[R][M];
            for (int l = 0; l < R; l++) {
                if (n[l] <= 0) {
                    continue;
                }
                double sumw = 0.0;
                for (int i = 0; i < M; i++) {
                    sumw += w[i][l];
                }
                double den = Z.get(0, l) + sumw;
                lam[l] = n[l] / den;
                for (int h = 0; h < M; h++) {
                    double dden = 0.0;
                    double d2den = 0.0;
                    for (int i = 0; i < M; i++) {
                        dden += dw[i][l][h];
                        d2den += d2w[i][l][h];
                    }
                    dlam[l][h] = -n[l] * dden / (den * den);
                    d2lam[l][h] = -n[l] * d2den / (den * den)
                            + 2 * n[l] * dden * dden / (den * den * den);
                }
            }
            m = new double[M][R];
            dm = new double[M][R][M];
            d2m = new double[M][R][M];
            for (int l = 0; l < R; l++) {
                if (n[l] <= 0) {
                    continue;
                }
                for (int i = 0; i < M; i++) {
                    m[i][l] = lam[l] * w[i][l];
                    for (int h = 0; h < M; h++) {
                        dm[i][l][h] = dlam[l][h] * w[i][l] + lam[l] * dw[i][l][h];
                        d2m[i][l][h] = d2lam[l][h] * w[i][l] + 2 * dlam[l][h] * dw[i][l][h]
                                + lam[l] * d2w[i][l][h];
                    }
                }
            }

            // ---- termination test of the reference ------------------------------
            double dev = 0.0;
            for (int l = 0; l < R; l++) {
                if (n[l] <= 0) {
                    continue;
                }
                double mx = 0.0;
                for (int i = 0; i < M; i++) {
                    mx = Math.max(mx, Math.abs(m[i][l] - mprev[i][l]));
                }
                dev = Math.max(dev, mx / n[l]);
            }
            double[] varnow = new double[M];
            double sv = 0.0;
            for (int i = 0; i < M; i++) {
                double acc = 0.0;
                for (int l = 0; l < R; l++) {
                    acc += dm[i][l][i];
                }
                varnow[i] = acc;
                sv += acc;
            }
            double vdev = 0.0;
            if (sv > 0) {
                double mx = 0.0;
                for (int i = 0; i < M; i++) {
                    mx = Math.max(mx, Math.abs(varnow[i] - varprev[i]));
                }
                vdev = mx / sv;
            }
            varprev = varnow;
            if (dev <= tolm && vdev <= tolv) {
                break;
            }
        }

        Core c = new Core();
        c.m = m;
        c.dm = dm;
        c.d2m = d2m;
        c.lam = lam;
        c.w = w;
        c.it = it;
        return c;
    }

    // =========================================================================
    private static Ret.pfqnSensLinearizer pack(Matrix X, Matrix Q, Matrix U, Matrix W, Matrix m,
                                               Matrix dm, Matrix d2m, int it) {
        int M = Q.getNumRows();

        // see _kb/03-api-layer.md for rationale
        double covAsym = 0.0;
        Matrix Cov = new Matrix(M, M);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                covAsym = Math.max(covAsym, Math.abs(dm.get(i, j) - dm.get(j, i)));
                Cov.set(i, j, (dm.get(i, j) + dm.get(j, i)) / 2.0);
            }
        }

        Matrix Var = new Matrix(M, 1);
        Matrix M2 = new Matrix(M, 1);
        Matrix M3 = new Matrix(M, 1);
        Matrix Skew = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            double mi = m.get(i, 0);
            double dii = dm.get(i, i);
            Var.set(i, 0, dii);
            M2.set(i, 0, dii + mi * mi);
            M3.set(i, 0, d2m.get(i, 0) + (1 + 3 * mi) * dii + mi * mi * mi);
            double mu3 = M3.get(i, 0) - 3 * mi * M2.get(i, 0) + 2 * mi * mi * mi;
            if (Var.get(i, 0) > 0) {
                Skew.set(i, 0, mu3 / Math.pow(Var.get(i, 0), 1.5));
            } else {
                Skew.set(i, 0, Double.NaN);
            }
        }
        return new Ret.pfqnSensLinearizer(X, Q, U, W, m, dm, d2m, Var, Cov, M2, M3, Skew,
                covAsym, it);
    }
}
