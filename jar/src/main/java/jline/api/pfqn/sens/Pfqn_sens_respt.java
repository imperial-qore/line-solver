/**
 * Exact moments of the sojourn time of a job at FCFS multiserver centers of a
 * closed product-form queueing network.
 *
 * <p>Mirrors the MATLAB reference {@code pfqn_sens_respt.m}.</p>
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.sens;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_sens_respt {
    private Pfqn_sens_respt() {}

    /**
     * Exact raw moments E[W_(i,l)^t], t = 1..tmax, of the sojourn time of a class-l
     * job at an FCFS b-server center i of a closed product-form queueing network,
     * together with the variance of that sojourn time.
     *
     * <p>This is Theorem 4.1 of the reference. Its mechanism is the arrival theorem
     * of Lavenberg-Reiser and Sevcik-Mitrani: a class-l job arriving at center i
     * finds j jobs already there with probability p_i(j, N - 1_l). Conditioning the
     * sojourn time on j and inverting the Laplace transform of the conditional
     * density gives</p>
     *
     * <pre>
     *   E[W_(i,l)^t] = t!/mu^t + sum_{tau=0..t} a_(t,tau)(0) E[Qt_i^tau]
     *                  - sum_{j=0..b-1} p_i(j,N-1_l) sum_{tau=0..t} a_(t,tau)(0) j^tau
     * </pre>
     *
     * <p>where mu = 1/S(i) is the rate of each of the b servers, Qt_i is the total
     * queue length at center i at population N - 1_l (so its moments are those of
     * {@link Pfqn_sens_mom} evaluated one job down in class l), and the coefficients
     * a_(t,tau)(0) depend only on b and mu, not on the network (Remark 4.3 of the
     * reference). The moments E[Qt_i^tau] up to tau = 3 need the second derivative
     * of the MVA recursion, so this routine carries a second-order forward-mode pass
     * exactly as {@link Pfqn_sens_mom} does, but over the b-server recursion
     * (4.1)-(4.2) rather than the single-server one.</p>
     *
     * <p>For b = 1 the coefficients a_(t,0)(0) vanish identically and the double-sum
     * correction disappears, so no marginal probabilities are needed (Remark 4.2 of
     * the reference); the routine still evaluates the general expression, which
     * reduces to that case on its own.</p>
     *
     * <p>Only FCFS centers are covered. The reference is explicit that the
     * sojourn-time distribution at PS and LCFS centers is in general not known, so
     * no analogue exists there. FCFS in a BCMP network further requires the service
     * time to be exponential and class-independent, which is why this routine takes
     * a per-station service time S(i) and a separate visit-ratio matrix V rather
     * than a demand matrix: the sojourn time is per visit, so the per-visit rate
     * mu = 1/S(i) must be known and cannot be recovered from the demand
     * L(i,l) = S(i)*V(i,l) alone.</p>
     *
     * <p>Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks and
     * its Linearizer", Performance Evaluation 11:127-142, 1990, Theorem 4.1 with
     * equations (4.1)-(4.5) and Remarks 4.2-4.3.</p>
     *
     * <p>Restricted to closed populations. Load-dependent rates are not covered
     * here; the b-server dependence is the only state dependence, and it is carried
     * exactly by (4.1)-(4.2).</p>
     *
     * @param S    service time at each station (M x 1), common to all classes
     * @param V    visit ratio matrix (M x R); the demand is L(i,r) = S(i)*V(i,r)
     * @param N    population vector (1 x R)
     * @param Z    think time vector (1 x R), null or empty for zeros
     * @param b    number of servers at each station (M x 1), null for ones
     * @param tmax highest sojourn-time moment to return, 1..3; the coefficients
     *             a_(t,tau)(0) are tabulated in the reference up to t = 3
     * @return the base measures together with the exact sojourn-time moments
     */
    public static Ret.pfqnSensRespt pfqn_sens_respt(Matrix S, Matrix V, Matrix N, Matrix Z,
                                                    Matrix b, int tmax) {
        int M = V.getNumRows();
        int R = V.getNumCols();

        if (S.getNumCols() > 1) {
            S = S.transpose();
        }
        N = N.copy().ceil();
        if (N.getNumRows() > 1) {
            N = N.transpose();
        }
        if (Z == null || Z.isEmpty()) {
            Z = new Matrix(1, R);
        } else if (Z.getNumRows() > 1) {
            Z = Z.transpose();
        }
        int[] bs = new int[M];
        if (b == null || b.isEmpty()) {
            for (int i = 0; i < M; i++) {
                bs[i] = 1;
            }
        } else {
            Matrix bb = b.getNumCols() > 1 ? b.transpose() : b;
            for (int i = 0; i < M; i++) {
                bs[i] = (int) Math.round(bb.get(i, 0));
            }
        }
        if (tmax < 1 || tmax > 3) {
            throw new RuntimeException("pfqn_sens_respt: tmax must be 1, 2 or 3: the coefficients "
                    + "a_{t,tau}(0) are tabulated in the reference up to order three.");
        }
        if (N.length() != R) {
            throw new RuntimeException("pfqn_sens_respt: visit matrix and population vector have "
                    + "different number of classes");
        }
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(N.get(0, r))) {
                throw new RuntimeException("pfqn_sens_respt: requires a closed population");
            }
        }
        for (int i = 0; i < M; i++) {
            if (bs[i] < 1) {
                throw new RuntimeException("pfqn_sens_respt: the number of servers must be at "
                        + "least one at every station");
            }
            if (S.get(i, 0) <= 0) {
                throw new RuntimeException("pfqn_sens_respt: every FCFS station must have a "
                        + "strictly positive service time");
            }
        }

        double[][] rho = new double[M][R];
        double[] mu = new double[M];
        int bmax = 1;
        for (int i = 0; i < M; i++) {
            for (int l = 0; l < R; l++) {
                rho[i][l] = S.get(i, 0) * V.get(i, l);
            }
            mu[i] = 1.0 / S.get(i, 0);
            bmax = Math.max(bmax, bs[i]);
        }

        Matrix X = new Matrix(1, R);
        Matrix Q = new Matrix(M, R);
        Matrix U = new Matrix(M, R);
        Matrix m = new Matrix(M, 1);
        Matrix W = new Matrix(M, R);
        Matrix[] WM = new Matrix[tmax];
        for (int t = 0; t < tmax; t++) {
            WM[t] = new Matrix(M, R);
        }
        Matrix Wresid = new Matrix(M, R);
        Matrix pN = new Matrix(M, bmax);

        if (!N.any()) {
            return pack(X, Q, U, m, new Matrix(M, 1), pN, W, WM, Wresid, tmax);
        }

        // ---- population lattice ---------------------------------------------
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
        double[][] Mrow = new double[TP][M];
        double[][][] D1m = new double[TP][M][M];
        double[][][] D2m = new double[TP][M][M];
        double[][][] Prow = new double[TP][M][bmax];
        double[][][][] D1p = new double[TP][M][bmax][M];
        double[][][][] D2p = new double[TP][M][bmax][M];
        for (int i = 0; i < M; i++) {
            Prow[0][i][0] = 1.0;   // empty population: every station holds zero jobs
        }

        int currentpop = 1;
        Matrix n = new Matrix(1, R);
        n.set(0, firstNonEmpty, 1);
        int[] rows = new int[R];

        double[][] wv = new double[M][R];
        double[][][] d1w = new double[M][R][M];
        double[][][] d2w = new double[M][R][M];
        double[] lam = new double[R];
        double[][] d1lam = new double[R][M];
        double[][] d2lam = new double[R][M];
        double[] d1ui = new double[M];
        double[] d2ui = new double[M];

        while (ctr > 0) {
            int hnvec = currentpop;

            // ---- residence times, eq. (4.1), and their first two derivatives ----
            for (int i = 0; i < M; i++) {
                for (int s = 0; s < R; s++) {
                    wv[i][s] = 0.0;
                    for (int h = 0; h < M; h++) {
                        d1w[i][s][h] = 0.0;
                        d2w[i][s][h] = 0.0;
                    }
                }
            }
            for (int s = 0; s < R; s++) {
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
                rows[s] = pos;
                if (n.get(0, s) == 0) {
                    continue;   // w and every derivative stay zero, as does X(s)
                }
                for (int i = 0; i < M; i++) {
                    // bracket = 1 + m_i(n-e_s) + sum_{j=0}^{b_i-2} (b_i-1-j) p_i(j,n-e_s)
                    double brk = 1.0 + Mrow[pos][i];
                    for (int j = 0; j <= bs[i] - 2; j++) {
                        brk += (bs[i] - 1 - j) * Prow[pos][i][j];
                    }
                    wv[i][s] = (rho[i][s] / bs[i]) * brk;
                    for (int h = 0; h < M; h++) {
                        double dbrk = D1m[pos][i][h];
                        double d2brk = D2m[pos][i][h];
                        for (int j = 0; j <= bs[i] - 2; j++) {
                            dbrk += (bs[i] - 1 - j) * D1p[pos][i][j][h];
                            d2brk += (bs[i] - 1 - j) * D2p[pos][i][j][h];
                        }
                        // w = y_i * (rho/b) * brk
                        if (i == h) {
                            d1w[i][s][h] = (rho[i][s] / bs[i]) * (brk + dbrk);
                            d2w[i][s][h] = (rho[i][s] / bs[i]) * (2 * dbrk + d2brk);
                        } else {
                            d1w[i][s][h] = (rho[i][s] / bs[i]) * dbrk;
                            d2w[i][s][h] = (rho[i][s] / bs[i]) * d2brk;
                        }
                    }
                }
            }

            // ---- throughputs and their derivatives ------------------------------
            for (int s = 0; s < R; s++) {
                lam[s] = 0.0;
                for (int h = 0; h < M; h++) {
                    d1lam[s][h] = 0.0;
                    d2lam[s][h] = 0.0;
                }
            }
            for (int s = 0; s < R; s++) {
                if (n.get(0, s) == 0) {
                    continue;
                }
                double sumw = 0.0;
                for (int i = 0; i < M; i++) {
                    sumw += wv[i][s];
                }
                double den = Z.get(0, s) + sumw;
                double ns = n.get(0, s);
                lam[s] = ns / den;
                for (int h = 0; h < M; h++) {
                    double dden = 0.0;
                    double d2den = 0.0;
                    for (int i = 0; i < M; i++) {
                        dden += d1w[i][s][h];
                        d2den += d2w[i][s][h];
                    }
                    d1lam[s][h] = -ns * dden / (den * den);
                    d2lam[s][h] = -ns * d2den / (den * den) + 2 * ns * dden * dden / (den * den * den);
                }
            }

            // ---- mean queue lengths ---------------------------------------------
            for (int i = 0; i < M; i++) {
                double acc = 0.0;
                for (int s = 0; s < R; s++) {
                    if (n.get(0, s) == 0) {
                        continue;
                    }
                    acc += lam[s] * wv[i][s];
                }
                Mrow[hnvec][i] = acc;
                for (int h = 0; h < M; h++) {
                    double d1acc = 0.0;
                    double d2acc = 0.0;
                    for (int s = 0; s < R; s++) {
                        if (n.get(0, s) == 0) {
                            continue;
                        }
                        d1acc += d1lam[s][h] * wv[i][s] + lam[s] * d1w[i][s][h];
                        d2acc += d2lam[s][h] * wv[i][s] + 2 * d1lam[s][h] * d1w[i][s][h]
                                + lam[s] * d2w[i][s][h];
                    }
                    D1m[hnvec][i][h] = d1acc;
                    D2m[hnvec][i][h] = d2acc;
                }
            }

            // ---- marginal probabilities, eq. (4.2), and their derivatives --------
            double nc = 0.0;
            for (int s = 0; s < R; s++) {
                nc += n.get(0, s);
            }
            for (int i = 0; i < M; i++) {
                // p_i(j,n) = (1/j) sum_l lam(l) * rho_i(l)*y_i * p_i(j-1, n-e_l)
                for (int j = 1; j <= bs[i] - 1; j++) {
                    if (j > nc) {
                        Prow[hnvec][i][j] = 0.0;
                        continue;
                    }
                    double acc = 0.0;
                    for (int l = 0; l < R; l++) {
                        if (n.get(0, l) == 0) {
                            continue;
                        }
                        acc += lam[l] * rho[i][l] * Prow[rows[l]][i][j - 1];
                    }
                    Prow[hnvec][i][j] = acc / j;
                    for (int h = 0; h < M; h++) {
                        double d1acc = 0.0;
                        double d2acc = 0.0;
                        for (int l = 0; l < R; l++) {
                            if (n.get(0, l) == 0) {
                                continue;
                            }
                            double pprev = Prow[rows[l]][i][j - 1];
                            double d1prev = D1p[rows[l]][i][j - 1][h];
                            double d2prev = D2p[rows[l]][i][j - 1][h];
                            // g = rho * u * v with u = lam, v = y_i * pprev
                            double v1;
                            double v2;
                            if (i == h) {
                                v1 = pprev + d1prev;
                                v2 = 2 * d1prev + d2prev;
                            } else {
                                v1 = d1prev;
                                v2 = d2prev;
                            }
                            d1acc += rho[i][l] * (d1lam[l][h] * pprev + lam[l] * v1);
                            d2acc += rho[i][l] * (d2lam[l][h] * pprev + 2 * d1lam[l][h] * v1
                                    + lam[l] * v2);
                        }
                        D1p[hnvec][i][j][h] = d1acc / j;
                        D2p[hnvec][i][j][h] = d2acc / j;
                    }
                }
                // u_i = sum_l lam(l) * rho_i(l)*y_i  (mean number of busy servers)
                double ui = 0.0;
                for (int h = 0; h < M; h++) {
                    d1ui[h] = 0.0;
                    d2ui[h] = 0.0;
                }
                for (int l = 0; l < R; l++) {
                    if (n.get(0, l) == 0) {
                        continue;
                    }
                    ui += lam[l] * rho[i][l];
                    for (int h = 0; h < M; h++) {
                        if (i == h) {
                            d1ui[h] += rho[i][l] * (d1lam[l][h] + lam[l]);
                            d2ui[h] += rho[i][l] * (d2lam[l][h] + 2 * d1lam[l][h]);
                        } else {
                            d1ui[h] += rho[i][l] * d1lam[l][h];
                            d2ui[h] += rho[i][l] * d2lam[l][h];
                        }
                    }
                }
                // p_i(0,n) = 1 - (1/b)(u_i + sum_{j=1}^{b-1} (b-j) p_i(j,n))
                double acc0 = ui;
                for (int j = 1; j <= bs[i] - 1; j++) {
                    acc0 += (bs[i] - j) * Prow[hnvec][i][j];
                }
                Prow[hnvec][i][0] = 1.0 - acc0 / bs[i];
                for (int h = 0; h < M; h++) {
                    double d1acc0 = d1ui[h];
                    double d2acc0 = d2ui[h];
                    for (int j = 1; j <= bs[i] - 1; j++) {
                        d1acc0 += (bs[i] - j) * D1p[hnvec][i][j][h];
                        d2acc0 += (bs[i] - j) * D2p[hnvec][i][j][h];
                    }
                    D1p[hnvec][i][0][h] = -d1acc0 / bs[i];
                    D2p[hnvec][i][0][h] = -d2acc0 / bs[i];
                }
            }

            // keep the measures of the last (full) population
            for (int s = 0; s < R; s++) {
                X.set(0, s, lam[s]);
            }
            for (int i = 0; i < M; i++) {
                for (int s = 0; s < R; s++) {
                    Wresid.set(i, s, wv[i][s]);
                    Q.set(i, s, lam[s] * wv[i][s]);
                    U.set(i, s, lam[s] * rho[i][s]);
                }
            }

            // ---- odometer advance ------------------------------------------------
            int s = R - 1;
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

        int lastrow = currentpop;
        for (int i = 0; i < M; i++) {
            m.set(i, 0, Mrow[lastrow][i]);
            for (int j = 0; j < bmax; j++) {
                pN.set(i, j, Prow[lastrow][i][j]);
            }
        }
        Matrix Var = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            Var.set(i, 0, D1m[lastrow][i][i]);
        }

        // ---- sojourn-time moments, eq. (4.5) -------------------------------------
        // index of N - e_l on the lattice
        int[] rowsN = new int[R];
        for (int l = 0; l < R; l++) {
            if (N.get(0, l) > 0) {
                Matrix nn = N.copy();
                nn.set(0, l, nn.get(0, l) - 1);
                int pos = (int) nn.get(0, R - 1);
                for (int w = 0; w < R - 1; w++) {
                    pos = (int) (pos + nn.get(0, w) * prods.get(0, w));
                }
                rowsN[l] = pos;
            }
        }

        for (int i = 0; i < M; i++) {
            for (int l = 0; l < R; l++) {
                if (N.get(0, l) == 0 || V.get(i, l) <= 0) {
                    continue;
                }
                int rl = rowsN[l];
                // moments of the queue length seen by an arriving class-l job, i.e.
                // of the total queue at station i at population N - e_l, from (3.2)
                double mt = Mrow[rl][i];
                double d1t = D1m[rl][i][i];
                double d2t = D2m[rl][i][i];
                double[] EQ = new double[4];        // EQ[tau] = E[Qt_i^tau]
                EQ[0] = 1.0;
                EQ[1] = mt;
                EQ[2] = d1t + mt * mt;
                EQ[3] = d2t + (1 + 3 * mt) * d1t + mt * mt * mt;
                double[][] acoef = strelenA(bs[i], mu[i], tmax);
                for (int t = 1; t <= tmax; t++) {
                    double val = factorial(t) / Math.pow(mu[i], t);
                    for (int tau = 0; tau <= t; tau++) {
                        val += acoef[t - 1][tau] * EQ[tau];
                    }
                    // correction over the states in which a server is idle
                    for (int j = 0; j <= bs[i] - 1; j++) {
                        double inner = 0.0;
                        for (int tau = 0; tau <= t; tau++) {
                            inner += acoef[t - 1][tau] * jpow(j, tau);
                        }
                        val -= Prow[rl][i][j] * inner;
                    }
                    WM[t - 1].set(i, l, val);
                }
                W.set(i, l, WM[0].get(i, l));
            }
        }

        return pack(X, Q, U, m, Var, pN, W, WM, Wresid, tmax);
    }

    public static Ret.pfqnSensRespt pfqn_sens_respt(Matrix S, Matrix V, Matrix N, Matrix Z,
                                                    Matrix b) {
        return pfqn_sens_respt(S, V, N, Z, b, 3);
    }

    public static Ret.pfqnSensRespt pfqn_sens_respt(Matrix S, Matrix V, Matrix N, Matrix Z) {
        return pfqn_sens_respt(S, V, N, Z, null, 3);
    }

    public static Ret.pfqnSensRespt pfqn_sens_respt(Matrix S, Matrix V, Matrix N) {
        return pfqn_sens_respt(S, V, N, null, null, 3);
    }

    // =========================================================================
    /** j^tau with the convention j^0 = 1, so that 0^0 = 1 as the reference states. */
    private static double jpow(int j, int tau) {
        if (tau == 0) {
            return 1.0;
        }
        return Math.pow(j, tau);
    }

    private static double factorial(int t) {
        double v = 1.0;
        for (int i = 2; i <= t; i++) {
            v *= i;
        }
        return v;
    }

    // =========================================================================
    /**
     * Coefficients a_{t,tau}(0) of Remark 4.3 of the reference. They depend only on
     * the number of servers b and on the per-server rate mu, not on the network.
     */
    private static double[][] strelenA(int b, double mu, int tmax) {
        double[][] a = new double[3][4];
        a[0][0] = (1 - b) / (b * mu);
        a[0][1] = 1.0 / (b * mu);
        if (tmax >= 2) {
            a[1][0] = (2 - b - b * b) / (b * b * mu * mu);
            a[1][1] = 3.0 / (b * b * mu * mu);
            a[1][2] = 1.0 / (b * b * mu * mu);
        }
        if (tmax >= 3) {
            a[2][0] = (6 - 5 * b + 3 * b * b - 4 * b * b * b) / (b * b * b * mu * mu * mu);
            a[2][1] = (11 - 3 * b + 3 * b * b) / (b * b * b * mu * mu * mu);
            a[2][2] = 6.0 / (b * b * b * mu * mu * mu);
            a[2][3] = 1.0 / (b * b * b * mu * mu * mu);
        }
        return a;
    }

    // =========================================================================
    private static Ret.pfqnSensRespt pack(Matrix X, Matrix Q, Matrix U, Matrix m, Matrix Var,
                                          Matrix p, Matrix W, Matrix[] WM, Matrix Wresid,
                                          int tmax) {
        int M = Q.getNumRows();
        int R = Q.getNumCols();
        Matrix WVar = new Matrix(M, R);
        Matrix WSkew = new Matrix(M, R);
        if (tmax >= 2) {
            for (int i = 0; i < M; i++) {
                for (int l = 0; l < R; l++) {
                    WVar.set(i, l, WM[1].get(i, l) - WM[0].get(i, l) * WM[0].get(i, l));
                }
            }
        }
        if (tmax >= 3) {
            for (int i = 0; i < M; i++) {
                for (int l = 0; l < R; l++) {
                    double w1 = WM[0].get(i, l);
                    double mu3 = WM[2].get(i, l) - 3 * w1 * WM[1].get(i, l) + 2 * w1 * w1 * w1;
                    if (WVar.get(i, l) > 0) {
                        WSkew.set(i, l, mu3 / Math.pow(WVar.get(i, l), 1.5));
                    } else {
                        WSkew.set(i, l, Double.NaN);
                    }
                }
            }
        }
        return new Ret.pfqnSensRespt(X, Q, U, W, WM, WVar, WSkew, m, Var, p, Wresid);
    }
}
