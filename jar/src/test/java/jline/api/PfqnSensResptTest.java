/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.mva.Pfqn_mva;
import jline.api.pfqn.sens.Pfqn_sens_respt;
import jline.io.Ret;
import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates {@link Pfqn_sens_respt#pfqn_sens_respt}, the FCFS sojourn-time moment
 * analysis of Strelen (1990), Theorem 4.1, mirroring the MATLAB harness
 * {@code pfqn_sens_respt_validate.m}:
 *
 * <ul>
 *   <li>A. brute-force enumeration. This is the strongest check because it shares
 *       none of Theorem 4.1's algebra. The equilibrium product form is enumerated
 *       to get the exact arrival-theorem marginals p_i(j, N-1_l); the sojourn time
 *       conditioned on finding j jobs is known in closed form (Exp(mu) if j &lt; b,
 *       otherwise an Erlang(j-b+1, b*mu) queueing delay plus an Exp(mu) service),
 *       so its moments are formed directly and mixed over j. Neither the
 *       coefficients a_(t,tau)(0) nor the recursion (4.2) enter, so the agreement
 *       tests both;</li>
 *   <li>B. the published table of Example 3.4 (continued) of the reference, which
 *       prints E(W_i) and sigma^2_(W_i) for the Kobayashi model;</li>
 *   <li>C. the internal identity W(i,l) = w_i(l)/V(i,l): the t = 1 case of (4.5)
 *       must reproduce the MVA residence time divided by the visit ratio, which is
 *       a completely different expression;</li>
 *   <li>D. {@link Pfqn_mva} for the base measures in the single-server case;</li>
 *   <li>E. the single-job network, where an arriving job always finds an empty
 *       station, so W is exactly Exp(mu) and every moment is known in closed
 *       form.</li>
 * </ul>
 *
 * <p>Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks and
 * its Linearizer", Performance Evaluation 11:127-142, 1990.</p>
 */
public class PfqnSensResptTest {

    private static final double TOL_BRUTE = 1e-9;
    /** the paper prints 5 decimals on small sojourn times */
    private static final double TOL_PAPER = 5e-4;
    private static final double TOL_IDENT = 1e-10;
    private static final double TOL_MVA = 1e-10;
    private static final double TOL_EXP = 1e-12;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    private static Matrix col(double... v) {
        Matrix m = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) {
            m.set(i, 0, v[i]);
        }
        return m;
    }

    /** max_k |a_k - b_k| / max(1, |a_k|, |b_k|), matching the MATLAB relerr. */
    private static double relerr(double[] a, double[] b) {
        double e = 0.0;
        for (int k = 0; k < a.length; k++) {
            double scale = Math.max(1.0, Math.max(Math.abs(a[k]), Math.abs(b[k])));
            e = Math.max(e, Math.abs(a[k] - b[k]) / scale);
        }
        return e;
    }

    private static double relerr(Matrix a, Matrix b) {
        return relerr(flat(a), flat(b));
    }

    private static double[] flat(Matrix m) {
        double[] v = new double[m.getNumRows() * m.getNumCols()];
        int k = 0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                v[k++] = m.get(i, j);
            }
        }
        return v;
    }

    /** A/C/D. random models, single- and multi-server. */
    @Test
    public void sensResptMatchesBruteForceIdentityAndMva() {
        Random rng = new Random(5);
        double errBrute = 0.0;
        double errIdent = 0.0;
        double errMva = 0.0;
        int nBrute = 0;

        for (int trial = 1; trial <= 36; trial++) {
            int M = 1 + rng.nextInt(3);
            int R = 1 + rng.nextInt(2);
            Matrix S = new Matrix(M, 1);
            for (int i = 0; i < M; i++) {
                S.set(i, 0, 0.2 + rng.nextDouble());
            }
            Matrix V = new Matrix(M, R);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    V.set(i, r, 0.3 + rng.nextDouble());
                }
            }
            if (trial % 4 == 0 && M > 1) {
                V.set(0, 0, 0.0);      // a class that skips a station
            }
            Matrix N = new Matrix(1, R);
            for (int r = 0; r < R; r++) {
                N.set(0, r, 1 + rng.nextInt(3));
            }
            Matrix Z = new Matrix(1, R);
            if (trial % 2 == 0) {
                for (int r = 0; r < R; r++) {
                    Z.set(0, r, 0.3 + rng.nextDouble());
                }
            }
            Matrix b = Matrix.ones(M, 1);
            boolean single = true;
            if (trial % 3 == 0) {
                for (int i = 0; i < M; i++) {
                    b.set(i, 0, 1 + rng.nextInt(3));
                }
                single = false;
            }
            int[] bs = new int[M];
            for (int i = 0; i < M; i++) {
                bs[i] = (int) b.get(i, 0);
                single = single && bs[i] == 1;
            }

            Ret.pfqnSensRespt res = Pfqn_sens_respt.pfqn_sens_respt(S, V, N, Z, b, 3);

            // ---- C. W = w/V ------------------------------------------------
            for (int i = 0; i < M; i++) {
                for (int l = 0; l < R; l++) {
                    if (V.get(i, l) > 0 && N.get(0, l) > 0) {
                        errIdent = Math.max(errIdent, relerr(new double[]{res.W.get(i, l)},
                                new double[]{res.Wresid.get(i, l) / V.get(i, l)}));
                    }
                }
            }

            // ---- D. base measures, single-server only ----------------------
            if (single) {
                Matrix L = new Matrix(M, R);
                for (int i = 0; i < M; i++) {
                    for (int r = 0; r < R; r++) {
                        L.set(i, r, S.get(i, 0) * V.get(i, r));
                    }
                }
                Ret.pfqnMVA mva = Pfqn_mva.pfqn_mva(L, N, Z);
                errMva = Math.max(errMva, relerr(res.X, mva.X));
                errMva = Math.max(errMva, relerr(res.Q, mva.Q));
                errMva = Math.max(errMva, relerr(res.U, mva.U));
            }

            // ---- A. brute force --------------------------------------------
            double totpop = 1.0;
            for (int r = 0; r < R; r++) {
                totpop *= (N.get(0, r) + 1);
            }
            if (totpop <= 24 && M <= 3) {
                Brute bt = bruteRespt(S, V, N, Z, bs, 3);
                for (int t = 0; t < 3; t++) {
                    errBrute = Math.max(errBrute, relerr(res.WM[t], bt.WM[t]));
                }
                // .p is ragged: station i only defines j = 0..b_i-1, the range the
                // b-server recursion needs, and the rest of the row is zero padding
                // out to max(b). Comparing the padding against the true marginal
                // would be comparing against something the algorithm never claims to
                // compute.
                for (int i = 0; i < M; i++) {
                    double[] got = new double[bs[i]];
                    double[] ref = new double[bs[i]];
                    for (int j = 0; j < bs[i]; j++) {
                        got[j] = res.p.get(i, j);
                        ref[j] = bt.p.get(i, j);
                    }
                    errBrute = Math.max(errBrute, relerr(got, ref));
                }
                nBrute++;
            }
        }

        assertTrue(nBrute > 0, "no model was small enough for brute-force enumeration");
        assertTrue(errBrute <= TOL_BRUTE,
                "brute force, arrival theorem (" + nBrute + " models): " + errBrute);
        assertTrue(errIdent <= TOL_IDENT, "identity W(i,l) = w_i(l)/V(i,l): " + errIdent);
        assertTrue(errMva <= TOL_MVA, "pfqn_mva base measures (b=1): " + errMva);
    }

    /**
     * E. one job: the arriving job always finds the station empty, so W ~ Exp(mu)
     * and every moment is known in closed form.
     */
    @Test
    public void sensResptIsExponentialWithASingleJob() {
        Matrix S1 = col(0.4, 0.25);
        Matrix V1 = col(1, 2);
        Ret.pfqnSensRespt r1 = Pfqn_sens_respt.pfqn_sens_respt(S1, V1, row(1), row(0.7),
                col(1, 1), 3);
        double errExp = 0.0;
        for (int i = 0; i < 2; i++) {
            double mu = 1.0 / S1.get(i, 0);
            errExp = Math.max(errExp, relerr(new double[]{r1.WM[0].get(i, 0)},
                    new double[]{1 / mu}));
            errExp = Math.max(errExp, relerr(new double[]{r1.WM[1].get(i, 0)},
                    new double[]{2 / (mu * mu)}));
            errExp = Math.max(errExp, relerr(new double[]{r1.WM[2].get(i, 0)},
                    new double[]{6 / (mu * mu * mu)}));
            errExp = Math.max(errExp, relerr(new double[]{r1.WVar.get(i, 0)},
                    new double[]{1 / (mu * mu)}));
        }
        assertTrue(errExp <= TOL_EXP, "single job, W ~ Exp(mu) exactly: " + errExp);
    }

    /**
     * B. Example 3.4 (continued): Kobayashi central-server model, n = 3. The paper
     * prints E(W_i) and sigma^2_(W_i).
     */
    @Test
    public void sensResptMatchesStrelenExample34Sojourn() {
        Matrix xs = new Matrix(12, 1);
        Matrix es = new Matrix(12, 1);
        for (int i = 0; i < 9; i++) {
            xs.set(i, 0, 0.0215);
            es.set(i, 0, 9.333);
        }
        xs.set(9, 0, 0.104);
        es.set(9, 0, 10.5);
        xs.set(10, 0, 0.104);
        es.set(10, 0, 10.5);
        xs.set(11, 0, 0.019);
        es.set(11, 0, 105);

        Ret.pfqnSensRespt rk = Pfqn_sens_respt.pfqn_sens_respt(xs, es, row(3), row(0.0),
                Matrix.ones(12, 1), 3);
        double[] paperW = {0.02275, 0.14178, 0.03322};
        double[] paperWV = {0.00052, 0.01846, 0.00083};
        double[] gotW = {rk.W.get(0, 0), rk.W.get(9, 0), rk.W.get(11, 0)};
        double[] gotWV = {rk.WVar.get(0, 0), rk.WVar.get(9, 0), rk.WVar.get(11, 0)};
        double errPaper = Math.max(relerr(gotW, paperW), relerr(gotWV, paperWV));
        assertTrue(errPaper <= TOL_PAPER,
                "Strelen Example 3.4 published sojourn: " + errPaper
                        + " (paper E(W_12)=" + paperW[2] + " got " + gotW[2]
                        + " ; sigma2=" + paperWV[2] + " got " + gotWV[2] + ")");
    }

    // =========================================================================
    private static final class Brute {
        Matrix[] WM;
        Matrix p;
    }

    /**
     * Sojourn-time moments from first principles: enumerate the product form to get
     * the exact arrival-theorem marginals p_i(j, N-e_l), then mix the conditional
     * sojourn-time moments over j. Uses none of Theorem 4.1.
     */
    private static Brute bruteRespt(Matrix S, Matrix V, Matrix N, Matrix Z, int[] b, int tmax) {
        int M = V.getNumRows();
        int R = V.getNumCols();
        int bmax = 1;
        for (int i = 0; i < M; i++) {
            bmax = Math.max(bmax, b[i]);
        }
        Brute out = new Brute();
        out.WM = new Matrix[tmax];
        for (int t = 0; t < tmax; t++) {
            out.WM[t] = new Matrix(M, R);
        }
        for (int l = 0; l < R; l++) {
            if (N.get(0, l) == 0) {
                continue;
            }
            Matrix Nl = N.copy();
            Nl.set(0, l, Nl.get(0, l) - 1);
            double[][] pj = bruteMarginals(S, V, Nl, Z, b);   // at population N - e_l
            for (int i = 0; i < M; i++) {
                if (V.get(i, l) <= 0) {
                    continue;
                }
                double mu = 1.0 / S.get(i, 0);
                for (int t = 1; t <= tmax; t++) {
                    double acc = 0.0;
                    for (int j = 0; j < pj[i].length; j++) {
                        acc += pj[i][j] * condMoment(j, b[i], mu, t);
                    }
                    out.WM[t - 1].set(i, l, acc);
                }
            }
        }
        double[][] pAll = bruteMarginals(S, V, N, Z, b);
        out.p = new Matrix(M, bmax);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < bmax; j++) {
                out.p.set(i, j, j < pAll[i].length ? pAll[i][j] : 0.0);
            }
        }
        return out;
    }

    /**
     * E[(W|j)^t] where a job arriving to find j jobs at an FCFS b-server station
     * waits an Erlang(max(0,j-b+1), b*mu) and is then served for an Exp(mu).
     */
    private static double condMoment(int j, int b, double mu, int t) {
        int k = Math.max(0, j - b + 1);
        double v = 0.0;
        for (int s = 0; s <= t; s++) {
            // E[X^s] with X ~ Exp(mu)
            double EX = factorial(s) / Math.pow(mu, s);
            int p = t - s;
            // E[Y^p] with Y ~ Erlang(k, b*mu), and Y = 0 when k = 0
            double EY;
            if (k == 0) {
                EY = (p == 0) ? 1.0 : 0.0;
            } else {
                double theta = b * mu;
                EY = 1.0;
                for (int a = 0; a <= p - 1; a++) {
                    EY *= (k + a);
                }
                EY /= Math.pow(theta, p);
            }
            v += nchoosek(t, s) * EX * EY;
        }
        return v;
    }

    /**
     * P[Q_i = j] for every station, by enumerating the closed product form of a
     * network of FCFS b-server stations:
     * f_i(q_i) = q_i! prod_l (a(i,l)^q_il / q_il!) prod_{j=1}^{q_i} 1/min(j,b_i)
     * with a(i,l) = S(i)*V(i,l), plus the delay term for the think times.
     */
    private static double[][] bruteMarginals(Matrix S, Matrix V, Matrix N, Matrix Z, int[] b) {
        int M = V.getNumRows();
        int R = V.getNumCols();
        double[][] a = new double[M][R];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                a[i][r] = S.get(i, 0) * V.get(i, r);
            }
        }
        List<int[][]> states = enumerateStates(N, M, R);
        int K = states.size();
        double[] w = new double[K];
        int[][] tot = new int[K][M];
        for (int k = 0; k < K; k++) {
            int[][] nir = states.get(k);
            double lw = 0.0;
            boolean ok = true;
            for (int i = 0; i < M && ok; i++) {
                int ni = 0;
                for (int r = 0; r < R; r++) {
                    ni += nir[i][r];
                }
                lw += logFactorial(ni);
                for (int j = 1; j <= ni; j++) {
                    lw -= Math.log(Math.min(j, b[i]));
                }
                for (int r = 0; r < R; r++) {
                    if (nir[i][r] > 0) {
                        if (a[i][r] <= 0) {
                            ok = false;
                            break;
                        }
                        lw += nir[i][r] * Math.log(a[i][r]) - logFactorial(nir[i][r]);
                    }
                }
            }
            if (ok) {
                for (int r = 0; r < R; r++) {
                    int sum = 0;
                    for (int i = 0; i < M; i++) {
                        sum += nir[i][r];
                    }
                    int n0r = (int) N.get(0, r) - sum;
                    if (n0r > 0) {
                        if (Z.get(0, r) <= 0) {
                            ok = false;
                            break;
                        }
                        lw += n0r * Math.log(Z.get(0, r)) - logFactorial(n0r);
                    }
                }
            }
            w[k] = ok ? Math.exp(lw) : 0.0;
            for (int i = 0; i < M; i++) {
                int s = 0;
                for (int r = 0; r < R; r++) {
                    s += nir[i][r];
                }
                tot[k][i] = s;
            }
        }
        double sum = 0.0;
        for (int k = 0; k < K; k++) {
            sum += w[k];
        }
        for (int k = 0; k < K; k++) {
            w[k] /= sum;
        }
        int npop = 0;
        for (int r = 0; r < R; r++) {
            npop += (int) N.get(0, r);
        }
        int maxj = Math.max(npop, 1);
        double[][] pj = new double[M][1 + maxj];
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < M; i++) {
                pj[i][tot[k][i]] += w[k];
            }
        }
        return pj;
    }

    /** All allocations of N(r) class-r jobs over M stations (remainder in the delay). */
    private static List<int[][]> enumerateStates(Matrix N, int M, int R) {
        List<List<int[]>> per = new ArrayList<List<int[]>>();
        for (int r = 0; r < R; r++) {
            per.add(PfqnSensMvaTest.compositionsLeq((int) N.get(0, r), M));
        }
        List<int[][]> states = new ArrayList<int[][]>();
        int[] idx = new int[R];
        while (true) {
            int[][] row = new int[M][R];
            for (int r = 0; r < R; r++) {
                int[] alloc = per.get(r).get(idx[r]);
                for (int i = 0; i < M; i++) {
                    row[i][r] = alloc[i];
                }
            }
            states.add(row);
            int r = R - 1;
            while (r >= 0) {
                idx[r]++;
                if (idx[r] < per.get(r).size()) {
                    break;
                }
                idx[r] = 0;
                r--;
            }
            if (r < 0) {
                break;
            }
        }
        return states;
    }

    private static double factorial(int n) {
        double v = 1.0;
        for (int i = 2; i <= n; i++) {
            v *= i;
        }
        return v;
    }

    private static double nchoosek(int n, int k) {
        return factorial(n) / (factorial(k) * factorial(n - k));
    }

    private static double logFactorial(int n) {
        return PfqnSensMvaTest.logFactorial(n);
    }
}
