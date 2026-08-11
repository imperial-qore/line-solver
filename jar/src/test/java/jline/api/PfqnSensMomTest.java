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
import jline.api.pfqn.sens.Pfqn_sens_mom;
import jline.api.pfqn.sens.Pfqn_sens_mva;
import jline.io.Ret;
import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates {@link Pfqn_sens_mom#pfqn_sens_mom}, the higher-moment analysis of
 * Strelen (1990), against five independent references, mirroring the MATLAB
 * harness {@code pfqn_sens_mom_validate.m}:
 *
 * <ul>
 *   <li>A. brute-force enumeration of the closed product-form distribution, which
 *       is ground truth for m, Var, Cov, E[Q^2] and E[Q^3];</li>
 *   <li>B. {@link Pfqn_sens_mva}. Summing its per-class covariance matrix at
 *       station i over all class pairs must give Var[Q_i], since
 *       Var[sum_r n(i,r)] = sum_{r,s} Cov[n(i,r),n(i,s)]. This ties the
 *       per-station-total moments of Strelen to the finer per-class moments of de
 *       Souza e Silva and Muntz;</li>
 *   <li>C. {@link Pfqn_mva} for the base measures;</li>
 *   <li>D. Cov symmetry: x_j dm_i/dx_j and x_i dm_j/dx_i are computed by different
 *       derivative tracks and must agree;</li>
 *   <li>E. the published table of Example 3.4 of the reference (the Kobayashi
 *       central-server model), which pins the second derivative against numbers
 *       the author printed rather than against our own code.</li>
 * </ul>
 *
 * <p>Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks and
 * its Linearizer", Performance Evaluation 11:127-142, 1990.</p>
 */
public class PfqnSensMomTest {

    private static final double TOL_BRUTE = 1e-9;
    private static final double TOL_MVA = 1e-10;
    private static final double TOL_TOT = 1e-9;
    private static final double TOL_SYM = 1e-9;
    /**
     * The per-class grouping must reproduce brute force and Pfqn_sens_mva's second
     * moments; the MATLAB validator measures 1.55e-15 on the latter.
     */
    private static final double TOL_GRP = 1e-9;
    /** the paper prints 5 significant digits */
    private static final double TOL_PAPER = 5e-5;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
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

    @Test
    public void sensMomMatchesBruteForceSensMvaAndMva() {
        Random rng = new Random(3);
        double errBrute = 0.0;
        double errGrp = 0.0;
        double errMva = 0.0;
        double errTot = 0.0;
        double errSym = 0.0;
        int nBrute = 0;
        int nGrp = 0;

        for (int trial = 1; trial <= 40; trial++) {
            int M = 1 + rng.nextInt(3);
            int R = 1 + rng.nextInt(2);
            Matrix L = new Matrix(M, R);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    L.set(i, r, 0.2 + rng.nextDouble());
                }
            }
            Matrix N = new Matrix(1, R);
            boolean anyPop = false;
            for (int r = 0; r < R; r++) {
                double nr = rng.nextInt(4);
                N.set(0, r, nr);
                anyPop = anyPop || nr > 0;
            }
            if (!anyPop) {
                N.set(0, 0, 2);
            }
            Matrix Z = new Matrix(1, R);
            if (trial % 2 == 0) {
                for (int r = 0; r < R; r++) {
                    Z.set(0, r, 0.3 + rng.nextDouble());
                }
            }
            Matrix mi = Matrix.ones(1, M);
            boolean plainMi = true;
            if (trial % 3 == 0) {
                for (int i = 0; i < M; i++) {
                    mi.set(0, i, 1 + rng.nextInt(3));
                }
                plainMi = false;
            }

            Ret.pfqnSensMom mom = Pfqn_sens_mom.pfqn_sens_mom(L, N, Z, mi);

            // ---- C. base measures ------------------------------------------
            Ret.pfqnMVA mva = Pfqn_mva.pfqn_mva(L, N, Z, mi);
            errMva = Math.max(errMva, relerr(mom.X, mva.X));
            errMva = Math.max(errMva, relerr(mom.Q, mva.Q));
            errMva = Math.max(errMva, relerr(mom.U, mva.U));
            // The residence time returned by Pfqn_mva is CN straight off the
            // recursion, CN(i,s) = L(i,s)*(mi(i) + Qtot(i|N-e_s)), matching
            // MATLAB pfqn_mva and pfqn_sens_mom. It was previously re-derived as
            // QN/XN, which agrees wherever XN(r) > 0 but evaluated 0/0 on an
            // empty class and returned NaN where MATLAB is finite. Both must now
            // agree everywhere, empty classes included.
            for (int r = 0; r < R; r++) {
                for (int i = 0; i < M; i++) {
                    assertTrue(!Double.isNaN(mva.R.get(i, r)),
                            "Pfqn_mva residence time must be finite even on an empty class; class "
                                    + r + " population " + N.get(0, r) + " gave " + mva.R.get(i, r));
                    assertTrue(!Double.isNaN(mom.R.get(i, r)),
                            "pfqn_sens_mom residence time must stay finite on empty classes");
                    errMva = Math.max(errMva, relerr(new double[]{mom.R.get(i, r)},
                            new double[]{mva.R.get(i, r)}));
                }
            }

            // ---- D. symmetry -----------------------------------------------
            errSym = Math.max(errSym, mom.CovAsym);

            // ---- B. total variance against the per-class covariances -------
            Ret.pfqnSensMva ref = Pfqn_sens_mva.pfqn_sens_mva(L, N, Z, mi);
            errTot = Math.max(errTot, relerr(mom.Var, ref.QTotVar));

            // ---- A. brute force --------------------------------------------
            double totpop = 1.0;
            for (int r = 0; r < R; r++) {
                totpop *= (N.get(0, r) + 1);
            }
            if (totpop <= 32 && M <= 3 && plainMi) {
                BruteTotals bt = bruteTotals(L, N, Z);
                errBrute = Math.max(errBrute, relerr(mom.m, bt.m));
                errBrute = Math.max(errBrute, relerr(mom.Var, bt.Var));
                errBrute = Math.max(errBrute, relerr(mom.Cov, bt.Cov));
                errBrute = Math.max(errBrute, relerr(mom.M2, bt.M2));
                errBrute = Math.max(errBrute, relerr(mom.M3, bt.M3));
                nBrute++;

                // ---- F. the per-class grouping ------------------------------
                // groups = 1..R scales one class at a time, which is Akyildiz and
                // Strelen's Theorem 1 with T = {r}. It must reproduce the per-class
                // moments of the brute-force distribution, INCLUDING the third, and
                // its second moments must equal Pfqn_sens_mva's exactly.
                Ret.pfqnSensMom momc = Pfqn_sens_mom.pfqn_sens_mom(L, N, Z, mi,
                        Pfqn_sens_mom.perClassGroups(R));
                int[] perClass = new int[R];
                for (int r = 0; r < R; r++) {
                    perClass[r] = r;
                }
                BruteTotals bc = bruteGroups(L, N, Z, perClass, R);
                errGrp = Math.max(errGrp, relerr(momc.m, bc.m));
                errGrp = Math.max(errGrp, relerr(momc.Var, bc.Var));
                errGrp = Math.max(errGrp, relerr(momc.M2, bc.M2));
                errGrp = Math.max(errGrp, relerr(momc.M3, bc.M3));
                errGrp = Math.max(errGrp, relerr(momc.Var, ref.QVar));
                assertEquals(R, momc.G);
                assertTrue(momc.CovG != null, "CovG must always be populated");
                // the collapsed (M x M) Cov view exists only where the group index
                // carries no information, which for the per-class grouping means the
                // single-class case; a random trial may well draw R == 1
                assertTrue((momc.Cov != null) == (R == 1),
                        "the collapsed Cov view is populated exactly when G == 1");
                nGrp++;

                // ---- G. an intermediate grouping ----------------------------
                // A grouping that puts every class in one group must reproduce the
                // default (station totals).
                if (R == 2) {
                    Ret.pfqnSensMom momg = Pfqn_sens_mom.pfqn_sens_mom(L, N, Z, mi, row(1, 1));
                    errGrp = Math.max(errGrp, relerr(momg.m, mom.m));
                    errGrp = Math.max(errGrp, relerr(momg.Var, mom.Var));
                    errGrp = Math.max(errGrp, relerr(momg.M3, mom.M3));
                }
            }
        }

        assertTrue(nBrute > 0, "no model was small enough for brute-force enumeration");
        assertTrue(nGrp > 0, "no model exercised the groups argument");
        assertTrue(errGrp <= TOL_GRP, "grouped moments (" + nGrp + " models): " + errGrp);
        assertTrue(errBrute <= TOL_BRUTE,
                "brute-force product form (" + nBrute + " models): " + errBrute);
        assertTrue(errTot <= TOL_TOT, "Var vs pfqn_sens_mva QTotVar: " + errTot);
        assertTrue(errMva <= TOL_MVA, "pfqn_mva base measures: " + errMva);
        assertTrue(errSym <= TOL_SYM, "Cov symmetry (raw, pre-symmetrize): " + errSym);
    }

    /**
     * E. Example 3.4 of the reference: Kobayashi central-server model. 12 type-1
     * queues, one class. Queues 1-9: x=0.0215, e=9.333; queues 10,11: x=0.104,
     * e=10.5; queue 12: x=0.019, e=105. The paper prints E(Q_i) and
     * sigma^2_{Q_i} for n=3, 2, 1.
     */
    @Test
    public void sensMomMatchesStrelenExample34() {
        double[] xs = new double[12];
        double[] es = new double[12];
        for (int i = 0; i < 9; i++) {
            xs[i] = 0.0215;
            es[i] = 9.333;
        }
        xs[9] = 0.104;
        es[9] = 10.5;
        xs[10] = 0.104;
        es[10] = 10.5;
        xs[11] = 0.019;
        es[11] = 105;
        Matrix Lk = new Matrix(12, 1);
        for (int i = 0; i < 12; i++) {
            Lk.set(i, 0, xs[i] * es[i]);
        }
        // rows: queues 1-9, 10-11, 12; columns: n = 3, 2, 1
        double[][] paperM = {
            {0.07606, 0.05835, 0.03353},
            {0.53316, 0.36327, 0.18246},
            {1.24917, 0.74835, 0.33334}};
        double[][] paperVar = {
            {0.07893, 0.05873, 0.03240},
            {0.57689, 0.34341, 0.14917},
            {1.02546, 0.56250, 0.22222}};

        double errPaper = 0.0;
        for (int c = 0; c < 3; c++) {
            int nJobs = 3 - c;   // c = 0 -> n = 3, c = 1 -> n = 2, c = 2 -> n = 1
            Ret.pfqnSensMom mk = Pfqn_sens_mom.pfqn_sens_mom(Lk, row(nJobs), row(0.0));
            double[] got = {mk.m.get(0, 0), mk.m.get(9, 0), mk.m.get(11, 0)};
            double[] gotV = {mk.Var.get(0, 0), mk.Var.get(9, 0), mk.Var.get(11, 0)};
            double[] refM = {paperM[0][c], paperM[1][c], paperM[2][c]};
            double[] refV = {paperVar[0][c], paperVar[1][c], paperVar[2][c]};
            errPaper = Math.max(errPaper, relerr(got, refM));
            errPaper = Math.max(errPaper, relerr(gotV, refV));
            // the nine identical queues must be identical, and so must 10 and 11
            double[] nine = new double[9];
            double[] first = new double[9];
            for (int i = 0; i < 9; i++) {
                nine[i] = mk.m.get(i, 0);
                first[i] = mk.m.get(0, 0);
            }
            errPaper = Math.max(errPaper, relerr(nine, first));
            errPaper = Math.max(errPaper,
                    relerr(new double[]{mk.m.get(9, 0)}, new double[]{mk.m.get(10, 0)}));
        }
        assertTrue(errPaper <= TOL_PAPER,
                "Strelen Example 3.4 published table: " + errPaper);
    }

    // =====================================================================
    private static final class BruteTotals {
        Matrix m;
        Matrix Var;
        Matrix Cov;
        Matrix M2;
        Matrix M3;
    }

    /**
     * Moments of the per-station total queue lengths by enumerating the closed
     * product-form equilibrium distribution.
     */
    private static BruteTotals bruteTotals(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        List<int[][]> states = enumerateStates(N, M, R);
        int K = states.size();
        double[] w = stateWeights(L, N, Z, states);
        double[][] q = new double[K][M];
        for (int k = 0; k < K; k++) {
            int[][] nir = states.get(k);
            for (int i = 0; i < M; i++) {
                int s = 0;
                for (int r = 0; r < R; r++) {
                    s += nir[i][r];
                }
                q[k][i] = s;
            }
        }

        BruteTotals out = new BruteTotals();
        out.m = new Matrix(M, 1);
        out.M2 = new Matrix(M, 1);
        out.M3 = new Matrix(M, 1);
        out.Var = new Matrix(M, 1);
        out.Cov = new Matrix(M, M);
        for (int i = 0; i < M; i++) {
            double m1 = 0.0;
            double m2 = 0.0;
            double m3 = 0.0;
            for (int k = 0; k < K; k++) {
                m1 += w[k] * q[k][i];
                m2 += w[k] * q[k][i] * q[k][i];
                m3 += w[k] * q[k][i] * q[k][i] * q[k][i];
            }
            out.m.set(i, 0, m1);
            out.M2.set(i, 0, m2);
            out.M3.set(i, 0, m3);
            out.Var.set(i, 0, m2 - m1 * m1);
        }
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                double c = 0.0;
                for (int k = 0; k < K; k++) {
                    c += w[k] * q[k][i] * q[k][j];
                }
                out.Cov.set(i, j, c - out.m.get(i, 0) * out.m.get(j, 0));
            }
        }
        return out;
    }

    /**
     * The normalized closed product-form equilibrium probability of every state.
     */
    private static double[] stateWeights(Matrix L, Matrix N, Matrix Z, List<int[][]> states) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        int K = states.size();
        double[] w = new double[K];
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
                for (int r = 0; r < R; r++) {
                    if (nir[i][r] > 0) {
                        if (L.get(i, r) <= 0) {
                            ok = false;
                            break;
                        }
                        lw += nir[i][r] * Math.log(L.get(i, r)) - logFactorial(nir[i][r]);
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
        }
        double tot = 0.0;
        for (int k = 0; k < K; k++) {
            tot += w[k];
        }
        for (int k = 0; k < K; k++) {
            w[k] /= tot;
        }
        return w;
    }

    /**
     * Moments of Q_(i,g) = sum_(r in group g) n(i,r) by enumerating the closed
     * product-form equilibrium distribution. The grouped generalization of
     * {@link #bruteTotals}: {@code grp} all-zero with G == 1 reproduces it, and
     * {@code grp[r] == r} with G == R gives the per-class moments.
     */
    private static BruteTotals bruteGroups(Matrix L, Matrix N, Matrix Z, int[] grp, int G) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        List<int[][]> states = enumerateStates(N, M, R);
        int K = states.size();
        double[] w = stateWeights(L, N, Z, states);

        BruteTotals out = new BruteTotals();
        out.m = new Matrix(M, G);
        out.M2 = new Matrix(M, G);
        out.M3 = new Matrix(M, G);
        out.Var = new Matrix(M, G);
        for (int i = 0; i < M; i++) {
            for (int g = 0; g < G; g++) {
                double m1 = 0.0;
                double m2 = 0.0;
                double m3 = 0.0;
                for (int k = 0; k < K; k++) {
                    int[][] nir = states.get(k);
                    double q = 0.0;
                    for (int r = 0; r < R; r++) {
                        if (grp[r] == g) {
                            q += nir[i][r];
                        }
                    }
                    m1 += w[k] * q;
                    m2 += w[k] * q * q;
                    m3 += w[k] * q * q * q;
                }
                out.m.set(i, g, m1);
                out.M2.set(i, g, m2);
                out.M3.set(i, g, m3);
                out.Var.set(i, g, m2 - m1 * m1);
            }
        }
        return out;
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

    private static double logFactorial(int n) {
        return PfqnSensMvaTest.logFactorial(n);
    }
}
