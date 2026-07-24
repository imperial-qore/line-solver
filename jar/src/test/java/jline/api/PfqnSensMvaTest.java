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
import jline.api.pfqn.sens.Pfqn_sens_mva;
import jline.api.pfqn.sens.Pfqn_sens;
import jline.io.Ret;
import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates {@link Pfqn_sens_mva#pfqn_sens_mva}, the MVA-like moment recursion of
 * de Souza e Silva and Muntz (1988), Corollary 1, against three independent
 * references, mirroring the MATLAB harness {@code pfqn_sens_mva_validate.m}:
 *
 * <ul>
 *   <li>A. brute-force enumeration of the closed product-form equilibrium
 *       distribution (ground truth, mi==1);</li>
 *   <li>B. the differentiated-MVA Jacobian of {@link Pfqn_sens}, via the identity
 *       Cov[n(i,r),n(i,s)] = L(i,s) * dQ(i,r)/dL(i,s) (covers mi&gt;1);</li>
 *   <li>C. {@link Pfqn_mva} for the base measures X, Q, U, R.</li>
 * </ul>
 *
 * <p>It also checks the raw asymmetry of QCov before symmetrization: the
 * recursion computes W(k,j;t,j) and W(t,j;k,j) by numerically distinct
 * expressions, so their agreement is a nontrivial check of the formula.</p>
 */
public class PfqnSensMvaTest {

    private static final double TOL_BRUTE = 1e-9;
    private static final double TOL_SENS = 1e-9;
    private static final double TOL_MVA = 1e-10;
    private static final double TOL_SYM = 1e-9;

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

    private static double[] flat(Matrix[] arr) {
        int n = arr.length * arr[0].getNumRows() * arr[0].getNumCols();
        double[] v = new double[n];
        int k = 0;
        for (int i = 0; i < arr.length; i++) {
            for (int r = 0; r < arr[i].getNumRows(); r++) {
                for (int s = 0; s < arr[i].getNumCols(); s++) {
                    v[k++] = arr[i].get(r, s);
                }
            }
        }
        return v;
    }

    @Test
    public void sensMvaMatchesBruteForceSensAndMva() {
        Random rng = new Random(0);
        double errBrute = 0.0;
        double errSens = 0.0;
        double errMva = 0.0;
        double errSym = 0.0;
        int nBrute = 0;
        int nSens = 0;

        for (int trial = 1; trial <= 40; trial++) {
            int M = 1 + rng.nextInt(3);
            int R = 1 + rng.nextInt(3);
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
            // exercise a zero-demand column now and then: class r never visits station i
            if (trial % 5 == 0 && M > 1) {
                L.set(0, 0, 0.0);
            }

            Ret.pfqnSensMva mom = Pfqn_sens_mva.pfqn_sens_mva(L, N, Z);

            // ---- C. base measures against pfqn_mva --------------------------
            Ret.pfqnMVA mva = Pfqn_mva.pfqn_mva(L, N, Z);
            errMva = Math.max(errMva, relerr(mom.X, mva.X));
            errMva = Math.max(errMva, relerr(mom.Q, mva.Q));
            errMva = Math.max(errMva, relerr(mom.U, mva.U));
            // The residence time returned by Pfqn_mva is CN straight off the
            // recursion, CN(i,s) = L(i,s)*(mi(i) + Qtot(i|N-e_s)), matching
            // MATLAB pfqn_mva and pfqn_sens_mva. It was previously re-derived as
            // QN/XN, which is equivalent wherever XN(r) > 0 (QN = XN*CN is the
            // very equation that produced QN) but evaluated 0/0 on an empty
            // class and returned NaN where MATLAB is finite. Both must now agree
            // everywhere, empty classes included, and neither may be NaN.
            for (int r = 0; r < R; r++) {
                for (int i = 0; i < M; i++) {
                    assertTrue(!Double.isNaN(mva.R.get(i, r)),
                            "Pfqn_mva residence time must be finite even on an empty class; class "
                                    + r + " population " + N.get(0, r) + " gave " + mva.R.get(i, r));
                    assertTrue(!Double.isNaN(mom.R.get(i, r)),
                            "pfqn_sens_mva residence time must stay finite on empty classes");
                    errMva = Math.max(errMva, relerr(new double[]{mom.R.get(i, r)},
                            new double[]{mva.R.get(i, r)}));
                }
            }

            // ---- A. brute force ---------------------------------------------
            double totpop = 1.0;
            for (int r = 0; r < R; r++) {
                totpop *= (N.get(0, r) + 1);
            }
            if (totpop <= 64 && M <= 3) {
                BruteMoments bm = bruteMoments(L, N, Z);
                errBrute = Math.max(errBrute, relerr(mom.Q, bm.Q));
                errBrute = Math.max(errBrute, relerr(flat(mom.QCov), flat(bm.QCov)));
                nBrute++;
            }

            // ---- B. pfqn_sens Jacobian, and the raw asymmetry ---------------
            for (int micase = 1; micase <= 2; micase++) {
                Matrix mi = new Matrix(1, M);
                for (int i = 0; i < M; i++) {
                    mi.set(0, i, micase == 1 ? 1.0 : 1 + rng.nextInt(3));
                }
                Ret.pfqnSensMva momi = Pfqn_sens_mva.pfqn_sens_mva(L, N, Z, mi);
                Ret.pfqnSens sens = Pfqn_sens.pfqn_sens(L, N, Z, mi);
                // Read the reference off the raw Jacobian, NOT off sens.QCov:
                // Pfqn_sens now sources its same-station blocks from
                // Pfqn_sens_mva, so comparing against sens.QCov would compare the
                // recursion with itself.
                int[][] pL = new int[M][R];
                boolean[][] pLset = new boolean[M][R];
                for (int p = 0; p < sens.paramType.length; p++) {
                    if (sens.paramType[p] == 0) {
                        pL[sens.paramStation[p]][sens.paramClass[p]] = p;
                        pLset[sens.paramStation[p]][sens.paramClass[p]] = true;
                    }
                }
                Matrix[] covRef = new Matrix[M];
                for (int i = 0; i < M; i++) {
                    covRef[i] = new Matrix(R, R);
                    for (int r = 0; r < R; r++) {
                        for (int s = 0; s < R; s++) {
                            if (pLset[i][s]) {
                                covRef[i].set(r, s, L.get(i, s) * sens.dQ[pL[i][s]].get(i, r));
                            }
                        }
                    }
                }
                errSens = Math.max(errSens, relerr(flat(momi.QCov), flat(covRef)));
                errSym = Math.max(errSym, momi.QCovAsym);
                // Theorem 3: variance of the total queue length at a station
                Matrix totRef = new Matrix(M, 1);
                for (int i = 0; i < M; i++) {
                    double t = 0.0;
                    for (int r = 0; r < R; r++) {
                        for (int s = 0; s < R; s++) {
                            t += covRef[i].get(r, s);
                        }
                    }
                    totRef.set(i, 0, t);
                }
                errSens = Math.max(errSens, relerr(momi.QTotVar, totRef));
                nSens++;
            }
        }

        assertTrue(nBrute > 0, "no model was small enough for brute-force enumeration");
        assertTrue(nSens > 0, "no model was checked against the Jacobian");
        assertTrue(errBrute <= TOL_BRUTE,
                "brute-force product form (" + nBrute + " models): " + errBrute);
        assertTrue(errSens <= TOL_SENS,
                "pfqn_sens Jacobian (" + nSens + " models): " + errSens);
        assertTrue(errMva <= TOL_MVA, "pfqn_mva base measures: " + errMva);
        assertTrue(errSym <= TOL_SYM, "QCov symmetry (raw, pre-symmetrize): " + errSym);
    }

    /**
     * pfqn_sens must expose the same-station covariance blocks of the recursion
     * and the Jacobian-sourced cross-station ones, and its QVar/QTotVar/QCovAsym
     * must come from the recursion.
     */
    @Test
    public void sensSourcesSameStationBlocksFromSensMva() {
        Matrix L = new Matrix(2, 2);
        L.set(0, 0, 0.9);
        L.set(0, 1, 0.4);
        L.set(1, 0, 0.3);
        L.set(1, 1, 0.7);
        Matrix N = row(3, 2);
        Matrix Z = row(0.5, 0.8);
        Matrix mi = row(1, 2);

        Ret.pfqnSens sens = Pfqn_sens.pfqn_sens(L, N, Z, mi);
        Ret.pfqnSensMva mom = Pfqn_sens_mva.pfqn_sens_mva(L, N, Z, mi);

        double errSame = 0.0;
        double errCross = 0.0;
        int[][] pL = new int[2][2];
        for (int p = 0; p < sens.paramType.length; p++) {
            if (sens.paramType[p] == 0) {
                pL[sens.paramStation[p]][sens.paramClass[p]] = p;
            }
        }
        for (int i = 0; i < 2; i++) {
            for (int r = 0; r < 2; r++) {
                for (int j = 0; j < 2; j++) {
                    for (int s = 0; s < 2; s++) {
                        double got = sens.QCov[i][r].get(j, s);
                        if (i == j) {
                            errSame = Math.max(errSame, Math.abs(got - mom.QCov[i].get(r, s)));
                        } else {
                            double ref = L.get(j, s) * sens.dQ[pL[j][s]].get(i, r);
                            errCross = Math.max(errCross, Math.abs(got - ref));
                        }
                    }
                }
            }
        }
        assertTrue(errSame == 0.0, "same-station blocks not sourced from pfqn_sens_mva: " + errSame);
        assertTrue(errCross == 0.0, "cross-station blocks not read off the Jacobian: " + errCross);
        assertTrue(relerr(sens.QVar, mom.QVar) == 0.0, "QVar not sourced from pfqn_sens_mva");
        assertTrue(relerr(sens.QTotVar, mom.QTotVar) == 0.0, "QTotVar not sourced from pfqn_sens_mva");
        assertTrue(sens.QCovAsym == mom.QCovAsym, "QCovAsym not sourced from pfqn_sens_mva");
        assertTrue(sens.QCovAsym <= TOL_SYM, "QCov symmetry residual: " + sens.QCovAsym);

        // The full covariance matrix over the M*R queue-length variables must be
        // symmetric: Cov[n(i,r),n(j,s)] = Cov[n(j,s),n(i,r)].
        double errFullSym = 0.0;
        for (int i = 0; i < 2; i++) {
            for (int r = 0; r < 2; r++) {
                for (int j = 0; j < 2; j++) {
                    for (int s = 0; s < 2; s++) {
                        errFullSym = Math.max(errFullSym,
                                Math.abs(sens.QCov[i][r].get(j, s) - sens.QCov[j][s].get(i, r)));
                    }
                }
            }
        }
        assertTrue(errFullSym <= TOL_SYM, "full QCov symmetry: " + errFullSym);
    }

    // =====================================================================
    private static final class BruteMoments {
        Matrix Q;
        Matrix[] QCov;
    }

    /**
     * Exact moments by enumerating the closed product-form equilibrium
     * distribution. Stations 1..M are single-server fixed-rate centers; the think
     * time Z is an infinite-server station indexed 0 and carries no moment.
     * p(n) ~ prod_i [ n_i! prod_r L(i,r)^n(i,r)/n(i,r)! ] * prod_r Z(r)^n(0,r)/n(0,r)!
     */
    private static BruteMoments bruteMoments(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        List<int[][]> states = enumerateStates(N, M, R);
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

        BruteMoments out = new BruteMoments();
        out.Q = new Matrix(M, R);
        for (int k = 0; k < K; k++) {
            int[][] nir = states.get(k);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    out.Q.set(i, r, out.Q.get(i, r) + w[k] * nir[i][r]);
                }
            }
        }
        out.QCov = new Matrix[M];
        for (int i = 0; i < M; i++) {
            out.QCov[i] = new Matrix(R, R);
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    double m2 = 0.0;
                    for (int k = 0; k < K; k++) {
                        int[][] nir = states.get(k);
                        m2 += w[k] * nir[i][r] * nir[i][s];
                    }
                    out.QCov[i].set(r, s, m2 - out.Q.get(i, r) * out.Q.get(i, s));
                }
            }
        }
        return out;
    }

    /** All allocations of N(r) class-r jobs over M stations (remainder in the delay). */
    private static List<int[][]> enumerateStates(Matrix N, int M, int R) {
        List<List<int[]>> per = new ArrayList<List<int[]>>();
        for (int r = 0; r < R; r++) {
            per.add(compositionsLeq((int) N.get(0, r), M));
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

    /** All nonnegative integer vectors of length M summing to at most n. */
    static List<int[]> compositionsLeq(int n, int M) {
        List<int[]> out = new ArrayList<int[]>();
        if (M == 1) {
            for (int i = 0; i <= n; i++) {
                out.add(new int[]{i});
            }
            return out;
        }
        for (int first = 0; first <= n; first++) {
            List<int[]> sub = compositionsLeq(n - first, M - 1);
            for (int j = 0; j < sub.size(); j++) {
                int[] s = sub.get(j);
                int[] v = new int[M];
                v[0] = first;
                System.arraycopy(s, 0, v, 1, M - 1);
                out.add(v);
            }
        }
        return out;
    }

    static double logFactorial(int n) {
        double s = 0.0;
        for (int i = 2; i <= n; i++) {
            s += Math.log(i);
        }
        return s;
    }
}
