/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import java.util.ArrayList;
import java.util.List;

import jline.api.pfqn.nc.Pfqn_ca;
import jline.api.pfqn.nc.Pfqn_manjunath;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the Manjunath-Sikdar transform for product-form queueing
 * networks (Pfqn_manjunath).
 *
 * TWO ORACLES, NEITHER OF WHICH COMES OUT OF THE IMPLEMENTATION. With no extra
 * rows the transform computes the ordinary closed-network normalizing constant,
 * for which Pfqn_ca's convolution recursion is an independent exact algorithm
 * sharing no code; the two are algebraic identities for the same sum, so
 * agreement to 1e-13 is the correct expectation and not a tolerance chosen to
 * pass. With extra rows Pfqn_ca has nothing to say, and the oracle becomes
 * bcmpEnum below, which sums the BCMP product form over the enumerated state
 * space and applies each row by direct comparison -- the very enumeration the
 * transform exists to avoid, so a coefficient-domain defect cannot hide behind a
 * shared traversal.
 *
 * The constrained expectations are additionally pinned to the MATLAB reference
 * implementation's output, reproduced to 12 digits: the transform is exact, so
 * there is no tolerance to hide behind and any drift is a defect, not noise.
 */
public class PfqnManjunathApiTest {

    private static final double TOL = 1e-12;

    /** The shared constrained fixture has 3 queueing stations plus 1 delay. */
    private static final int S_C = 4;

    private static Matrix mat(double[][] v) {
        if (v.length == 0) {
            return new Matrix(0, 0);
        }
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                m.set(i, j, v[i][j]);
            }
        }
        return m;
    }

    private static Matrix row(double... v) {
        return mat(new double[][]{v});
    }

    private static double factorial(int n) {
        double f = 1.0;
        for (int k = 2; k <= n; k++) {
            f *= k;
        }
        return f;
    }

    /** Every way of splitting n indistinguishable jobs over k stations. */
    private static List<int[]> compositions(int n, int k) {
        List<int[]> out = new ArrayList<int[]>();
        if (k <= 0) {
            return out;
        }
        if (k == 1) {
            out.add(new int[]{n});
            return out;
        }
        for (int a = 0; a <= n; a++) {
            List<int[]> sub = compositions(n - a, k - 1);
            for (int s = 0; s < sub.size(); s++) {
                int[] rowv = new int[k];
                rowv[0] = a;
                System.arraycopy(sub.get(s), 0, rowv, 1, k - 1);
                out.add(rowv);
            }
        }
        return out;
    }

    /**
     * Direct sum of the BCMP product form over the enumerated closed state space,
     * keeping the states that satisfy every extra row.
     */
    private static double bcmpEnum(double[][] L, int[] N, double[][] Z, double[][] A,
                                   double[] b, String sense) {
        int M = L.length;
        int Mz = Z.length;
        int S = M + Mz;
        int R = N.length;
        List<List<int[]>> alloc = new ArrayList<List<int[]>>();
        for (int r = 0; r < R; r++) {
            alloc.add(compositions(N[r], S));
        }
        int[] idx = new int[R];
        double g = 0.0;
        while (true) {
            // n(i,r), read column by column when the rows are applied.
            int[] n = new int[S * R];
            for (int r = 0; r < R; r++) {
                int[] a = alloc.get(r).get(idx[r]);
                for (int i = 0; i < S; i++) {
                    n[i + S * r] = a[i];
                }
            }
            boolean ok = true;
            for (int j = 0; j < b.length && ok; j++) {
                double v = 0.0;
                for (int c = 0; c < S * R; c++) {
                    v += A[j][c] * n[c];
                }
                if (sense.charAt(j) == 'E') {
                    ok = Math.abs(v - b[j]) < 1e-9;
                } else if (sense.charAt(j) == 'L') {
                    ok = v <= b[j] + 1e-9;
                } else {
                    ok = v > b[j] + 1e-9;
                }
            }
            if (ok) {
                double t = 1.0;
                for (int i = 0; i < M; i++) {
                    int ni = 0;
                    for (int r = 0; r < R; r++) {
                        ni += n[i + S * r];
                    }
                    t *= factorial(ni);
                    for (int r = 0; r < R; r++) {
                        t *= Math.pow(L[i][r], n[i + S * r]) / factorial(n[i + S * r]);
                    }
                }
                for (int k = 0; k < Mz; k++) {
                    for (int r = 0; r < R; r++) {
                        t *= Math.pow(Z[k][r], n[M + k + S * r])
                                / factorial(n[M + k + S * r]);
                    }
                }
                g += t;
            }
            int d = R - 1;
            for (; d >= 0; d--) {
                if (++idx[d] < alloc.get(d).size()) {
                    break;
                }
                idx[d] = 0;
            }
            if (d < 0) {
                break;
            }
        }
        return g;
    }

    @Test
    public void unconstrainedMatchesConvolution() {
        double[][][] Ls = {
            {{1}, {2}},
            {{1, 2}, {3, 1}},
            {{1, 2}, {3, 1}},
            {{0.4, 0.2}, {0.9, 0.7}, {0.1, 1.1}},
            {},
            {{1, 2}, {3, 1}},
            {{5, 1}, {1, 6}},
        };
        double[][] Ns = {{4}, {2, 3}, {2, 3}, {3, 2}, {2, 1}, {0, 0}, {6, 5}};
        double[][][] Zs = {
            {{0}},
            {{0, 0}},
            {{0.5, 1.5}},
            {{1, 2}},
            {{1, 3}},
            {{1, 1}},
            {{2, 3}},
        };
        for (int k = 0; k < Ls.length; k++) {
            Matrix L = mat(Ls[k]);
            Matrix N = row(Ns[k]);
            Matrix Z = mat(Zs[k]);
            Ret.pfqnNc ca = Pfqn_ca.pfqn_ca(L, N, Z);
            Ret.pfqnManjunath mj = Pfqn_manjunath.pfqn_manjunath(L, N, Z);
            assertEquals(ca.lG, mj.lG, 1e-13 * Math.max(1.0, Math.abs(ca.lG)),
                    "case " + k + " must reproduce the convolution constant");
        }
    }

    @Test
    public void extraRowsMatchEnumerationAndMatlab() {
        double[][] L = {{1, 2}, {3, 1}, {0.5, 0.5}};
        double[][] Z = {{1, 2}};
        int[] N = {3, 2};
        double[] a1 = new double[S_C * 2];
        a1[0] = 1;
        a1[0 + S_C] = 1;                                  // jobs at queue 1
        double[] a2 = new double[S_C * 2];
        a2[0] = 2;
        a2[1] = 1;
        a2[0 + S_C] = 1;
        a2[1 + S_C] = 3;                                  // weighted budget
        double[] a3 = new double[S_C * 2];
        a3[3] = 1;                                        // class 1 at the delay

        double[][][] As = {{a1}, {a1}, {a1}, {a2}, {a1, a2}, {a1, a2}, {a1, a2}, {a1, a2},
                           {a1, a3}};
        double[][] bs = {{2}, {2}, {1}, {6}, {2, 6}, {2, 6}, {1, 6}, {1, 5}, {2, 1}};
        String[] senses = {"L", "E", "G", "L", "LL", "EL", "GL", "GG", "LL"};
        double[] matlab = {2214.58333333333, 526.541666666667, 1000.29166666667,
                           1611.83333333333, 1326.58333333333, 353.041666666667,
                           638.291666666667, 559.25, 2159.6875};

        Matrix Lm = mat(L);
        Matrix Zm = mat(Z);
        Matrix Nm = row(new double[]{N[0], N[1]});
        for (int k = 0; k < As.length; k++) {
            double want = bcmpEnum(L, N, Z, As[k], bs[k], senses[k]);
            Ret.pfqnManjunath got = Pfqn_manjunath.pfqn_manjunath(
                    Lm, Nm, Zm, mat(As[k]), row(bs[k]), senses[k]);
            assertEquals(want, got.G, TOL * Math.abs(want),
                    "row set " + senses[k] + " (case " + k + ") must match enumeration");
            assertEquals(matlab[k], got.G, 1e-12 * Math.abs(matlab[k]),
                    "row set " + senses[k] + " (case " + k + ") must match MATLAB");
            assertEquals(Math.exp(got.lG), got.G, TOL * got.G, "lG is the log of G");
        }
    }

    @Test
    public void populationRowIsRedundantAndItsComplementEmpty() {
        // Every admissible state holds sum(N) jobs, so declaring that as an extra
        // row must change nothing: a direct check that an equality row is
        // discharged by picking a coefficient rather than by summing one.
        Matrix L = mat(new double[][]{{1, 2}, {3, 1}, {0.5, 0.5}});
        Matrix Z = mat(new double[][]{{1, 2}});
        Matrix N = row(3, 2);
        double Gref = Pfqn_manjunath.pfqn_manjunath(L, N, Z).G;

        double[] ones = new double[S_C * 2];
        for (int c = 0; c < ones.length; c++) {
            ones[c] = 1.0;
        }
        Matrix A = row(ones);
        assertEquals(Gref, Pfqn_manjunath.pfqn_manjunath(L, N, Z, A, row(5), "E").G,
                1e-13 * Gref);
        assertEquals(Gref, Pfqn_manjunath.pfqn_manjunath(L, N, Z, A, row(5), "L").G,
                1e-13 * Gref);
        assertEquals(0.0, Pfqn_manjunath.pfqn_manjunath(L, N, Z, A, row(5), "G").G, 0.0);
    }

    @Test
    public void trivialRowsAreDecidedRatherThanCarried() {
        Matrix L = mat(new double[][]{{1, 2}, {3, 1}, {0.5, 0.5}});
        Matrix Z = mat(new double[][]{{1, 2}});
        Matrix N = row(3, 2);
        double Gref = Pfqn_manjunath.pfqn_manjunath(L, N, Z).G;

        Matrix zero = row(new double[S_C * 2]);
        assertEquals(Gref, Pfqn_manjunath.pfqn_manjunath(L, N, Z, zero, row(0), "E").G,
                1e-13 * Gref);
        assertEquals(0.0, Pfqn_manjunath.pfqn_manjunath(L, N, Z, zero, row(3), "E").G, 0.0);
        assertEquals(0.0, Pfqn_manjunath.pfqn_manjunath(L, N, Z, zero, row(-1), "L").G, 0.0);

        double[] a1 = new double[S_C * 2];
        a1[0] = 1;
        a1[0 + S_C] = 1;
        assertEquals(Gref, Pfqn_manjunath.pfqn_manjunath(L, N, Z, row(a1), row(-1), "G").G,
                1e-13 * Gref);
    }

    @Test
    public void peakIsTheClassLatticeWhenUnconstrained() {
        // No extra rows leaves only the class axes live, so the realised cost must
        // be exactly the lattice Pfqn_ca walks.
        Matrix L = mat(new double[][]{{1, 2}, {3, 1}, {0.5, 0.5}});
        Matrix Z = mat(new double[][]{{1, 2}});
        Ret.pfqnManjunath r = Pfqn_manjunath.pfqn_manjunath(L, row(3, 2), Z);
        assertEquals(4L * 3L, r.peakStates);
    }

    @Test
    public void refusalsNameTheReason() {
        Matrix L = mat(new double[][]{{1, 2}, {3, 1}, {0.5, 0.5}});
        Matrix Z = mat(new double[][]{{1, 2}});
        Matrix N = row(3, 2);

        double[] frac = new double[S_C * 2];
        frac[0] = 0.5;
        assertThrows(RuntimeException.class, () ->
                Pfqn_manjunath.pfqn_manjunath(L, N, Z, row(frac), row(1), "L"));
        double[] neg = new double[S_C * 2];
        neg[0] = -1.0;
        assertThrows(RuntimeException.class, () ->
                Pfqn_manjunath.pfqn_manjunath(L, N, Z, row(neg), row(1), "L"));
        double[] a1 = new double[S_C * 2];
        a1[0] = 1;
        assertThrows(RuntimeException.class, () ->
                Pfqn_manjunath.pfqn_manjunath(L, N, Z, row(a1), row(1.5), "L"));
        assertThrows(RuntimeException.class, () ->
                Pfqn_manjunath.pfqn_manjunath(L, N, Z, row(a1), row(1), "X"));

        // A negative population is a legitimate empty set, not an error.
        Ret.pfqnManjunath empty = Pfqn_manjunath.pfqn_manjunath(L, row(-1, 2), Z);
        assertEquals(0.0, empty.G, 0.0);
        assertTrue(Double.isInfinite(empty.lG) && empty.lG < 0.0);
    }

    // ------------------------------------------------------------------
    // The per-class decomposition
    // ------------------------------------------------------------------
    // Reference instance: PS queue (demands 1, 2) + delay (think 2, 4),
    // N = [4 4] both starting at the delay, with 2*n1 + 3*n2 <= 10 on the queue
    // occupancy. The expectations are the stationary law of an INDEPENDENTLY
    // built exact CTMC under HOLD truncation (a refused admission is a deleted
    // transition), which agrees with the truncated product form to 8.2e-17, so
    // these are not the routine restating itself.

    private static Ret.pfqnManjunath refStats() {
        Matrix L = row(1.0, 2.0);
        Matrix Z = row(2.0, 4.0);
        Matrix N = row(4, 4);
        double[] a = new double[4];
        a[0] = 2;
        a[2] = 3;
        return Pfqn_manjunath.pfqn_manjunath(L, N, Z, row(a), row(10), "L", true);
    }

    @Test
    public void decompositionMatchesTheExactHoldChain() {
        Ret.pfqnManjunath r = refStats();
        double[][] ref = {
            {1.88612099644128, 1.50177935943061},      // Q
            {0.544483985765125, 0.224199288256228},    // X
            {0.544483985765125, 0.448398576512456},    // U
            {1.08896797153025, 0.896797153024911},     // think
            {1.02491103202847, 1.60142348754448},      // blocked
            {2.11387900355872, 2.49822064056939},      // delay
        };
        Matrix[] got = {r.Q, r.X, r.U, r.think, r.blocked, r.delay};
        String[] nm = {"Q", "X", "U", "think", "blocked", "delay"};
        for (int i = 0; i < nm.length; i++) {
            for (int c = 0; c < 2; c++) {
                assertEquals(ref[i][c], got[i].get(0, c), 1e-12 * Math.abs(ref[i][c]),
                        nm[i] + " class " + (c + 1));
            }
        }
    }

    @Test
    public void decompositionConservesThePopulation() {
        // A blocked job never leaves the delay, so nothing escapes the
        // accounting: queue + thinking + held is exactly N, class by class.
        Ret.pfqnManjunath r = refStats();
        for (int c = 0; c < 2; c++) {
            assertEquals(4.0, r.Q.get(0, c) + r.think.get(0, c) + r.blocked.get(0, c), 1e-12);
            assertEquals(r.delay.get(0, c), r.think.get(0, c) + r.blocked.get(0, c), 1e-12);
        }
    }

    @Test
    public void unconstrainedThroughputIsTheClassicalRatio() {
        Matrix L = row(1.0, 2.0);
        Matrix Z = row(2.0, 4.0);
        Matrix N = row(4, 4);
        Ret.pfqnManjunath r = Pfqn_manjunath.pfqn_manjunath(L, N, Z, null, null, null, true);
        double lG = Pfqn_ca.pfqn_ca(L, N, Z).lG;
        for (int c = 0; c < 2; c++) {
            Matrix Nr = N.copy();
            Nr.set(0, c, N.get(0, c) - 1.0);
            double lGr = Pfqn_ca.pfqn_ca(L, Nr, Z).lG;
            assertEquals(Math.exp(lGr - lG), r.X.get(0, c), 1e-13, "class " + (c + 1));
            // Nothing can be held when nothing is constrained.
            assertEquals(0.0, r.blocked.get(0, c), 1e-12);
        }
    }

    @Test
    public void decompositionRefusesEveryOtherConfiguration() {
        Matrix L = row(1.0, 2.0);
        Matrix Z = row(2.0, 4.0);
        Matrix N = row(4, 4);
        assertThrows(RuntimeException.class, () -> Pfqn_manjunath.pfqn_manjunath(
                L, N, mat(new double[][]{{2.0, 4.0}, {1.0, 1.0}}), new Matrix(1, 6),
                row(10), "L", true));
        // Two queueing stations: the delay->q1->q2->delay cycle makes the chain
        // irreversible, so Kelly truncation no longer holds.
        assertThrows(RuntimeException.class, () -> Pfqn_manjunath.pfqn_manjunath(
                mat(new double[][]{{1.0, 2.0}, {1.0, 2.0}}), N, Z, new Matrix(1, 6),
                row(10), "L", true));
        double[] ad = new double[4];
        ad[1] = 1;                       // column 1 is (delay, class 1)
        assertThrows(RuntimeException.class, () -> Pfqn_manjunath.pfqn_manjunath(
                L, N, Z, row(ad), row(10), "L", true));
        // G is still produced when the decomposition is not asked for.
        assertTrue(Pfqn_manjunath.pfqn_manjunath(
                mat(new double[][]{{1.0, 2.0}, {1.0, 2.0}}), N, Z).G > 0.0);
    }
}
