/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.ld;

import java.util.Random;

import jline.api.pfqn.mva.Pfqn_mvac;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the Section V extension of MVAC to queue-length dependent (QLD) service
 * centers (Conway, de Souza e Silva and Lavenberg, IEEE TC 38(3):432-442, 1989).
 *
 * <p>Four independent oracles are used, since MVAC-LD is exact and must not merely be
 * self-consistent:</p>
 * <ol>
 *   <li>mu = 1 must reduce EXACTLY to {@link Pfqn_mvac} (Sections II-IV), a structurally
 *       different recursion that propagates means rather than distributions.</li>
 *   <li>Load-dependent results must match {@link Pfqn_mvald}, the classic load-dependent
 *       MVA population recursion.</li>
 *   <li>Eq. (25) is self-normalizing, so sum_n P_j(n) = 1 identically; a deviation means
 *       the implementation is broken, not the model.</li>
 *   <li>The marginals must reproduce the means: Q_j = sum_n n P_j(n).</li>
 * </ol>
 * <p>The three chain-resolution paths (part 1 / part 2 / label interchange) are covered
 * explicitly, as in {@code Pfqn_mvacTest}.</p>
 */
public class Pfqn_mvacldTest {

    private static final double EXACT_TOL = 1e-9;

    private static Matrix mat(double[][] v) {
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                m.set(i, j, v[i][j]);
            }
        }
        return m;
    }

    /** Rates of a c-server queue at each of M centers, mu(n) = min(n,c). */
    private static Matrix msRates(int M, int Nt, int c) {
        Matrix mu = new Matrix(M, Nt);
        for (int i = 0; i < M; i++) {
            for (int n = 1; n <= Nt; n++) {
                mu.set(i, n - 1, Math.min(n, c));
            }
        }
        return mu;
    }

    private static Matrix ones(int M, int Nt) {
        return Matrix.ones(M, Nt);
    }

    // ---- oracle 1: mu = 1 must reduce exactly to Pfqn_mvac ----

    @Test
    public void testFixedRateReducesToMvacNoDelay() {
        assertReducesToMvac(mat(new double[][]{{0.4, 0.7, 0.2}, {0.9, 0.3, 0.5}}),
                mat(new double[][]{{1, 2, 1}}), new Matrix(1, 3)); // S=D
    }

    @Test
    public void testFixedRateReducesToMvacAllDelay() {
        assertReducesToMvac(mat(new double[][]{{0.4, 0.7, 0.2}, {0.9, 0.3, 0.5}}),
                mat(new double[][]{{1, 2, 1}}), mat(new double[][]{{1.0, 0.5, 2.0}})); // S=0
    }

    @Test
    public void testFixedRateReducesToMvacMixed() {
        assertReducesToMvac(mat(new double[][]{{0.4, 0.7, 0.2}, {0.9, 0.3, 0.5}}),
                mat(new double[][]{{2, 1, 2}}), mat(new double[][]{{1.0, 0.0, 0.5}})); // 0<S<D
    }

    @Test
    public void testFixedRateReducesToMvacIdenticalClasses() {
        assertReducesToMvac(mat(new double[][]{{0.4, 0.4, 0.2}, {0.9, 0.9, 0.5}}),
                mat(new double[][]{{2, 3, 1}}), mat(new double[][]{{1.0, 1.0, 0.0}}));
    }

    @Test
    public void testFixedRateReducesToMvacTwoDelayCenters() {
        assertReducesToMvac(mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}}),
                mat(new double[][]{{2, 2}}), mat(new double[][]{{1.0, 0.0}, {0.5, 2.0}}));
    }

    /** A null rate matrix must default to single-server fixed rate. */
    @Test
    public void testDefaultRateMatrixIsFixedRate() {
        Matrix L = mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}});
        Matrix N = mat(new double[][]{{2, 2}});
        Matrix Z = mat(new double[][]{{1.0, 0.5}});
        Ret.pfqnMVACLD got = Pfqn_mvacld.pfqn_mvacld(L, N, Z, null);
        Ret.pfqnMVAC want = Pfqn_mvac.pfqn_mvac(L, N, Z);
        for (int r = 0; r < 2; r++) {
            assertEquals(want.X.get(0, r), got.X.get(0, r), EXACT_TOL);
            for (int i = 0; i < 2; i++) {
                assertEquals(want.Q.get(i, r), got.Q.get(i, r), EXACT_TOL);
            }
        }
    }

    // ---- oracle 2: genuine load dependence vs Pfqn_mvald ----

    @Test
    public void testMultiserverTwoServers() {
        assertMatchesMvald(mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}}),
                mat(new double[][]{{2, 2}}), mat(new double[][]{{1.0, 0.5}}), msRates(2, 4, 2));
    }

    /** S=D with load dependence: every chain resolved by a label interchange. */
    @Test
    public void testMultiserverNoDelay() {
        assertMatchesMvald(mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}}),
                mat(new double[][]{{2, 3}}), new Matrix(1, 2), msRates(2, 5, 3));
    }

    /** mu(j,n)=n turns a QLD center into an infinite server; c_i = 1 falls out of (21). */
    @Test
    public void testDelayEmulatedByRateN() {
        assertMatchesMvald(mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}}),
                mat(new double[][]{{2, 2}}), new Matrix(1, 2), msRates(2, 4, 99));
    }

    /** Nothing in (21)-(25) requires mu to be monotone or concave. */
    @Test
    public void testNonmonotoneRates() {
        int Nt = 4;
        Matrix mu = new Matrix(2, Nt);
        for (int i = 0; i < 2; i++) {
            for (int n = 1; n <= Nt; n++) {
                mu.set(i, n - 1, 1.0 + 0.5 * Math.sin(n));
            }
        }
        assertMatchesMvald(mat(new double[][]{{0.5, 0.2}, {0.3, 0.8}}),
                mat(new double[][]{{2, 2}}), mat(new double[][]{{1.0, 1.0}}), mu);
    }

    @Test
    public void testSingleClassMultiserver() {
        assertMatchesMvald(mat(new double[][]{{0.4}, {0.9}}), mat(new double[][]{{5}}),
                mat(new double[][]{{1.5}}), msRates(2, 5, 2));
    }

    @Test
    public void testThreeClassesMultiserver() {
        assertMatchesMvald(mat(new double[][]{{0.3, 0.5, 0.7}, {0.6, 0.1, 0.4}}),
                mat(new double[][]{{2, 1, 2}}), mat(new double[][]{{1.0, 0.0, 0.5}}),
                msRates(2, 5, 2));
    }

    /** A fixed-rate queue and a 2-server queue in the same network. */
    @Test
    public void testPerStationRatesDiffer() {
        int Nt = 4;
        Matrix mu = new Matrix(2, Nt);
        for (int n = 1; n <= Nt; n++) {
            mu.set(0, n - 1, 1.0);
            mu.set(1, n - 1, Math.min(n, 2));
        }
        assertMatchesMvald(mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}}),
                mat(new double[][]{{2, 2}}), mat(new double[][]{{1.0, 0.5}}), mu);
    }

    @Test
    public void testRandomizedAgainstMvald() {
        Random rng = new Random(5);
        int checked = 0;
        for (int tc = 0; tc < 30; tc++) {
            int M = 1 + rng.nextInt(2);
            int R = 1 + rng.nextInt(2);
            Matrix L = new Matrix(M, R);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    L.set(i, r, Math.round(10 * rng.nextDouble()) / 5.0);
                }
            }
            Matrix N = new Matrix(1, R);
            int tot = 0;
            for (int r = 0; r < R; r++) {
                int nr = rng.nextInt(3);
                N.set(0, r, nr);
                tot += nr;
            }
            Matrix Z = new Matrix(1, R);
            if (rng.nextBoolean()) {
                for (int r = 0; r < R; r++) {
                    Z.set(0, r, Math.round(10 * rng.nextDouble()) / 5.0);
                }
            }
            if (tot == 0) {
                continue;
            }
            boolean feasible = true;
            for (int r = 0; r < R; r++) {
                if (N.get(0, r) == 0) {
                    continue;
                }
                double s = Z.get(0, r);
                for (int i = 0; i < M; i++) {
                    s += L.get(i, r);
                }
                if (s == 0) {
                    feasible = false;
                }
            }
            if (!feasible) {
                continue;
            }
            Matrix mu = new Matrix(M, tot);
            for (int i = 0; i < M; i++) {
                int c = 1 + rng.nextInt(3);
                for (int n = 1; n <= tot; n++) {
                    mu.set(i, n - 1, Math.min(n, c));
                }
            }
            assertMatchesMvald(L, N, Z, mu);
            checked++;
        }
        assertTrue(checked > 10, "too few feasible random cases were generated");
    }

    /**
     * An empty class has no cycle time, reported as 0 as in Pfqn_dac (and Python
     * pfqn_mvald). MATLAB/JAR Pfqn_mvald instead leave an unguarded 0/0 = NaN there, so
     * this value is deliberately NOT taken from Pfqn_mvald.
     */
    @Test
    public void testEmptyClassHasNoCycleTime() {
        Ret.pfqnMVACLD got = Pfqn_mvacld.pfqn_mvacld(mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}}),
                mat(new double[][]{{3, 0}}), mat(new double[][]{{1.0, 2.0}}), msRates(2, 3, 2));
        assertEquals(0.0, got.X.get(0, 1), EXACT_TOL);
        assertEquals(0.0, got.R.get(0, 1), EXACT_TOL);
        assertTrue(got.R.get(0, 0) > 0, "the populated class must have a positive cycle time");
    }

    /**
     * Regression: Pfqn_mvald left an empty class's cycle time as an unguarded N/XN =
     * 0/0, so MATLAB/JAR returned NaN there (Python already guarded), which propagated
     * to callers such as Pfqn_mvams. An empty class has zero throughput, hence no cycle
     * time: all three codebases now report 0, as Pfqn_dac and Pfqn_mvacld do.
     */
    @Test
    public void testMvaldEmptyClassHasNoCycleTime() {
        Matrix L = mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}});
        Matrix N = mat(new double[][]{{3, 0}});
        Matrix Z = mat(new double[][]{{1.0, 2.0}});
        Matrix mu = msRates(2, 3, 2);
        Ret.pfqnMVALD ref = Pfqn_mvald.pfqn_mvald(L, N, Z, mu);
        assertTrue(!Double.isNaN(ref.R.get(0, 1)),
                "an empty class must not yield a NaN cycle time");
        assertEquals(0.0, ref.R.get(0, 1), EXACT_TOL);
        assertTrue(ref.R.get(0, 0) > 0);
        // and it must agree with Pfqn_mvacld on every class, empty ones included
        Ret.pfqnMVACLD got = Pfqn_mvacld.pfqn_mvacld(L, N, Z, mu);
        for (int r = 0; r < 2; r++) {
            assertEquals(ref.R.get(0, r), got.R.get(0, r), EXACT_TOL, "C class " + r);
        }
    }

    // ---- input validation ----

    @Test
    public void testRejectsBadRateMatrix() {
        Matrix L = mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}});
        Matrix N = mat(new double[][]{{2, 2}});
        assertThrows(() -> Pfqn_mvacld.pfqn_mvacld(L, N, null, ones(3, 4)));
        assertThrows(() -> Pfqn_mvacld.pfqn_mvacld(L, N, null, ones(2, 2)));
        assertThrows(() -> Pfqn_mvacld.pfqn_mvacld(L, N, null, new Matrix(2, 4)));
    }

    private static void assertThrows(Runnable r) {
        try {
            r.run();
        } catch (IllegalArgumentException e) {
            return;
        }
        throw new AssertionError("expected an IllegalArgumentException");
    }

    // ---- helpers ----

    /** With mu = 1 the QLD recursion must reproduce the SSFR/IS recursion exactly. */
    private static void assertReducesToMvac(Matrix L, Matrix N, Matrix Z) {
        int Nt = 0;
        for (int r = 0; r < N.getNumCols(); r++) {
            Nt += (int) N.get(0, r);
        }
        Ret.pfqnMVACLD got = Pfqn_mvacld.pfqn_mvacld(L, N, Z, ones(L.getNumRows(), Nt));
        Ret.pfqnMVAC want = Pfqn_mvac.pfqn_mvac(L, N, Z);
        for (int r = 0; r < L.getNumCols(); r++) {
            assertEquals(want.X.get(0, r), got.X.get(0, r), EXACT_TOL, "X class " + r);
            for (int i = 0; i < L.getNumRows(); i++) {
                assertEquals(want.Q.get(i, r), got.Q.get(i, r), EXACT_TOL,
                        "Q station " + i + " class " + r);
                // C is NOT compared: Pfqn_mvac returns an (M x R) residence time
                // while the LD family returns a (1 x R) cycle time
            }
        }
        checkInternalOracles(got, L, N, Z);
    }

    /** MVAC-LD and the load-dependent population recursion are both exact. */
    private static void assertMatchesMvald(Matrix L, Matrix N, Matrix Z, Matrix mu) {
        Ret.pfqnMVACLD got = Pfqn_mvacld.pfqn_mvacld(L, N, Z, mu);
        Matrix Zagg = new Matrix(1, Z.getNumCols());
        for (int r = 0; r < Z.getNumCols(); r++) {
            double s = 0;
            for (int i = 0; i < Z.getNumRows(); i++) {
                s += Z.get(i, r);
            }
            Zagg.set(0, r, s);
        }
        Ret.pfqnMVALD want = Pfqn_mvald.pfqn_mvald(L, N, Zagg, mu);
        // every output must match, not just X and Q: the LD family contract is
        // X (1xR), Q (MxR), U (Mx1) = 1-P_j(0), R (1xR) cycle time
        assertEquals(want.U.getNumRows(), got.U.getNumRows(), "U rows");
        assertEquals(want.U.getNumCols(), got.U.getNumCols(), "U cols");
        assertEquals(want.R.getNumRows(), got.R.getNumRows(), "C rows");
        assertEquals(want.R.getNumCols(), got.R.getNumCols(), "C cols");
        for (int r = 0; r < L.getNumCols(); r++) {
            assertEquals(want.X.get(0, r), got.X.get(0, r), EXACT_TOL, "X class " + r);
            if (N.get(0, r) > 0) {
                // Only populated classes: Pfqn_mvald leaves an empty class's cycle time
                // as an unguarded 0/0 = NaN, whereas Pfqn_dac (and this routine) report
                // 0. That divergence is Pfqn_mvald's; see the empty-class test below.
                assertEquals(want.R.get(0, r), got.R.get(0, r), EXACT_TOL, "C class " + r);
            }
            for (int i = 0; i < L.getNumRows(); i++) {
                assertEquals(want.Q.get(i, r), got.Q.get(i, r), EXACT_TOL,
                        "Q station " + i + " class " + r);
            }
        }
        for (int i = 0; i < L.getNumRows(); i++) {
            assertEquals(want.U.get(i, 0), got.U.get(i, 0), EXACT_TOL, "U station " + i);
        }
        checkInternalOracles(got, L, N, Z);
    }

    private static void checkInternalOracles(Ret.pfqnMVACLD got, Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        int Nt = 0;
        for (int r = 0; r < R; r++) {
            Nt += (int) N.get(0, r);
        }
        for (int i = 0; i < M; i++) {
            // (3) eq. (25) is self-normalizing
            double s = 0;
            double mean = 0;
            for (int n = 0; n <= Nt; n++) {
                s += got.pij.get(i, n);
                mean += n * got.pij.get(i, n);
            }
            assertEquals(1.0, s, EXACT_TOL, "marginal of station " + i + " must sum to 1");
            // (4) the marginals reproduce the means
            double q = 0;
            for (int r = 0; r < R; r++) {
                q += got.Q.get(i, r);
            }
            assertEquals(q, mean, EXACT_TOL, "Q must equal E[n] at station " + i);
            // utilization is 1-P_j(0), per station
            assertEquals(1.0 - got.pij.get(i, 0), got.U.get(i, 0), EXACT_TOL,
                    "U station " + i);
        }
        // every customer is accounted for
        double tot = 0;
        for (int r = 0; r < R; r++) {
            for (int i = 0; i < Z.getNumRows(); i++) {
                tot += got.X.get(0, r) * Z.get(i, r);
            }
            for (int i = 0; i < M; i++) {
                tot += got.Q.get(i, r);
            }
        }
        assertEquals(Nt, tot, 1e-8, "population conservation");
    }
}
