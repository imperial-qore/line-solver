/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.mva;

import java.util.Random;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the MVAC (mean value analysis by chain) algorithm of Conway, de Souza e
 * Silva and Lavenberg, IEEE Trans. Computers 38(3):432-442, 1989.
 *
 * <p>MVAC is exact, so it is checked against two independent references: the worked
 * example of Section III of the paper, whose closed-form values are quoted there, and the
 * classic population-recursion MVA of {@link Pfqn_mva}, which computes the same measures
 * by a structurally unrelated recursion. The three chain-resolution paths of the algorithm
 * (part 1 for chain K, part 2 for the chains visiting an IS center, and the label
 * interchanges for the chains visiting only SSFR centers) are covered explicitly, since a
 * network exercising only one of them would leave the others untested.</p>
 */
public class Pfqn_mvacTest {

    private static final double EXACT_TOL = 1e-10;

    private static Matrix mat(double[][] v) {
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                m.set(i, j, v[i][j]);
            }
        }
        return m;
    }

    /**
     * Section III example: J = 2 SSFR centers, K = 2 single-customer chains,
     * a11 = 1, a21 = 2, a12 = 2, a22 = 3. The paper reports lambda_2 = 3/23,
     * L_12 = 8/23, L_22 = 15/23, L_1 = 15/23 and L_2 = 31/23.
     */
    @Test
    public void testPaperSectionIIIExample() {
        Matrix L = mat(new double[][]{{1, 2}, {2, 3}});
        Matrix N = mat(new double[][]{{1, 1}});
        Ret.pfqnMVAC ret = Pfqn_mvac.pfqn_mvac(L, N);

        assertEquals(3.0 / 23.0, ret.X.get(0, 1), EXACT_TOL);
        assertEquals(8.0 / 23.0, ret.Q.get(0, 1), EXACT_TOL);
        assertEquals(15.0 / 23.0, ret.Q.get(1, 1), EXACT_TOL);
        // total mean number of customers at each center, chains 1 and 2 summed
        assertEquals(15.0 / 23.0, ret.Q.get(0, 0) + ret.Q.get(0, 1), EXACT_TOL);
        assertEquals(31.0 / 23.0, ret.Q.get(1, 0) + ret.Q.get(1, 1), EXACT_TOL);
    }

    /** Every chain visits an IS center (S = 0): all chains but K are resolved by part 2. */
    @Test
    public void testAllChainsVisitDelay() {
        assertAgreesWithMva(mat(new double[][]{{0.4, 0.7, 0.2}, {0.9, 0.3, 0.5}}),
                mat(new double[][]{{1, 2, 1}}), mat(new double[][]{{1.0, 0.5, 2.0}}));
    }

    /** No IS center (S = D): every chain is resolved by a label interchange. */
    @Test
    public void testNoDelay() {
        assertAgreesWithMva(mat(new double[][]{{0.4, 0.7, 0.2}, {0.9, 0.3, 0.5}}),
                mat(new double[][]{{1, 2, 1}}), new Matrix(1, 3));
    }

    /** Mixed network (0 &lt; S &lt; D): part 2 and the label interchanges are both used. */
    @Test
    public void testMixedDelayAndQueueOnlyChains() {
        assertAgreesWithMva(mat(new double[][]{{0.4, 0.7, 0.2}, {0.9, 0.3, 0.5}}),
                mat(new double[][]{{2, 1, 2}}), mat(new double[][]{{1.0, 0.0, 0.5}}));
    }

    /** Classes with identical demand columns collapse into one subset of identical chains. */
    @Test
    public void testIdenticalClassesCollapse() {
        assertAgreesWithMva(mat(new double[][]{{0.4, 0.4, 0.2}, {0.9, 0.9, 0.5}}),
                mat(new double[][]{{2, 3, 1}}), mat(new double[][]{{1.0, 1.0, 0.0}}));
    }

    /** A single class with several customers, i.e. D = 1 and no label interchange. */
    @Test
    public void testSingleClass() {
        assertAgreesWithMva(mat(new double[][]{{0.4}, {0.9}}),
                mat(new double[][]{{6}}), mat(new double[][]{{1.5}}));
    }

    /** More than one IS center, each carried as its own cell of the multiplicity vector. */
    @Test
    public void testMultipleDelayCenters() {
        assertAgreesWithMva(mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}}),
                mat(new double[][]{{2, 2}}), mat(new double[][]{{1.0, 0.0}, {0.5, 2.0}}));
    }

    /** An empty class must carry no throughput and no customers. */
    @Test
    public void testEmptyClass() {
        Matrix L = mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}});
        Ret.pfqnMVAC ret = Pfqn_mvac.pfqn_mvac(L, mat(new double[][]{{3, 0}}),
                mat(new double[][]{{1.0, 2.0}}));
        assertEquals(0.0, ret.X.get(0, 1), EXACT_TOL);
        assertEquals(0.0, ret.Q.get(0, 1), EXACT_TOL);
        assertEquals(0.0, ret.Q.get(1, 1), EXACT_TOL);
        assertAgreesWithMva(L, mat(new double[][]{{3, 0}}), mat(new double[][]{{1.0, 2.0}}));
    }

    /** The regime MVAC targets: few centers, many chains. */
    @Test
    public void testManyChainsFewCenters() {
        assertAgreesWithMva(
                mat(new double[][]{{0.3, 0.5, 0.7, 0.2, 0.9}, {0.6, 0.1, 0.4, 0.8, 0.3}}),
                mat(new double[][]{{2, 2, 2, 2, 2}}),
                mat(new double[][]{{1.0, 0.0, 2.0, 0.0, 0.5}}));
    }

    /** Randomized cross-check against the population-recursion MVA. */
    @Test
    public void testRandomizedAgainstMva() {
        Random rng = new Random(42);
        int checked = 0;
        for (int tc = 0; tc < 60; tc++) {
            int M = 1 + rng.nextInt(3);
            int R = 1 + rng.nextInt(3);
            Matrix L = new Matrix(M, R);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    L.set(i, r, Math.round(10 * rng.nextDouble()) / 5.0);
                }
            }
            Matrix N = new Matrix(1, R);
            int tot = 0;
            for (int r = 0; r < R; r++) {
                int nr = rng.nextInt(4);
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
            // every populated class must have a nonzero demand somewhere
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
            assertAgreesWithMva(L, N, Z);
            checked++;
        }
        assertTrue(checked > 20, "too few feasible random cases were generated");
    }

    /**
     * MVAC and MVA are both exact, so they must agree to numerical precision on the
     * throughput, queue-length and utilization of every class.
     */
    private static void assertAgreesWithMva(Matrix L, Matrix N, Matrix Z) {
        Ret.pfqnMVAC got = Pfqn_mvac.pfqn_mvac(L, N, Z);
        // Pfqn_mva takes a single aggregated delay, so the IS centers are summed
        Matrix Zagg = new Matrix(1, Z.getNumCols());
        for (int r = 0; r < Z.getNumCols(); r++) {
            double s = 0;
            for (int i = 0; i < Z.getNumRows(); i++) {
                s += Z.get(i, r);
            }
            Zagg.set(0, r, s);
        }
        Ret.pfqnMVA want = Pfqn_mva.pfqn_mva(L, N, Zagg, null);
        for (int r = 0; r < L.getNumCols(); r++) {
            assertEquals(want.X.get(0, r), got.X.get(0, r), EXACT_TOL, "X class " + r);
            for (int i = 0; i < L.getNumRows(); i++) {
                assertEquals(want.Q.get(i, r), got.Q.get(i, r), EXACT_TOL,
                        "Q station " + i + " class " + r);
                assertEquals(want.U.get(i, r), got.U.get(i, r), EXACT_TOL,
                        "U station " + i + " class " + r);
            }
        }
        // the customers must all be accounted for
        double tot = 0;
        for (int r = 0; r < L.getNumCols(); r++) {
            tot += got.X.get(0, r) * Zagg.get(0, r);
            for (int i = 0; i < L.getNumRows(); i++) {
                tot += got.Q.get(i, r);
            }
        }
        double want_tot = 0;
        for (int r = 0; r < L.getNumCols(); r++) {
            want_tot += N.get(0, r);
        }
        assertEquals(want_tot, tot, 1e-8, "population conservation");
    }
}
