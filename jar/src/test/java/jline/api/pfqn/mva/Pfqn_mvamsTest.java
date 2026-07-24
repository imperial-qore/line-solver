/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.mva;

import jline.api.pfqn.ld.Pfqn_mvald;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Contract tests for {@link Pfqn_mvams}, the general-purpose MVA dispatcher for mixed
 * networks with multiserver nodes.
 *
 * <p>Pfqn_mvams routes to one of four algorithms depending on whether the model has open
 * classes and whether any station is multiserver:</p>
 *
 * <pre>
 *   single-server + closed  -&gt; Pfqn_mva
 *   single-server + mixed   -&gt; Pfqn_mvamx
 *   multiserver   + mixed   -&gt; Pfqn_mvaldms
 *   multiserver   + closed  -&gt; Pfqn_mvald
 * </pre>
 *
 * <p>Those four do NOT share a return contract of their own: the load-dependent family
 * ({@link Pfqn_mvald}) reports a (1 x R) cycle time. Pfqn_mvams must nonetheless present
 * ONE contract on every branch, namely the {@link Ret.pfqnMVA} one it declares: an
 * (M x R) per-station residence time. These tests pin that, because a branch returning a
 * different QUANTITY under the same name is invisible to any test that exercises only one
 * branch.</p>
 */
public class Pfqn_mvamsTest {

    private static final double TOL = 1e-9;

    private static Matrix mat(double[][] v) {
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                m.set(i, j, v[i][j]);
            }
        }
        return m;
    }

    private static final Matrix L = mat(new double[][]{{0.4, 0.7}, {0.9, 0.3}});
    private static final Matrix Z = mat(new double[][]{{1.0, 0.5}});

    /** CN must be an (M x R) residence time satisfying QN = XN*CN, on every branch. */
    private void assertResidenceTimeContract(String name, double[] lam, double[] N,
                                             double[] S) {
        Ret.pfqnMVA got = Pfqn_mvams.pfqn_mvams(mat(new double[][]{lam}), L,
                mat(new double[][]{N}), Z, Matrix.ones(2, 1), mat(new double[][]{S}).transpose());
        assertEquals(L.getNumRows(), got.R.getNumRows(), "CN rows on branch: " + name);
        assertEquals(L.getNumCols(), got.R.getNumCols(), "CN cols on branch: " + name);
        for (int r = 0; r < L.getNumCols(); r++) {
            for (int i = 0; i < L.getNumRows(); i++) {
                assertTrue(!Double.isNaN(got.R.get(i, r)), "CN NaN on branch: " + name);
            }
            if (Double.isInfinite(N[r]) || N[r] == 0) {
                continue;
            }
            for (int i = 0; i < L.getNumRows(); i++) {
                // the defining identity of a residence time
                assertEquals(got.Q.get(i, r), got.X.get(0, r) * got.R.get(i, r), TOL,
                        "QN != XN*CN on branch: " + name);
            }
        }
    }

    @Test
    public void testClosedSingleServerBranch() {
        assertResidenceTimeContract("closed single-server (Pfqn_mva)",
                new double[]{0, 0}, new double[]{2, 2}, new double[]{1, 1});
    }

    @Test
    public void testClosedMultiserverBranch() {
        assertResidenceTimeContract("closed multiserver (Pfqn_mvald)",
                new double[]{0, 0}, new double[]{2, 2}, new double[]{2, 2});
    }

    @Test
    public void testMixedSingleServerBranch() {
        assertResidenceTimeContract("mixed single-server (Pfqn_mvamx)",
                new double[]{0.2, 0}, new double[]{Double.POSITIVE_INFINITY, 2},
                new double[]{1, 1});
    }

    @Test
    public void testMixedMultiserverBranch() {
        assertResidenceTimeContract("mixed multiserver (Pfqn_mvaldms)",
                new double[]{0.2, 0}, new double[]{Double.POSITIVE_INFINITY, 2},
                new double[]{2, 2});
    }

    @Test
    public void testClosedMultiserverEmptyClass() {
        assertResidenceTimeContract("closed multiserver, empty class",
                new double[]{0, 0}, new double[]{3, 0}, new double[]{2, 2});
    }

    @Test
    public void testClosedMultiserverThreeServers() {
        assertResidenceTimeContract("closed multiserver, 3 servers",
                new double[]{0, 0}, new double[]{3, 2}, new double[]{3, 2});
    }

    /** With S=1 the dispatcher takes the Pfqn_mva branch and must agree with it. */
    @Test
    public void testSingleServerLimitEqualsMva() {
        Matrix N = mat(new double[][]{{2, 2}});
        Ret.pfqnMVA got = Pfqn_mvams.pfqn_mvams(new Matrix(1, 2), L, N, Z,
                Matrix.ones(2, 1), mat(new double[][]{{1}, {1}}));
        Ret.pfqnMVA want = Pfqn_mva.pfqn_mva(L, N, Z, null);
        for (int r = 0; r < 2; r++) {
            assertEquals(want.X.get(0, r), got.X.get(0, r), TOL, "X class " + r);
            for (int i = 0; i < 2; i++) {
                assertEquals(want.Q.get(i, r), got.Q.get(i, r), TOL, "Q " + i + "," + r);
                assertEquals(want.R.get(i, r), got.R.get(i, r), TOL, "C " + i + "," + r);
            }
        }
    }

    /** The fix must not disturb X and Q, which come straight from Pfqn_mvald. */
    @Test
    public void testMultiserverBranchStillMatchesMvald() {
        Matrix N = mat(new double[][]{{2, 2}});
        int Nt = 4;
        Matrix mu = new Matrix(2, Nt);
        for (int i = 0; i < 2; i++) {
            for (int n = 1; n <= Nt; n++) {
                mu.set(i, n - 1, Math.min(n, 2));
            }
        }
        Ret.pfqnMVA got = Pfqn_mvams.pfqn_mvams(new Matrix(1, 2), L, N, Z,
                Matrix.ones(2, 1), mat(new double[][]{{2}, {2}}));
        Ret.pfqnMVALD want = Pfqn_mvald.pfqn_mvald(L, N, Z, mu);
        for (int r = 0; r < 2; r++) {
            assertEquals(want.X.get(0, r), got.X.get(0, r), TOL, "X class " + r);
            double sumC = 0;
            for (int i = 0; i < 2; i++) {
                assertEquals(want.Q.get(i, r), got.Q.get(i, r), TOL, "Q " + i + "," + r);
                sumC += got.R.get(i, r);
            }
            // the residence times must decompose Pfqn_mvald's cycle time
            assertEquals(want.R.get(0, r), sumC, TOL, "sum_i C(i,r) != cycle time, class " + r);
        }
    }
}
