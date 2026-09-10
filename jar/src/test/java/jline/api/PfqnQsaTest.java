package jline.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.mva.Pfqn_mva;
import jline.api.pfqn.mva.Pfqn_qsa;
import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

/**
 * Queue-Shift Approximation (Pfqn_qsa).
 *
 * The three stress networks of Sect. 5.2 of Schweitzer, Serazzi and Broglia
 * (Tools'98, LNCS 1469) are fully specified in print together with the QSA
 * errors they produce under the paper's own metric eq. (17), so they pin the
 * algorithm rather than this port: any drift in the shift definition or in the
 * extrapolation (15) moves them.
 */
public class PfqnQsaTest {

    private static Matrix mat(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                m.set(i, j, a[i][j]);
            }
        }
        return m;
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    private static final double[][] TABLE2 = {{16, 50}, {19, 28}, {0, 42}, {12, 25}, {39, 29}};
    private static final double[][] TABLE3 = {{16, 73}, {59, 15}, {6, 26}, {2, 36}, {10, 0}};
    private static final double[][] TABLE4 = {{10, 1, 1}, {1, 10, 1}, {1, 1, 10}};

    /** eq. (17): err(Q) = max |Q_appr - Q_MVA| / K_r, err(U) = max |U_appr - U_MVA|. */
    private static double[] eq17(double[][] demands, double[] pop) {
        Matrix L = mat(demands);
        Matrix N = row(pop);
        int M = L.getNumRows();
        int R = L.getNumCols();
        Matrix Z = new Matrix(1, R);
        Ret.pfqnMVA ex = Pfqn_mva.pfqn_mva(L, N, Z, Matrix.ones(1, M));
        Ret.pfqnAMVA q = Pfqn_qsa.pfqn_qsa(L, N, Z);
        double eQ = 0.0;
        double eU = 0.0;
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                eQ = Math.max(eQ, Math.abs(q.Q.get(i, r) - ex.Q.get(i, r)) / pop[r]);
                eU = Math.max(eU, Math.abs(q.U.get(i, r) - ex.U.get(i, r)));
            }
        }
        return new double[]{eQ, eU};
    }

    @Test
    public void table2StressCaseMatchesThePublishedErrors() {
        // K = (21,29), the case with the largest Linearizer error.
        // Published: QSA err(Q) = 0.001812, err(U) = 0.000639.
        double[] e = eq17(TABLE2, new double[]{21, 29});
        assertEquals(0.001812, e[0], 1e-6);
        assertEquals(0.000639, e[1], 1e-6);
    }

    @Test
    public void table3StressCaseIsExactToMachinePrecision() {
        // K = (111,89), maximizing the Linearizer-to-QSA error ratio. The
        // published 4.18e-10 / 5.54e-10 is that paper's own residual tolerance;
        // the Newton solve here reaches machine precision, consistent with its
        // claim of "errors of zero (at least 6 digits)".
        double[] e = eq17(TABLE3, new double[]{111, 89});
        assertTrue(e[0] < 4.18e-10, "err(Q) = " + e[0]);
        assertTrue(e[1] < 5.54e-10, "err(U) = " + e[1]);
    }

    @Test
    public void table4ThreeClassStressCaseMatchesThePublishedErrors() {
        // K = (2,2,2) on the Chandy-Neuse Example 2 loadings. Published:
        // err(Q) = 0.000185, err(U) = 0.000415. The tiny population is what
        // makes this discriminating: (15) has no 1/Ksum to hide in.
        double[] e = eq17(TABLE4, new double[]{2, 2, 2});
        assertEquals(0.000185, e[0], 1e-6);
        assertEquals(0.000415, e[1], 1e-6);
    }

    @Test
    public void qsaBeatsLinearizerOnTable4() {
        // Sect. 5.2: on this model QSA has roughly half the Linearizer error.
        double[] e = eq17(TABLE4, new double[]{2, 2, 2});
        assertTrue(e[0] < 0.000318, "err(Q) = " + e[0]);   // published LIN error
        assertTrue(e[1] < 0.000716, "err(U) = " + e[1]);
    }

    @Test
    public void populationIsConservedAndQueuesStayNonNegative() {
        Matrix L = mat(TABLE2);
        Matrix N = row(21, 29);
        Ret.pfqnAMVA q = Pfqn_qsa.pfqn_qsa(L, N, new Matrix(1, 2));
        double total = 0.0;
        for (int r = 0; r < 2; r++) {
            double perClass = 0.0;
            for (int i = 0; i < 5; i++) {
                assertTrue(q.Q.get(i, r) >= -1e-12, "negative queue length");
                perClass += q.Q.get(i, r);
            }
            assertEquals(N.get(0, r), perClass, 1e-10);   // eq. (6)
            total += perClass;
        }
        assertEquals(50.0, total, 1e-10);
    }

    @Test
    public void thinkTimeIsHeldOutsideTheQueues() {
        Matrix L = mat(new double[][]{{1.0, 0.5}, {0.7, 1.2}, {0.3, 0.9}});
        Matrix N = row(6, 4);
        Matrix Z = row(2, 1);
        Ret.pfqnAMVA q = Pfqn_qsa.pfqn_qsa(L, N, Z);
        double held = 0.0;
        for (int i = 0; i < 3; i++) {
            for (int r = 0; r < 2; r++) {
                held += q.Q.get(i, r);
            }
        }
        for (int r = 0; r < 2; r++) {
            held += q.X.get(0, r) * Z.get(0, r);
        }
        assertEquals(10.0, held, 1e-10);
        Ret.pfqnMVA ex = Pfqn_mva.pfqn_mva(L, N, Z, Matrix.ones(1, 3));
        for (int i = 0; i < 3; i++) {
            for (int r = 0; r < 2; r++) {
                assertTrue(Math.abs(q.Q.get(i, r) - ex.Q.get(i, r)) < 2e-2);
            }
        }
    }

    @Test
    public void delayCentreDeclaredThroughTypeHasNoQueueingTerm() {
        Matrix L = mat(new double[][]{{1.0, 0.5}, {0.7, 1.2}, {0.3, 0.9}});
        Matrix N = row(6, 4);
        SchedStrategy[] type = {SchedStrategy.PS, SchedStrategy.INF, SchedStrategy.PS};
        Ret.pfqnAMVA q = Pfqn_qsa.pfqn_qsa(L, N, new Matrix(1, 2), type, 1e-10, 100);
        for (int r = 0; r < 2; r++) {
            assertEquals(q.X.get(0, r) * L.get(1, r), q.Q.get(1, r), 1e-12);   // eq. (13b)
            assertEquals(L.get(1, r), q.R.get(1, r), 1e-12);
        }
    }

    @Test
    public void twoLevelVariantIsTheCruderOne() {
        // Sect. 3.1 gives eq. (14) an error of O(1/Ksum); on Table 2
        // (Ksum = 50) that is the order observed.
        Matrix L = mat(TABLE2);
        Matrix N = row(21, 29);
        Matrix Z = new Matrix(1, 2);
        Ret.pfqnMVA ex = Pfqn_mva.pfqn_mva(L, N, Z, Matrix.ones(1, 5));
        Ret.pfqnAMVA q2 = Pfqn_qsa.pfqn_qsa(L, N, Z, null, 1e-10, 100, 2, null);
        Ret.pfqnAMVA q3 = Pfqn_qsa.pfqn_qsa(L, N, Z, null, 1e-10, 100, 3, null);
        double e2 = 0.0;
        double e3 = 0.0;
        for (int i = 0; i < 5; i++) {
            for (int r = 0; r < 2; r++) {
                e2 = Math.max(e2, Math.abs(q2.Q.get(i, r) - ex.Q.get(i, r)) / N.get(0, r));
                e3 = Math.max(e3, Math.abs(q3.Q.get(i, r) - ex.Q.get(i, r)) / N.get(0, r));
            }
        }
        assertTrue(e3 < e2, "three-level " + e3 + " should beat two-level " + e2);
        assertTrue(e2 > 1.0 / 50.0 && e2 < 0.1, "two-level error out of the O(1/Ksum) band: " + e2);
    }

    @Test
    public void exactWhenOneJobPerClass() {
        // Every arrival sees the other classes only, so the arrival-instant
        // queue length is exact and so is the shift.
        Matrix L = mat(new double[][]{{1.0, 2.0}, {3.0, 1.0}, {2.0, 2.0}});
        Matrix N = row(1, 1);
        Matrix Z = new Matrix(1, 2);
        Ret.pfqnMVA ex = Pfqn_mva.pfqn_mva(L, N, Z, Matrix.ones(1, 3));
        Ret.pfqnAMVA q = Pfqn_qsa.pfqn_qsa(L, N, Z);
        for (int i = 0; i < 3; i++) {
            for (int r = 0; r < 2; r++) {
                assertEquals(ex.Q.get(i, r), q.Q.get(i, r), 1e-10);
            }
        }
    }

    @Test
    public void emptyClassContributesNothing() {
        Matrix L = mat(new double[][]{{1.0, 2.0}, {3.0, 1.0}, {2.0, 2.0}});
        Matrix N = row(3, 0);
        Ret.pfqnAMVA q = Pfqn_qsa.pfqn_qsa(L, N, new Matrix(1, 2));
        double sum = 0.0;
        for (int i = 0; i < 3; i++) {
            assertEquals(0.0, q.Q.get(i, 1), 0.0);
            sum += q.Q.get(i, 0);
        }
        assertEquals(0.0, q.X.get(0, 1), 0.0);
        assertEquals(3.0, sum, 1e-10);
    }
}
