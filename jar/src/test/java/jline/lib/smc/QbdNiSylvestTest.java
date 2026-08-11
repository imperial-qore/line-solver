package jline.lib.smc;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Map;
import java.util.Random;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;

/**
 * Regression tests for the Sylvester-based Newton iteration of the SMC solver.
 *
 * These guard three defects that made the whole QBD_NI/'Sylvest' path return
 * either garbage or an empty result: QBD_NI_Sylvest did not solve AXB+CX=D at
 * all, QBD_NI returned early exactly when the Newton iteration was needed, and
 * it overwrote the iterated R instead of deriving G from it. Expected values
 * are cross-checked against MATLAB QBD_NI_Sylvest.m / QBD_NI.m.
 */
public class QbdNiSylvestTest {

    private static Matrix mk(double[][] v) {
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++) {
                m.set(i, j, v[i][j]);
            }
        }
        return m;
    }

    /** The solver must actually satisfy its defining equation AXB + CX = D. */
    @Test
    public void sylvesterResidualIsNegligible() {
        Random rnd = new Random(11);
        int[] orders = {2, 3, 5, 8};
        for (int idx = 0; idx < orders.length; idx++) {
            int n = orders[idx];
            Matrix A = new Matrix(n, n);
            Matrix B = new Matrix(n, n);
            Matrix C = new Matrix(n, n);
            Matrix D = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    A.set(i, j, rnd.nextGaussian());
                    B.set(i, j, rnd.nextGaussian());
                    C.set(i, j, rnd.nextGaussian());
                    D.set(i, j, rnd.nextGaussian());
                }
            }
            Matrix X = QBD_NI_Sylvest.QBD_NI_Sylvest(A, B, C, D);
            Matrix res = A.mult(X).mult(B).add(1.0, C.mult(X)).sub(D);
            double scale = Math.max(1.0, Matrix.infNorm(X));
            assertTrue(Matrix.infNorm(res) < 1e-9 * scale,
                    "AXB+CX=D residual too large at order " + n + ": " + Matrix.infNorm(res));
        }
    }

    /** Ground truth from MATLAB QBD_NI_Sylvest(A,B,C,D). */
    @Test
    public void sylvesterMatchesMatlab() {
        Matrix A = mk(new double[][]{{0.8, 0.3, 0.1}, {0.2, 0.9, 0.4}, {0.5, 0.1, 0.7}});
        Matrix B = mk(new double[][]{{0.6, 0.2, 0.3}, {0.1, 0.5, 0.2}, {0.4, 0.3, 0.9}});
        Matrix C = mk(new double[][]{{-1.7, 0.4, 0.2}, {0.3, -1.9, 0.5}, {0.6, 0.2, -2.1}});
        Matrix D = mk(new double[][]{{1.0, 0.5, 0.25}, {0.75, 1.5, 0.5}, {0.3, 0.6, 1.2}});

        double[][] expected = {
            {0.827746630775, 0.541810047547, 2.127534745400},
            {1.180087678673, 0.194192059739, 2.183084850326},
            {1.036512776709, 0.488743672148, 1.476685277143}};

        Matrix X = QBD_NI_Sylvest.QBD_NI_Sylvest(A, B, C, D);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                assertEquals(expected[i][j], X.get(i, j), 1e-10,
                        "X[" + i + "][" + j + "] mismatch");
            }
        }
    }

    /**
     * The 'Sylvest' Newton path must reproduce the R and G of the reference
     * cyclic-reduction solver; both are cross-checked against MATLAB QBD_NI.
     */
    @Test
    public void qbdNiSylvestMatchesCyclicReduction() {
        Matrix A0 = mk(new double[][]{{0.2, 0.0}, {0.0, 0.1}});
        Matrix A1 = mk(new double[][]{{0.3, 0.1}, {0.1, 0.4}});
        Matrix A2 = mk(new double[][]{{0.3, 0.1}, {0.2, 0.2}});

        Map<String, Matrix> ni = QBD_NI.QBD_NI(A0, A1, A2, null, null, null, null);
        Matrix R = ni.get("R");
        Matrix G = ni.get("G");

        assertEquals(2, R.getNumRows(), "QBD_NI must return a populated R");
        assertEquals(2, R.getNumCols(), "QBD_NI must return a populated R");

        // MATLAB QBD_NI(A0,A1,A2)
        double[][] expectedR = {{0.648183878790, 0.339099965932}, {0.527724181815, 0.491350051101}};
        double[][] expectedG = {{0.384321787882, 0.046712470191}, {0.143402393933, 0.198962555360}};
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(expectedR[i][j], R.get(i, j), 1e-10, "R[" + i + "][" + j + "] mismatch");
                assertEquals(expectedG[i][j], G.get(i, j), 1e-10, "G[" + i + "][" + j + "] mismatch");
            }
        }

        // R must satisfy the QBD equation A2 + R*A1 + R^2*A0 = R.
        Matrix res = A2.add(1.0, R.mult(A1)).add(1.0, R.mult(R).mult(A0)).sub(R);
        assertTrue(Matrix.infNorm(res) < 1e-10, "QBD R residual too large: " + Matrix.infNorm(res));

        Map<String, Matrix> cr = QBD_CR.QBD_CR(A0, A1, A2, null, null, null, null);
        assertTrue(Matrix.infNorm(R.sub(cr.get("R"))) < 1e-10, "QBD_NI R disagrees with QBD_CR");
    }

    /**
     * QBD_NI must not rescale the caller's blocks in place. MATLAB rebinds
     * (A0=A0/lamb) and QBD_CR copies its arguments, so a caller reusing its
     * blocks after the call must see them unchanged. Uses a continuous-time QBD
     * (negative diagonal), which is what triggers the uniformization rescale.
     */
    @Test
    public void qbdNiDoesNotMutateItsInputs() {
        Matrix A0 = mk(new double[][]{{0.4, 0.0}, {0.0, 0.4}});
        Matrix A1 = mk(new double[][]{{-1.5, 0.1}, {0.1, -1.5}});
        Matrix A2 = mk(new double[][]{{1.0, 0.0}, {0.0, 1.0}});
        Matrix A0ref = A0.copy();
        Matrix A1ref = A1.copy();
        Matrix A2ref = A2.copy();

        QBD_NI.QBD_NI(A0, A1, A2, null, null, null, null);

        assertEquals(0.0, Matrix.infNorm(A0.sub(A0ref)), 1e-12, "QBD_NI modified its A0 argument");
        assertEquals(0.0, Matrix.infNorm(A1.sub(A1ref)), 1e-12, "QBD_NI modified its A1 argument");
        assertEquals(0.0, Matrix.infNorm(A2.sub(A2ref)), 1e-12, "QBD_NI modified its A2 argument");
    }

    /** W*A*V=LBAR (Hessenberg) and W*B*V=NBAR (triangular) with W,V orthogonal. */
    @Test
    public void generalizedHessenbergIdentitiesHold() {
        Random rnd = new Random(7);
        int[] orders = {2, 3, 5, 8};
        for (int idx = 0; idx < orders.length; idx++) {
            int n = orders[idx];
            Matrix A = new Matrix(n, n);
            Matrix B = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    A.set(i, j, rnd.nextGaussian());
                    B.set(i, j, rnd.nextGaussian());
                }
            }
            Map<String, Matrix> h = A.hess(B);
            Matrix lbar = h.get("LBAR");
            Matrix nbar = h.get("NBAR");
            Matrix w = h.get("W");
            Matrix v = h.get("V");

            assertTrue(Matrix.infNorm(w.mult(A).mult(v).sub(lbar)) < 1e-10, "W*A*V != LBAR at order " + n);
            assertTrue(Matrix.infNorm(w.mult(B).mult(v).sub(nbar)) < 1e-10, "W*B*V != NBAR at order " + n);
            assertTrue(Matrix.infNorm(w.mult(w.transpose()).sub(Matrix.eye(n))) < 1e-10, "W not orthogonal");
            assertTrue(Matrix.infNorm(v.mult(v.transpose()).sub(Matrix.eye(n))) < 1e-10, "V not orthogonal");

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i > j + 1) {
                        assertEquals(0.0, lbar.get(i, j), 1e-12, "LBAR not upper Hessenberg");
                    }
                    if (i > j) {
                        assertEquals(0.0, nbar.get(i, j), 1e-12, "NBAR not upper triangular");
                    }
                }
            }
        }
    }
}
