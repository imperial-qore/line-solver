/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mam;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * MAP/M/1-PS sojourn law (Masuyama and Takine, ORL 31(6), 2003, Theorem 1)
 * against MATLAB map_m1ps_sojourn.m, the reference.
 */
public class MapM1psSojournTest {

    private static Matrix mat(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                m.set(i, j, a[i][j]);
            }
        }
        return m;
    }

    private static final double[] X = {0.0, 0.5, 1.0, 2.0, 5.0, 10.0};

    @Test
    public void poissonArrivalsAgreeWithMatlab() {
        Matrix C = mat(new double[][]{{-0.8}});
        Matrix D = mat(new double[][]{{0.8}});
        Map_m1ps.SojournResult r = Map_m1ps.map_m1ps_sojourn(C, D, 1.0, X);
        double[] want = {0.99999999999104117, 0.83184081113551656, 0.71021946527499691,
                0.54353558515986722, 0.29479354171635369, 0.1381707689745244};
        for (int i = 0; i < X.length; i++) {
            assertEquals(want[i], r.Wbar.get(0, i), 1e-11, "Wbar at x=" + X[i]);
        }
    }

    @Test
    public void twoPhaseMapAgreesWithMatlab() {
        Matrix C = mat(new double[][]{{-1.5, 0.5}, {0.2, -0.9}});
        Matrix D = mat(new double[][]{{1.0, 0.0}, {0.0, 0.7}});
        Map_m1ps.SojournResult r = Map_m1ps.map_m1ps_sojourn(C, D, 1.5, X);
        double[] want = {0.99999999999419409, 0.63277184098269046, 0.42815112061027738,
                0.22175179952352916, 0.049008925370148254, 0.0072473936251071592};
        for (int i = 0; i < X.length; i++) {
            assertEquals(want[i], r.Wbar.get(0, i), 1e-11, "Wbar at x=" + X[i]);
        }
        double[] wn0 = {1.0, 0.50557478779816334, 0.28556847977845812, 0.11332068867788041,
                0.016058551427373502, 0.0016558319732668102};
        double[] wn1 = {1.0, 0.67052761049728682, 0.44694145119572903, 0.2114067657857209,
                0.03568083657759797, 0.003986736335512441};
        for (int i = 0; i < X.length; i++) {
            assertEquals(wn0[i], r.WbarN.get(0, i), 1e-11, "WbarN[0] at x=" + X[i]);
            assertEquals(wn1[i], r.WbarN.get(1, i), 1e-11, "WbarN[1] at x=" + X[i]);
        }
    }

    @Test
    public void rateMatrixSolvesItsQuadratic() {
        Matrix C = mat(new double[][]{{-1.5, 0.5}, {0.2, -0.9}});
        Matrix D = mat(new double[][]{{1.0, 0.0}, {0.0, 0.7}});
        double mu = 1.5;
        Matrix R = Map_m1ps.map_compute_R(C, D, mu);
        assertEquals(0.45686728693978579, R.get(0, 0), 1e-10);
        assertEquals(0.20979937962698852, R.get(0, 1), 1e-10);
        assertEquals(0.049413581442580963, R.get(1, 0), 1e-10);
        assertEquals(0.41725308518412879, R.get(1, 1), 1e-10);
        // D + R(C - mu I) + mu R^2 = 0, the defining equation
        Matrix CmI = C.copy();
        for (int i = 0; i < 2; i++) {
            CmI.set(i, i, CmI.get(i, i) - mu);
        }
        Matrix res = D.add(1.0, R.mult(CmI)).add(mu, R.mult(R));
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(0.0, res.get(i, j), 1e-9, "residual(" + i + "," + j + ")");
            }
        }
    }

    @Test
    public void hRecursionAgreesWithMatlab() {
        Matrix C = mat(new double[][]{{-1.5, 0.5}, {0.2, -0.9}});
        Matrix D = mat(new double[][]{{1.0, 0.0}, {0.0, 0.7}});
        double[][][] h = Map_m1ps.map_m1ps_h_recursive(C, D, 1.5, 2, 3);
        // h_{n,0} = e
        for (int n = 0; n <= 2; n++) {
            assertEquals(1.0, h[n][0][0], 0.0);
            assertEquals(1.0, h[n][0][1], 0.0);
        }
        assertEquals(0.75, h[1][1][0], 1e-14);
        assertEquals(0.75, h[1][1][1], 1e-14);
        assertEquals(0.19027777777777777, h[0][3][0], 1e-14);
        assertEquals(0.19238888888888886, h[0][3][1], 1e-14);
    }

    @Test
    public void isAComplementaryDistribution() {
        Matrix C = mat(new double[][]{{-1.5, 0.5}, {0.2, -0.9}});
        Matrix D = mat(new double[][]{{1.0, 0.0}, {0.0, 0.7}});
        Map_m1ps.SojournResult r = Map_m1ps.map_m1ps_sojourn(C, D, 1.5, X);
        for (int i = 0; i < X.length; i++) {
            double v = r.Wbar.get(0, i);
            assertTrue(v >= 0.0 && v <= 1.0 + 1e-9, "Wbar out of [0,1]: " + v);
            if (i > 0) {
                assertTrue(v <= r.Wbar.get(0, i - 1) + 1e-12, "Wbar is not nonincreasing");
            }
        }
        assertEquals(1.0, r.Wbar.get(0, 0), 1e-10);
    }
}
