/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.mam.Mmap3k_fit;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Exact-recovery test for the order-3 marking fit.
 *
 * Order two is determined by (p, F, B); order three needs the second-order
 * backward moment as well, because the canonical D1 has one more nonzero. The
 * independent characteristic set is derived in
 * io/sage/proofs/mmap3k_marking_inverse.py.
 */
public class Mmap3kFitTest {

    private static double[] chars(MatrixCell mm, int c, int n) {
        Matrix D0 = mm.get(0);
        Matrix D1 = mm.get(1);
        Matrix Dc = mm.get(2 + c);
        Matrix A = D0.scale(-1.0).inv();
        Matrix T = A.mult(D1).transpose().sub(1.0, Matrix.eye(n));
        for (int j = 0; j < n; j++) {
            T.set(n - 1, j, 1.0);
        }
        Matrix rhs = new Matrix(n, 1);
        rhs.set(n - 1, 0, 1.0);
        Matrix pc = new Matrix(n, 1);
        Matrix.solve(T, rhs, pc);
        Matrix pie = pc.transpose();
        Matrix one = Matrix.ones(n, 1);
        double p = pie.mult(A).mult(Dc).mult(one).get(0, 0);
        double f = pie.mult(A).mult(Dc).mult(A).mult(one).get(0, 0) / p;
        double b = pie.mult(A).mult(A).mult(Dc).mult(one).get(0, 0) / p;
        double b2 = pie.mult(A).mult(A).mult(A).mult(Dc).mult(one).get(0, 0) / p;
        return new double[]{p, f, b, b2};
    }

    @Test
    public void recoversTheMarkingExactlyAtOrdersTwoAndThree() {
        Random rng = new Random(23);
        for (int n : new int[]{2, 3}) {
            for (int K : new int[]{2, 3, 4, 6}) {
                double[] h = new double[n];
                double[] r = new double[Math.max(n - 1, 1)];
                for (int i = 0; i < n; i++) {
                    h[i] = 0.4 + rng.nextDouble();
                }
                for (int i = 0; i < n - 1; i++) {
                    r[i] = 0.2 + 0.5 * rng.nextDouble();
                }
                double s = 0.2 + 0.5 * rng.nextDouble();
                Matrix D0 = Matrix.zeros(n, n);
                Matrix D1 = Matrix.zeros(n, n);
                for (int i = 0; i < n; i++) {
                    D0.set(i, i, -1 / h[i]);
                    if (i + 1 < n) {
                        D0.set(i, i + 1, r[i] / h[i]);
                    }
                }
                for (int i = 0; i < n; i++) {
                    D1.set(i, 0, (1 - (i + 1 < n ? r[i] : s)) / h[i]);
                }
                D1.set(n - 1, n - 1, s / h[n - 1]);

                int z = 0;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (D1.get(i, j) != 0) {
                            z++;
                        }
                    }
                }
                double[][] q = new double[z][K];
                for (int j = 0; j < z; j++) {
                    double sum = 0;
                    for (int c = 0; c < K; c++) {
                        q[j][c] = rng.nextDouble();
                        sum += q[j][c];
                    }
                    for (int c = 0; c < K; c++) {
                        q[j][c] /= sum;
                    }
                }
                MatrixCell src = new MatrixCell(2 + K);
                src.set(0, D0);
                src.set(1, D1);
                for (int c = 0; c < K; c++) {
                    Matrix Dc = Matrix.zeros(n, n);
                    int jj = 0;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            if (D1.get(i, j) != 0) {
                                Dc.set(i, j, D1.get(i, j) * q[jj][c]);
                                jj++;
                            }
                        }
                    }
                    src.set(2 + c, Dc);
                }

                double[] P = new double[K];
                double[] F = new double[K];
                double[] B = new double[K];
                double[] B2 = new double[K];
                for (int c = 0; c < K; c++) {
                    double[] ch = chars(src, c, n);
                    P[c] = ch[0];
                    F[c] = ch[1];
                    B[c] = ch[2];
                    B2[c] = ch[3];
                }

                Mmap3k_fit.Result res = Mmap3k_fit.mmap3k_fit(D0, D1, P, F, B, B2);
                String tag = "n=" + n + " K=" + K;
                assertTrue(res.exact, tag + ": a marked source must be recovered exactly");
                for (int c = 0; c < K; c++) {
                    double[] g = chars(res.mmap, c, n);
                    assertEquals(P[c], g[0], Math.abs(P[c]) * 1e-9, tag + " p_" + c);
                    assertEquals(F[c], g[1], Math.abs(F[c]) * 1e-9, tag + " F_" + c);
                    assertEquals(B[c], g[2], Math.abs(B[c]) * 1e-9, tag + " B_" + c);
                    assertEquals(B2[c], g[3], Math.abs(B2[c]) * 1e-9, tag + " B2_" + c);
                }
            }
        }
    }
}
