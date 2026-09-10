/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.mam.Map_gamma2;
import jline.api.mam.Map_moment;
import jline.api.mam.Mmap2k_fit;
import jline.api.mam.Mmap_backward_moment;
import jline.api.mam.Mmap_forward_moment;
import jline.api.mam.Mmap_pc;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Exact-recovery test for the closed-form MMAP(2,K) fit.
 *
 * A random canonical MMAP(2,K) is built, its characteristics are measured, and
 * the fit must return the same process: three moments, decay rate, class
 * probabilities, forward and backward moments. The marking inverse is linear
 * and independent of K (derived in io/sage/proofs/mmap2k_marking_inverse.py), so
 * recovery is exact rather than approximate, at every K.
 */
public class Mmap2kFitTest {

    private static MatrixCell source(int K, int form, double h1, double h2,
                                     double r1, double r2, double[][] q) {
        Matrix D0 = new Matrix(new double[][]{{-1 / h1, r1 / h1}, {0, -1 / h2}});
        Matrix D1 = (form == 1)
                ? new Matrix(new double[][]{{(1 - r1) / h1, 0}, {(1 - r2) / h2, r2 / h2}})
                : new Matrix(new double[][]{{0, (1 - r1) / h1}, {(1 - r2) / h2, r2 / h2}});
        MatrixCell mm = new MatrixCell(2 + K);
        mm.set(0, D0);
        mm.set(1, D1);
        for (int c = 0; c < K; c++) {
            Matrix Dc = Matrix.zeros(2, 2);
            if (form == 1) {
                Dc.set(0, 0, D1.get(0, 0) * q[0][c]);
            } else {
                Dc.set(0, 1, D1.get(0, 1) * q[0][c]);
            }
            Dc.set(1, 0, D1.get(1, 0) * q[1][c]);
            Dc.set(1, 1, D1.get(1, 1) * q[2][c]);
            mm.set(2 + c, Dc);
        }
        return mm;
    }

    @Test
    public void recoversACanonicalMmapExactlyForEveryK() {
        Random rng = new Random(11);
        for (int K : new int[]{2, 3, 4, 6, 10}) {
            for (int form : new int[]{1, 2}) {
                double h1 = 0.4 + 0.6 * rng.nextDouble();
                double h2 = 1.2 + 0.8 * rng.nextDouble();
                double r1 = 0.25 + 0.5 * rng.nextDouble();
                double r2 = 0.2 + 0.5 * rng.nextDouble();
                double[][] q = new double[3][K];
                for (int j = 0; j < 3; j++) {
                    double s = 0;
                    for (int c = 0; c < K; c++) {
                        q[j][c] = rng.nextDouble();
                        s += q[j][c];
                    }
                    for (int c = 0; c < K; c++) {
                        q[j][c] /= s;
                    }
                }
                MatrixCell src = source(K, form, h1, h2, r1, r2, q);
                MatrixCell agg = new MatrixCell(src.get(0), src.get(1));
                double m1 = Map_moment.map_moment(agg, 1);
                double m2 = Map_moment.map_moment(agg, 2);
                double m3 = Map_moment.map_moment(agg, 3);
                double g = Map_gamma2.map_gamma2(agg)[0];
                Matrix P = Mmap_pc.mmap_pc(src);
                Matrix F = Mmap_forward_moment.mmap_forward_moment(src, Matrix.ones(1, 1));
                Matrix B = Mmap_backward_moment.mmap_backward_moment(src, Matrix.ones(1, 1));
                double[] Pv = new double[K];
                double[] Fv = new double[K];
                double[] Bv = new double[K];
                for (int c = 0; c < K; c++) {
                    Pv[c] = P.get(c);
                    Fv[c] = F.get(c);
                    Bv[c] = B.get(c);
                }

                Mmap2k_fit.Result r = Mmap2k_fit.mmap2k_fit(m1, m2, m3, g, Pv, Fv, Bv);
                String tag = "K=" + K + " form=" + form;
                assertTrue(r.exact, tag + ": the closed form must apply to a canonical source");

                MatrixCell fagg = new MatrixCell(r.mmap.get(0), r.mmap.get(1));
                assertEquals(m1, Map_moment.map_moment(fagg, 1), Math.abs(m1) * 1e-8, tag + " M1");
                assertEquals(m2, Map_moment.map_moment(fagg, 2), Math.abs(m2) * 1e-8, tag + " M2");
                assertEquals(m3, Map_moment.map_moment(fagg, 3), Math.abs(m3) * 1e-7, tag + " M3");
                assertEquals(g, Map_gamma2.map_gamma2(fagg)[0], Math.abs(g) * 1e-8, tag + " gamma");

                Matrix fP = Mmap_pc.mmap_pc(r.mmap);
                Matrix fF = Mmap_forward_moment.mmap_forward_moment(r.mmap, Matrix.ones(1, 1));
                Matrix fB = Mmap_backward_moment.mmap_backward_moment(r.mmap, Matrix.ones(1, 1));
                for (int c = 0; c < K; c++) {
                    assertEquals(Pv[c], fP.get(c), Math.abs(Pv[c]) * 1e-8, tag + " p_" + c);
                    assertEquals(Fv[c], fF.get(c), Math.abs(Fv[c]) * 1e-8, tag + " F_" + c);
                    assertEquals(Bv[c], fB.get(c), Math.abs(Bv[c]) * 1e-8, tag + " B_" + c);
                }
            }
        }
    }
}
