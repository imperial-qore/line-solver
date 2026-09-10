/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.lossn.Lossn_erlangfp;
import jline.api.lossn.Lossn_manjunath;
import jline.api.lossn.Lossn_rec;
import jline.api.lossn.Lossn_rec.LossnRecResult;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for Lossn_rec: the exact normalising constant of a loss
 * network by MDD-rec.
 *
 * <p>WHY THIS METHOD EXISTS. The admissible set of a Kelly loss network is
 * {n &gt;= 0 : A n &lt;= C} and the stationary law is independent Poisson counts
 * truncated to it, so the normalising constant is a sum of a product form over a
 * set that is finite and bounded per coordinate -- exactly what a decision
 * diagram holds and Mdd_rec walks. The route it complements is the
 * Manjunath-Sikdar residue transform, which is equally exact but whose residue
 * argument COUNTS WHOLE UNITS and so needs an integral A and C. On a fractional
 * region the analyzer used to fall back to the Erlang fixed point, an
 * approximation; these tests pin the size of the error that fallback was making
 * and show that MDD-rec removes it.</p>
 *
 * <p>THREE ORACLES: a brute-force sum over the admissible set that shares no
 * code path with the walk; Lossn_manjunath, which must agree on an INTEGRAL
 * region; and the other codebases, pinned at 12 decimals.</p>
 */
public class LossnRecApiTest {

    /** 2 links, 3 routes: routes 0 and 1 take one link each, route 2 takes both. */
    private static final double[][] A_INT = {{1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}};
    private static final double[] C_INT = {6.0, 5.0};
    private static final double[] NU_INT = {2.5, 1.8, 1.2};

    /**
     * One link with FRACTIONAL class sizes and capacity: the residue transform
     * cannot count these, and erlangfp is what the default used to fall back to.
     */
    private static final double[][] A_FRAC = {{1.5, 0.75, 2.25}};
    private static final double[] C_FRAC = {7.5};
    private static final double[] NU_FRAC = {2.0, 3.0, 1.0};

    /** MATLAB, native python and the C++ port, at %.12f. */
    private static final double[] INT_QLEN = {2.282838019229, 1.619031944010, 0.994469884949};
    private static final double[] INT_LOSS = {0.086864792309, 0.100537808883, 0.171275095875};
    private static final double INT_LG = 5.342866909910;
    private static final double[] FRAC_QLEN = {1.383844183191, 2.545477546770, 0.537948865597};
    private static final double[] FRAC_LOSS = {0.308077908404, 0.151507484410, 0.462051134403};
    private static final double FRAC_LG = 5.451139581448;

    private static Matrix mat(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                m.set(i, j, a[i][j]);
            }
        }
        return m;
    }

    /** Whether n, plus one class-`plus` call when non-negative, is admissible. */
    private static boolean feasible(int[] n, double[][] A, double[] C, int plus) {
        for (int j = 0; j < C.length; j++) {
            double s = 0;
            for (int r = 0; r < n.length; r++) {
                s += A[j][r] * (n[r] + (r == plus ? 1 : 0));
            }
            if (s > C[j] + 1e-12) {
                return false;
            }
        }
        return true;
    }

    /** The carried load, the blocking and log G, summed state by state. */
    private static double[][] brute(double[] nu, double[][] A, double[] C, int cap) {
        int K = nu.length;
        int[] n = new int[K];
        double tot = 0;
        double[] num = new double[K];
        double[] acc = new double[K];
        int total = 1;
        for (int r = 0; r < K; r++) {
            total *= cap;
        }
        for (int code = 0; code < total; code++) {
            int c = code;
            for (int r = 0; r < K; r++) {
                n[r] = c % cap;
                c /= cap;
            }
            if (!feasible(n, A, C, -1)) {
                continue;
            }
            double w = 1;
            for (int r = 0; r < K; r++) {
                double f = 1;
                for (int k = 2; k <= n[r]; k++) {
                    f *= k;
                }
                w *= Math.pow(nu[r], n[r]) / f;
            }
            tot += w;
            for (int r = 0; r < K; r++) {
                num[r] += w * n[r];
                if (feasible(n, A, C, r)) {
                    acc[r] += w;
                }
            }
        }
        double[] qlen = new double[K];
        double[] loss = new double[K];
        for (int r = 0; r < K; r++) {
            qlen[r] = num[r] / tot;
            loss[r] = 1 - acc[r] / tot;
        }
        return new double[][]{qlen, loss, {Math.log(tot)}};
    }

    @Test
    public void matchesTheBruteForceSumOnAnIntegralRegion() {
        LossnRecResult r = Lossn_rec.lossn_rec(new Matrix(NU_INT), mat(A_INT), new Matrix(C_INT));
        double[][] b = brute(NU_INT, A_INT, C_INT, 8);
        for (int k = 0; k < NU_INT.length; k++) {
            assertEquals(b[0][k], r.QLen[k], 1e-11);
            assertEquals(b[1][k], r.Loss[k], 1e-11);
        }
        assertEquals(b[2][0], r.lG, 1e-11);
        // one walk for G, one per class for the blocking ratios
        assertEquals(NU_INT.length + 1, r.niter);
    }

    @Test
    public void agreesWithTheResidueTransformWhereBothApply() {
        LossnRecResult r = Lossn_rec.lossn_rec(new Matrix(NU_INT), mat(A_INT), new Matrix(C_INT));
        Ret.lossnManjunath m = Lossn_manjunath.lossn_manjunath(
                new Matrix(NU_INT), mat(A_INT), new Matrix(C_INT));
        for (int k = 0; k < NU_INT.length; k++) {
            assertEquals(m.qLen.get(k), r.QLen[k], 1e-10);
            assertEquals(m.lossProb.get(k), r.Loss[k], 1e-10);
        }
        assertEquals(m.lG, r.lG, 1e-10);
    }

    @Test
    public void isExactWhereTheResidueTransformCannotCount() {
        LossnRecResult r =
                Lossn_rec.lossn_rec(new Matrix(NU_FRAC), mat(A_FRAC), new Matrix(C_FRAC));
        double[][] b = brute(NU_FRAC, A_FRAC, C_FRAC, 12);
        for (int k = 0; k < NU_FRAC.length; k++) {
            assertEquals(b[0][k], r.QLen[k], 1e-11);
            assertEquals(b[1][k], r.Loss[k], 1e-11);
        }
        assertEquals(b[2][0], r.lG, 1e-11);
    }

    @Test
    public void theErlangFallbackItReplacesWasNotExact() {
        // The point of the method: this is the error the old default was making
        // on a fractional region. Assert it is REAL, so a future change silently
        // routing back to erlangfp cannot pass unnoticed.
        double[][] b = brute(NU_FRAC, A_FRAC, C_FRAC, 12);
        Ret.lossnErlangFP e = Lossn_erlangfp.lossn_erlangfp(
                new Matrix(NU_FRAC), mat(A_FRAC), new Matrix(C_FRAC));
        double worst = 0;
        for (int k = 0; k < NU_FRAC.length; k++) {
            worst = Math.max(worst, Math.abs(e.lossProb.get(k) - b[1][k]));
        }
        assertTrue(worst > 1e-3, "erlangfp error was " + worst);
    }

    @Test
    public void valuesArePinnedAcrossTheCodebases() {
        LossnRecResult r = Lossn_rec.lossn_rec(new Matrix(NU_INT), mat(A_INT), new Matrix(C_INT));
        for (int k = 0; k < 3; k++) {
            assertEquals(INT_QLEN[k], r.QLen[k], 1e-10);
            assertEquals(INT_LOSS[k], r.Loss[k], 1e-10);
        }
        assertEquals(INT_LG, r.lG, 1e-10);

        LossnRecResult f =
                Lossn_rec.lossn_rec(new Matrix(NU_FRAC), mat(A_FRAC), new Matrix(C_FRAC));
        for (int k = 0; k < 3; k++) {
            assertEquals(FRAC_QLEN[k], f.QLen[k], 1e-10);
            assertEquals(FRAC_LOSS[k], f.Loss[k], 1e-10);
        }
        assertEquals(FRAC_LG, f.lG, 1e-10);
    }

    @Test
    public void aClassConsumingNothingIsRefusedByName() {
        RuntimeException ex = assertThrows(RuntimeException.class, () -> Lossn_rec.lossn_rec(
                new Matrix(new double[]{1.0, 1.0}),
                mat(new double[][]{{1.0, 0.0}}),
                new Matrix(new double[]{3.0})));
        assertTrue(ex.getMessage().contains("consumes no resource"), ex.getMessage());
    }
}
