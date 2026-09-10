/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.lossn.Lossn_mci;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of Monte Carlo importance-sampling summation for loss networks
 * (Lossn_mci), after Ross and Wang (1992).
 *
 * Assertions compare the Monte Carlo estimates against exact ground truth:
 * - single-link blocking against the Erlang-B recursion;
 * - a two-link multirate network against exact box enumeration of the
 *   product-form normalization constant and class blocking.
 * The estimate must lie within its own delta-method confidence interval,
 * and the log normalization constant must match the exact value.
 */
public class LossnMciApiTest {

    private static double erlangB(double rho, int C) {
        double inv = 1.0;
        for (int k = 1; k <= C; k++) {
            inv = 1.0 + inv * k / rho;
        }
        return 1.0 / inv;
    }

    /** Exact g(C) and class blocking by box enumeration. */
    private static double[] lossnExact(double[] nu, double[][] A, double[] C) {
        int R = nu.length;
        int J = C.length;
        int[] N = new int[R];
        long total = 1;
        for (int k = 0; k < R; k++) {
            double best = Double.POSITIVE_INFINITY;
            for (int j = 0; j < J; j++) {
                if (A[j][k] > 0) {
                    best = Math.min(best, Math.floor(C[j] / A[j][k]));
                }
            }
            N[k] = (int) best;
            total *= (N[k] + 1);
        }
        double g = 0.0;
        double[] gk = new double[R];
        int[] n = new int[R];
        for (long idx = 0; idx < total; idx++) {
            long rem = idx;
            double logq = 0.0;
            for (int k = 0; k < R; k++) {
                n[k] = (int) (rem % (N[k] + 1));
                rem /= (N[k] + 1);
                logq += n[k] * Math.log(nu[k]) - logFactorial(n[k]);
            }
            boolean feas = true;
            boolean[] feasK = new boolean[R];
            for (int k = 0; k < R; k++) {
                feasK[k] = true;
            }
            for (int j = 0; j < J; j++) {
                double val = 0.0;
                for (int k = 0; k < R; k++) {
                    val += A[j][k] * n[k];
                }
                if (val > C[j]) {
                    feas = false;
                }
                for (int k = 0; k < R; k++) {
                    if (val > C[j] - A[j][k]) {
                        feasK[k] = false;
                    }
                }
            }
            double q = Math.exp(logq);
            if (feas) {
                g += q;
            }
            for (int k = 0; k < R; k++) {
                if (feasK[k]) {
                    gk[k] += q;
                }
            }
        }
        double[] out = new double[R + 1];
        out[0] = g;
        for (int k = 0; k < R; k++) {
            out[k + 1] = 1.0 - gk[k] / g;   // blocking
        }
        return out;
    }

    private static double logFactorial(int n) {
        double s = 0.0;
        for (int i = 2; i <= n; i++) {
            s += Math.log(i);
        }
        return s;
    }

    @Test
    public void testSingleLinkErlangB() {
        double rho = 8.0;
        int Ccap = 10;
        Ret.lossnMCI r = Lossn_mci.lossn_mci(
                new Matrix(new double[]{rho}),
                new Matrix(new double[][]{{1.0}}),
                new Matrix(new double[]{Ccap}),
                200000, null, 1L, 0.05);
        double exactB = erlangB(rho, Ccap);
        double lo = r.lossCI.get(0, 0);
        double hi = r.lossCI.get(0, 1);
        assertTrue(lo <= exactB && exactB <= hi,
                "Erlang-B " + exactB + " outside CI [" + lo + "," + hi + "]");
        double[] ex = lossnExact(new double[]{rho}, new double[][]{{1.0}}, new double[]{Ccap});
        assertTrue(Math.abs(r.lG - Math.log(ex[0])) < 0.05,
                "lG " + r.lG + " vs exact " + Math.log(ex[0]));
    }

    @Test
    public void testTwoLinkMultirateExact() {
        double[] nu = {3.0, 1.5};
        double[][] A = {{1.0, 1.0}, {1.0, 2.0}};
        double[] C = {4.0, 5.0};
        double[] ex = lossnExact(nu, A, C);
        Ret.lossnMCI r = Lossn_mci.lossn_mci(
                new Matrix(nu), new Matrix(A), new Matrix(C),
                300000, null, 2L, 0.05);
        assertTrue(Math.abs(r.lG - Math.log(ex[0])) < 0.05,
                "lG " + r.lG + " vs exact " + Math.log(ex[0]));
        for (int k = 0; k < 2; k++) {
            double beta = ex[k + 1];
            double lo = r.lossCI.get(k, 0);
            double hi = r.lossCI.get(k, 1);
            assertTrue(lo <= beta && beta <= hi,
                    "class " + k + " beta " + beta + " outside CI [" + lo + "," + hi + "]");
            // carried-load identity QLen = nu*(1-beta)
            double expQ = nu[k] * (1.0 - r.lossProb.get(k));
            assertTrue(Math.abs(r.qLen.get(k) - expQ) < 1e-9, "QLen identity class " + k);
        }
    }

    @Test
    public void testReproducibleSeed() {
        double[] nu = {3.0, 1.5};
        double[][] A = {{1.0, 1.0}, {1.0, 2.0}};
        double[] C = {4.0, 5.0};
        Ret.lossnMCI r1 = Lossn_mci.lossn_mci(new Matrix(nu), new Matrix(A), new Matrix(C),
                50000, null, 7L, 0.05);
        Ret.lossnMCI r2 = Lossn_mci.lossn_mci(new Matrix(nu), new Matrix(A), new Matrix(C),
                50000, null, 7L, 0.05);
        assertTrue(r1.lossProb.get(0) == r2.lossProb.get(0)
                && r1.lossProb.get(1) == r2.lossProb.get(1), "fixed seed reproducible");
    }
}
