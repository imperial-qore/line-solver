/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.lossn.Lossn_manjunath;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the exact Manjunath-Sikdar transform for loss networks
 * (Lossn_manjunath).
 *
 * The primary oracle is an independently written enumeration of the admissible
 * set: the transform's whole point is that it never enumerates, so a
 * coefficient-domain bug cannot hide behind a shared traversal. Single-link cases
 * are additionally pinned to the Erlang-B recursion, which is a closed form the
 * enumeration and the transform have no code in common with.
 *
 * The heavy-load case (nu = C = 900) is a regression on the overflow the
 * reference itself once had: nu^n/n! peaks at exp(nu)/sqrt(2 pi nu), so forming
 * the terms before rescaling overflows to infinity, and a rescaling step guarded
 * on finiteness then declines to run, sending the infinity through to g and
 * returning NaN for every metric.
 */
public class LossnManjunathApiTest {

    private static final double TOL = 1e-10;

    private static double erlangB(double rho, int C) {
        double inv = 1.0;
        for (int k = 1; k <= C; k++) {
            inv = 1.0 + inv * k / rho;
        }
        return 1.0 / inv;
    }

    private static double logFactorial(int n) {
        double s = 0.0;
        for (int k = 2; k <= n; k++) {
            s += Math.log(k);
        }
        return s;
    }

    /**
     * Independent oracle: g(C), E[n_r] and blocking by direct enumeration of the
     * admissible box. A route in no constraint row has an untruncated Poisson
     * marginal, so it carries its full load and factors exp(nu_r) out of g(C);
     * summing a truncated box over it would be a different network.
     *
     * Returns {E[n_0..], loss_0.., lG} concatenated.
     */
    private static double[] enumerate(double[] nu, double[][] A, double[] C) {
        int R = nu.length;
        int J = C.length;
        boolean[] free = new boolean[R];
        for (int r = 0; r < R; r++) {
            free[r] = true;
            for (int j = 0; j < J; j++) {
                if (A[j][r] != 0) {
                    free[r] = false;
                }
            }
        }
        int[] N = new int[R];
        for (int r = 0; r < R; r++) {
            int cap = 0;
            boolean bounded = false;
            for (int j = 0; j < J; j++) {
                if (A[j][r] > 0) {
                    int v = (int) Math.floor(C[j] / A[j][r]);
                    cap = bounded ? Math.min(cap, v) : v;
                    bounded = true;
                }
            }
            N[r] = bounded ? cap : 0;
        }
        int total = 1;
        for (int r = 0; r < R; r++) {
            total *= (N[r] + 1);
        }
        double G = 0.0;
        double[] En = new double[R];
        double[] acc = new double[R];
        int[] n = new int[R];
        for (int i = 0; i < total; i++) {
            int rem = i;
            for (int r = 0; r < R; r++) {
                n[r] = rem % (N[r] + 1);
                rem /= (N[r] + 1);
            }
            boolean feas = true;
            for (int j = 0; j < J && feas; j++) {
                double s = 0.0;
                for (int r = 0; r < R; r++) {
                    s += A[j][r] * n[r];
                }
                if (s > C[j] + 1e-12) {
                    feas = false;
                }
            }
            if (!feas) {
                continue;
            }
            double lw = 0.0;
            for (int r = 0; r < R; r++) {
                if (n[r] > 0) {
                    lw += n[r] * Math.log(nu[r]) - logFactorial(n[r]);
                }
            }
            double w = Math.exp(lw);
            G += w;
            for (int r = 0; r < R; r++) {
                En[r] += w * n[r];
            }
            for (int r = 0; r < R; r++) {
                boolean fits = true;
                for (int j = 0; j < J && fits; j++) {
                    double s = A[j][r];
                    for (int q = 0; q < R; q++) {
                        s += A[j][q] * n[q];
                    }
                    if (s > C[j] + 1e-12) {
                        fits = false;
                    }
                }
                if (fits) {
                    acc[r] += w;
                }
            }
        }
        double[] out = new double[2 * R + 1];
        double lGfree = 0.0;
        for (int r = 0; r < R; r++) {
            if (free[r]) {
                out[r] = nu[r];
                out[R + r] = 0.0;
                lGfree += nu[r];
            } else {
                out[r] = En[r] / G;
                out[R + r] = 1.0 - acc[r] / G;
            }
        }
        out[2 * R] = Math.log(G) + lGfree;
        return out;
    }

    private static void checkAgainstEnumeration(String what, double[] nu, double[][] A, double[] C) {
        Ret.lossnManjunath r = Lossn_manjunath.lossn_manjunath(new Matrix(nu), new Matrix(A), new Matrix(C));
        double[] exp = enumerate(nu, A, C);
        int R = nu.length;
        for (int k = 0; k < R; k++) {
            assertEquals(exp[k], r.qLen.get(k), TOL, what + ": QLen[" + k + "]");
            assertEquals(exp[R + k], r.lossProb.get(k), TOL, what + ": Loss[" + k + "]");
        }
        assertEquals(exp[2 * R], r.lG, TOL, what + ": lG");
        assertEquals(1, r.niter, what + ": the transform is direct");
    }

    @Test
    public void testSingleLinkTwoClassesAgainstEnumeration() {
        checkAgainstEnumeration("single link, two classes",
                new double[]{1.5, 0.7}, new double[][]{{1.0, 2.0}}, new double[]{5.0});
    }

    @Test
    public void testGlobalPlusPerClassAgainstEnumeration() {
        // The shape a FiniteCapacityRegion produces: a global job cap plus one
        // per-class cap per class.
        checkAgainstEnumeration("global plus per class",
                new double[]{2.0, 1.0},
                new double[][]{{1.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}},
                new double[]{6.0, 4.0, 3.0});
    }

    @Test
    public void testThreeLinkChainAgainstEnumeration() {
        // Overlapping rows, so more than one link is live at once and the
        // interleaved elimination order actually does something.
        checkAgainstEnumeration("three link chain",
                new double[]{1.0, 2.0, 0.5},
                new double[][]{{1.0, 1.0, 0.0}, {0.0, 1.0, 1.0}, {1.0, 0.0, 1.0}},
                new double[]{4.0, 5.0, 3.0});
    }

    @Test
    public void testGcdReducibleRowAgainstEnumeration() {
        // (4, 8) <= 20 reduces to (1, 2) <= 5, shrinking the dimension from 21
        // coefficients to 6. The answer must not move.
        checkAgainstEnumeration("gcd reducible row",
                new double[]{3.0, 1.0}, new double[][]{{4.0, 8.0}}, new double[]{20.0});
    }

    @Test
    public void testFreeRouteAgainstEnumeration() {
        // Class 2 appears in no row: it never blocks, carries its full offered
        // load, and contributes exp(nu_2) to g(C).
        checkAgainstEnumeration("free route present",
                new double[]{1.0, 2.5}, new double[][]{{1.0, 0.0}}, new double[]{4.0});
    }

    @Test
    public void testZeroLoadRouteAgainstEnumeration() {
        // nu = 0 has no logarithm to shift by, which is why it is a separate
        // branch rather than a clamp to a tiny load.
        checkAgainstEnumeration("zero load route",
                new double[]{0.0, 2.0}, new double[][]{{1.0, 1.0}}, new double[]{3.0});
    }

    @Test
    public void testClassOverflowsAlone() {
        // A single class 2 call needs 5 units of a link with 3: it is blocked in
        // every state, including the empty one.
        double[] nu = {1.0, 1.0};
        double[][] A = {{1.0, 5.0}};
        double[] C = {3.0};
        Ret.lossnManjunath r = Lossn_manjunath.lossn_manjunath(new Matrix(nu), new Matrix(A), new Matrix(C));
        assertEquals(1.0, r.lossProb.get(1), TOL, "a class that cannot fit at all is fully blocked");
        assertEquals(0.0, r.qLen.get(1), TOL, "and carries nothing");
        checkAgainstEnumeration("class overflows alone", nu, A, C);
    }

    @Test
    public void testSingleLinkMatchesErlangB() {
        double[] loads = {1.0, 5.0, 30.0, 200.0};
        int[] caps = {1, 10, 25, 210};
        for (int i = 0; i < loads.length; i++) {
            Ret.lossnManjunath r = Lossn_manjunath.lossn_manjunath(new Matrix(new double[]{loads[i]}),
                    new Matrix(new double[][]{{1.0}}),
                    new Matrix(new double[]{(double) caps[i]}));
            double b = erlangB(loads[i], caps[i]);
            assertEquals(b, r.lossProb.get(0), 1e-11,
                    "M/M/C/C blocking must be Erlang B at nu=" + loads[i] + ", C=" + caps[i]);
            assertEquals(loads[i] * (1.0 - b), r.qLen.get(0), 1e-9,
                    "carried load is nu(1-B)");
        }
    }

    @Test
    public void testHeavyLoadDoesNotOverflow() {
        // nu = C = 900 is where the pre-fix reference returned NaN.
        Ret.lossnManjunath r = Lossn_manjunath.lossn_manjunath(new Matrix(new double[]{900.0}),
                new Matrix(new double[][]{{1.0}}), new Matrix(new double[]{900.0}));
        assertTrue(Double.isFinite(r.lossProb.get(0)), "blocking must be finite at nu=C=900");
        assertTrue(Double.isFinite(r.qLen.get(0)), "carried load must be finite at nu=C=900");
        assertTrue(Double.isFinite(r.lG), "log g(C) must be finite at nu=C=900");
        double b = erlangB(900.0, 900);
        assertEquals(b, r.lossProb.get(0), 1e-12, "and must still be Erlang B");
        assertEquals(900.0 * (1.0 - b), r.qLen.get(0), 1e-9, "carried load at nu=C=900");
    }

    @Test
    public void testFractionalInputRefusedByName() {
        // The residue argument counts whole units of capacity, so a fractional
        // entry is refused rather than rounded into a different network.
        RuntimeException eA = assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Lossn_manjunath.lossn_manjunath(new Matrix(new double[]{1.0}),
                        new Matrix(new double[][]{{1.5}}), new Matrix(new double[]{3.0}));
            }
        });
        assertTrue(eA.getMessage().contains("lossn_mci"),
                "the refusal must name the algorithm that does accept it");

        RuntimeException eC = assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Lossn_manjunath.lossn_manjunath(new Matrix(new double[]{1.0}),
                        new Matrix(new double[][]{{1.0}}), new Matrix(new double[]{2.5}));
            }
        });
        assertTrue(eC.getMessage().contains("lossn_mci"),
                "a fractional capacity is refused the same way");
    }

    @Test
    public void testLiveStateCapRefused() {
        // Peak memory is the product of (C_j+1) over the simultaneously live
        // links, so an oversized region must be refused rather than allowed to
        // exhaust the heap.
        RuntimeException e = assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            public void execute() {
                Lossn_manjunath.lossn_manjunath(new Matrix(new double[]{1.0, 1.0}),
                        new Matrix(new double[][]{{1.0, 1.0}}), new Matrix(new double[]{50.0}), 4L);
            }
        });
        assertTrue(e.getMessage().contains("lossn_mci"),
                "the refusal must point at the algorithm that is unbiased at any size");
    }
}
