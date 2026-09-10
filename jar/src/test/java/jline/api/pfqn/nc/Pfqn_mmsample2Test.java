/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.nc;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * McKenna-Mitra quadrature for the repairman normalizing constant.
 *
 * This had no test, and carried three independent defects that no smoke test
 * would have caught, since it returned a plausible finite number throughout:
 * the integrand read (Z_r + L_r)*u instead of (Z_r + L_r*u), the panel width
 * was added in the log domain rather than entering as log(du), and the panels
 * were combined with a max rather than summed. The expected values below are
 * pfqn_ca, which shares no code with this routine.
 */
public class Pfqn_mmsample2Test {

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    /**
     * The quadrature error is first order in the node spacing, so the tolerance
     * is loose by design; the defects it replaces were errors of 1.5 to 8 in lG,
     * three orders of magnitude outside this band.
     */
    private static void agreesWithConvolution(double[] L, double[] N, double[] Z) {
        double expected = Pfqn_ca.pfqn_ca(row(L), row(N), row(Z)).lG;
        Ret.pfqnNc got = Pfqn_mmsample2.pfqn_mmsample2(row(L), row(N), row(Z), 1000000);
        assertEquals(expected, got.lG, 1e-3);
    }

    @Test
    public void twoClassesAgreeWithConvolution() {
        agreesWithConvolution(new double[]{0.5, 0.3}, new double[]{3, 2}, new double[]{1, 1});
    }

    @Test
    public void twoClassesUnequalThinkTimes() {
        agreesWithConvolution(new double[]{0.6, 0.4}, new double[]{2, 1}, new double[]{1, 0.5});
    }

    @Test
    public void threeClasses() {
        agreesWithConvolution(new double[]{0.9, 0.2, 0.4}, new double[]{2, 3, 1},
                new double[]{0.5, 1, 2});
    }

    /**
     * Refining the grid must reduce the error. A max over panels does not
     * converge, so this is the assertion that separates a quadrature from the
     * largest-single-term estimate it replaced.
     */
    @Test
    public void errorShrinksWithMoreSamples() {
        Matrix L = row(0.5, 0.3), N = row(3, 2), Z = row(1, 1);
        double expected = Pfqn_ca.pfqn_ca(L, N, Z).lG;
        double coarse = Math.abs(Pfqn_mmsample2.pfqn_mmsample2(L, N, Z, 10000).lG - expected);
        double fine = Math.abs(Pfqn_mmsample2.pfqn_mmsample2(L, N, Z, 1000000).lG - expected);
        assertTrue(fine < coarse / 10.0,
                "error must fall with grid refinement: coarse=" + coarse + " fine=" + fine);
    }
}
