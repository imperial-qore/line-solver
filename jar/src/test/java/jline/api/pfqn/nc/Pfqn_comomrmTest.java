/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.nc;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Repairman-model CoMoM against the convolution algorithm.
 *
 * The zero-think-time classes drive a separate branch of the basis construction,
 * so the fixtures below deliberately include models with none, with some, and
 * with all classes at zero think time. That last case is the one that exercises
 * the "oner" decrement whose class index was wrong: it decremented the loop
 * counter rather than the class the loop counter selects, which agrees only
 * while the zero-think-time classes occupy the leading positions.
 */
public class Pfqn_comomrmTest {

    private static final double TOL = 1e-9;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    /** log G from CoMoM-RM must equal log G from convolution, which shares no intermediate. */
    private static void agreesWithConvolution(double[] L, double[] N, double[] Z) {
        double expected = Pfqn_ca.pfqn_ca(row(L), row(N), row(Z)).lG;
        double actual = Pfqn_comomrm.pfqn_comomrm(row(L), row(N), row(Z), 1, 1e-14).lG;
        assertEquals(expected, actual, TOL);
    }

    @Test
    public void allClassesHaveThinkTime() {
        agreesWithConvolution(new double[]{0.6, 0.4}, new double[]{2, 1}, new double[]{1, 0.5});
    }

    @Test
    public void threeClassesUnorderedThinkTime() {
        agreesWithConvolution(new double[]{0.6, 0.4, 0.5}, new double[]{2, 1, 3},
                new double[]{2, 0.5, 1});
    }

    /** One zero-think-time class among several: the mixed basis branch. */
    @Test
    public void oneClassHasZeroThinkTime() {
        agreesWithConvolution(new double[]{0.6, 0.4, 0.5}, new double[]{2, 2, 1},
                new double[]{0, 0.5, 1});
    }

    /** Two of three at zero think time, so the oner loop runs more than once. */
    @Test
    public void twoClassesHaveZeroThinkTime() {
        agreesWithConvolution(new double[]{0.6, 0.4, 0.5}, new double[]{2, 1, 2},
                new double[]{0, 0, 1});
    }

    /**
     * A zero-population class is dropped by pfqn_nc_sanitize, so the class count
     * read before the call overruns the shortened arrays afterwards.
     */
    @Test
    public void zeroPopulationClassIsDropped() {
        agreesWithConvolution(new double[]{0.6, 0.4, 0.5}, new double[]{2, 0, 1},
                new double[]{1, 0.5, 0.3});
    }

    /**
     * A think time below MATLAB's hardcoded 1e-6 cutoff but above FineTol. The
     * general basis must stay well conditioned here, otherwise the tighter
     * threshold buys accuracy on paper and loses it in the arithmetic.
     */
    @Test
    public void tinyButNonzeroThinkTime() {
        agreesWithConvolution(new double[]{0.6, 0.4}, new double[]{2, 1},
                new double[]{1e-7, 0.5});
    }

    /** Every class at zero think time: the pure closed-model branch. */
    @Test
    public void allClassesHaveZeroThinkTime() {
        agreesWithConvolution(new double[]{0.6, 0.4}, new double[]{2, 2}, new double[]{0, 0});
    }
}
