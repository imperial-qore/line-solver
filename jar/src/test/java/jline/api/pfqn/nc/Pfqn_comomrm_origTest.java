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
 * The three routines pfqn_comomrm_orig, pfqn_comom and pfqn_comomrm all return
 * log G for a single queueing station plus a delay, so pfqn_ca is an independent
 * oracle for every one of them: it shares no intermediate with the CoMoM
 * recursion. Reference values are the MATLAB ones, which agree to all printed
 * digits with pfqn_ca in each case.
 *
 * The third fixture is the load-bearing one for the port: R = 3 drives the P
 * expansion twice and its think times are neither ascending nor descending, so
 * a class permutation that is not carried through every array shows up here.
 */
public class Pfqn_comomrm_origTest {

    private static final double TOL = 1e-9;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    private static void check(double expected, double[] L, double[] N, double[] Z) {
        assertEquals(expected, Pfqn_comomrm_orig.pfqn_comomrm_orig(row(L), row(N), row(Z)), TOL);
        assertEquals(expected, Pfqn_ca.pfqn_ca(row(L), row(N), row(Z)).lG, TOL);
    }

    @Test
    public void twoClassesAscendingThinkTime() {
        check(0.4101209196, new double[]{0.6, 0.4}, new double[]{2, 1}, new double[]{0.5, 1});
    }

    @Test
    public void twoClassesDescendingThinkTime() {
        check(0.6108519378, new double[]{0.6, 0.4}, new double[]{2, 1}, new double[]{1, 0.5});
    }

    @Test
    public void threeClassesUnorderedThinkTime() {
        check(2.3339868561, new double[]{0.6, 0.4, 0.5}, new double[]{2, 1, 3},
                new double[]{2, 0.5, 1});
    }

    /**
     * F1r and F2r are now built once per class rather than once per job. This
     * asserts the numbers are unchanged when several jobs of one class follow the
     * first, which is the only case where the removed rebuild could have differed.
     */
    @Test
    public void manyJobsPerClassAgreeWithConvolution() {
        check(Pfqn_ca.pfqn_ca(row(new double[]{0.6, 0.4}), row(new double[]{5, 4}),
                        row(new double[]{1, 0.5})).lG,
                new double[]{0.6, 0.4}, new double[]{5, 4}, new double[]{1, 0.5});
    }
}
