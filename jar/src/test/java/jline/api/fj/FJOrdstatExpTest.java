/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.fj;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Validation of E[X_(k)] for independent exponential branches against the closed forms that hold
 * independently of any implementation: the minimum is 1/sum(lambda), the i.i.d. case is a partial
 * harmonic sum, and the maximum is the alternating inclusion-exclusion series the fork-join fixed
 * point used before the quorum existed.
 */
public class FJOrdstatExpTest {

    private static final double TOL = 1e-12;

    /** E[X_(1)] = 1 / sum_i lambda_i, whatever the branch means are. */
    @Test
    public void testMinimumIsReciprocalOfTheRateSum() {
        double[] ri = {1.0 / 3.0, 1.0 / 4.0, 1.0 / 2.0};
        assertEquals(1.0 / 9.0, FJ_ordstat_exp.fj_ordstat_exp(ri, 1), TOL);
    }

    /** For n i.i.d. Exp(1) branches, E[X_(k)] = sum_{j=n-k+1..n} 1/j. */
    @Test
    public void testIidCaseIsAPartialHarmonicSum() {
        double[] ri = {1.0, 1.0, 1.0, 1.0};
        assertEquals(1.0 / 4.0, FJ_ordstat_exp.fj_ordstat_exp(ri, 1), TOL);
        assertEquals(1.0 / 4.0 + 1.0 / 3.0, FJ_ordstat_exp.fj_ordstat_exp(ri, 2), TOL);
        assertEquals(1.0 / 4.0 + 1.0 / 3.0 + 1.0 / 2.0, FJ_ordstat_exp.fj_ordstat_exp(ri, 3), TOL);
        assertEquals(1.0 / 4.0 + 1.0 / 3.0 + 1.0 / 2.0 + 1.0,
                FJ_ordstat_exp.fj_ordstat_exp(ri, 4), TOL);
    }

    /**
     * At k = n the sum must reproduce the classical E[max] series TERM BY TERM, since that is what
     * keeps a standard join bit-identical to the pre-quorum fork-join fixed point.
     */
    @Test
    public void testMaximumMatchesTheClassicalSeries() {
        double[] ri = {0.7, 1.3, 2.9, 0.4};
        int n = ri.length;
        double[] lambda = new double[n];
        for (int i = 0; i < n; i++) {
            lambda[i] = 1.0 / ri[i];
        }
        double expected = 0;
        for (int mask = 1; mask < (1 << n); mask++) {
            double s = 0;
            int bits = 0;
            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    s += lambda[i];
                    bits++;
                }
            }
            expected += (bits % 2 == 1 ? 1 : -1) / s;
        }
        assertEquals(expected, FJ_ordstat_exp.fj_ordstat_exp(ri, n), 1e-10);
    }

    /** The order statistic is nondecreasing in k. */
    @Test
    public void testMonotoneInK() {
        double[] ri = {0.5, 1.5, 2.5, 3.5, 4.5};
        double prev = -1;
        for (int k = 1; k <= ri.length; k++) {
            double m = FJ_ordstat_exp.fj_ordstat_exp(ri, k);
            org.junit.jupiter.api.Assertions.assertTrue(m > prev, "E[X_(k)] must increase with k");
            prev = m;
        }
    }

    /**
     * A branch of zero mean completes instantly: it counts toward the quorum at once and never
     * delays the join.
     */
    @Test
    public void testZeroBranchCompletesInstantly() {
        double[] ri = {0.0, 2.0};
        assertEquals(0.0, FJ_ordstat_exp.fj_ordstat_exp(ri, 1), TOL);
        assertEquals(2.0, FJ_ordstat_exp.fj_ordstat_exp(ri, 2), TOL);
    }

    /** A single branch is its own order statistic, and an empty branch set takes no time. */
    @Test
    public void testDegenerateBranchSets() {
        assertEquals(3.0, FJ_ordstat_exp.fj_ordstat_exp(new double[]{3.0}, 1), TOL);
        assertEquals(0.0, FJ_ordstat_exp.fj_ordstat_exp(new double[]{}, 1), TOL);
    }

    /** k outside [1, n] is an input error, not a silently clamped answer. */
    @Test
    public void testQuorumOutOfRangeThrows() {
        double[] ri = {1.0, 2.0};
        assertThrows(IllegalArgumentException.class, () -> FJ_ordstat_exp.fj_ordstat_exp(ri, 0));
        assertThrows(IllegalArgumentException.class, () -> FJ_ordstat_exp.fj_ordstat_exp(ri, 3));
    }
}
