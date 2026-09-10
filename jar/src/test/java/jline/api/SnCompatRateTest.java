/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api;

import jline.api.sn.SnCompatRate;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The activated-server rate of a compatibility structure, and the scaling
 * SolverLN carries onto the layer station.
 *
 * <p>Pool layout, the 3x2 of Dorsman and Gardner Fig. 1 reduced to two
 * operands: s1 serves op0 only, s2 serves both, s3 serves op1 only.
 */
public class SnCompatRateTest {

    private static Matrix compat() {
        Matrix c = new Matrix(3, 2);
        c.set(0, 0, 1);
        c.set(1, 0, 1);
        c.set(1, 1, 1);
        c.set(2, 1, 1);
        return c;
    }

    private static final double[] COUNTS = {1, 1, 1};
    private static final double[] RATES = {1, 1, 1};
    private static final double TOL = 1e-12;

    @Test
    public void activatedRateAtIntegerStates() {
        assertEquals(2.0, SnCompatRate.snCompatRate(compat(), COUNTS, RATES, new double[]{2, 0}), TOL);
        assertEquals(3.0, SnCompatRate.snCompatRate(compat(), COUNTS, RATES, new double[]{2, 1}), TOL);
        assertEquals(2.0, SnCompatRate.snCompatRate(compat(), COUNTS, RATES, new double[]{0, 1}), TOL);
        assertEquals(0.0, SnCompatRate.snCompatRate(compat(), COUNTS, RATES, new double[]{0, 0}), TOL);
        assertEquals(3.0, SnCompatRate.snCompatPeak(COUNTS, RATES), TOL);
    }

    /**
     * min(1, .) and the hard indicator agree at EVERY integer state.
     *
     * <p>This is the invariant that keeps the order-independent law intact: the
     * fractional relaxation exists only so a mean-value solver can see the
     * structure, and CTMC and the simulators, which evaluate at integer states
     * alone, must not be able to tell the two apart.
     */
    @Test
    public void fractionalLawMatchesTheIndicatorOnTheLattice() {
        Matrix c = compat();
        for (int n0 = 0; n0 < 4; n0++) {
            for (int n1 = 0; n1 < 4; n1++) {
                double[] n = {n0, n1};
                double hard = 0;
                for (int t = 0; t < 3; t++) {
                    for (int j = 0; j < 2; j++) {
                        if (n[j] > 0 && c.get(t, j) != 0) {
                            hard += COUNTS[t] * RATES[t];
                            break;
                        }
                    }
                }
                assertEquals(hard, SnCompatRate.snCompatRate(c, COUNTS, RATES, n), TOL);
            }
        }
    }

    @Test
    public void poolScalesBelowOneCompatibleJob() {
        // only op0 present, at half a job: pools s1 and s2 reach it, each at 0.5
        assertEquals(1.0, SnCompatRate.snCompatRate(compat(), COUNTS, RATES, new double[]{0.5, 0.0}), TOL);
        // a pool sees the COMBINED load of the operands it can reach
        Matrix both = new Matrix(1, 2);
        both.set(0, 0, 1);
        both.set(0, 1, 1);
        double[] one = {1};
        assertEquals(1.0, SnCompatRate.snCompatRate(both, one, one, new double[]{0.4, 0.7}), TOL);
        assertEquals(0.6, SnCompatRate.snCompatRate(both, one, one, new double[]{0.4, 0.2}), TOL);
    }

    /**
     * A fully-compatible pool imposes NO compatibility penalty, at every
     * population.
     *
     * <p>This is what makes the lowering safe: the compatibility graph is the
     * only thing eta expresses. The statement of it is the EXACTNESS IDENTITY
     * {@code min(N,S) * (peak/S) * eta(n) == mu(n)} -- the solver applies
     * {@code min(N,S)} servers at the average server rate and eta corrects the
     * product to the activated-server rate, so a pool whose graph takes nothing
     * away is one where eta carries the whole remainder and nothing else.
     *
     * <p>eta ITSELF IS ONE ONLY AT SATURATION, and above one below it: a pool
     * of S servers facing one compatible job clears S, not 1, because every one
     * of them works on it and the first to finish cancels the rest. That
     * speed-up is exactly what the identity above has to carry, and it is why
     * the denominator damps by {@code min(1, N/S)} rather than by
     * {@code min(1, N)} -- the latter cancelled it, and LDES, which simulates
     * mu(n) directly, disagreed with the result by that factor.
     */
    @Test
    public void fullyCompatiblePoolIsNeutral() {
        Matrix both = new Matrix(1, 2);
        both.set(0, 0, 1);
        both.set(0, 1, 1);
        double[] three = {3};
        double[] one = {1};
        double[][] states = {{2, 1}, {0.61, 0.90}, {0.2, 0.3}, {0, 1}, {5, 5}};
        double peak = SnCompatRate.snCompatPeak(three, one);
        for (int k = 0; k < states.length; k++) {
            double total = states[k][0] + states[k][1];
            double eta = SnCompatRate.snCompatScaling(both, three, one, states[k]);
            // the identity, at every population
            assertEquals(SnCompatRate.snCompatRate(both, three, one, states[k]),
                    Math.min(total, three[0]) * (peak / three[0]) * eta, TOL);
            if (total >= three[0]) {
                assertEquals(1.0, eta, TOL); // saturated: nothing left to carry
            } else {
                assertTrue(eta > 1.0, "an unsaturated full pool carries the redundancy "
                        + "speed-up, got " + eta);
            }
        }
    }

    /** A partial graph must scale BELOW the fully-compatible pool it partitions. */
    @Test
    public void compatibilityGraphScalesBelowNeutral() {
        double[] n = {0.61, 0.90};
        // the same three servers, with every operand reachable from every one
        Matrix full = new Matrix(1, 2);
        full.set(0, 0, 1);
        full.set(0, 1, 1);
        double[] three = {3};
        double[] one = {1};
        double whole = SnCompatRate.snCompatScaling(full, three, one, n);
        double eta = SnCompatRate.snCompatScaling(compat(), COUNTS, RATES, n);
        assertTrue(eta < whole, "a partial compatibility graph must lose capacity against the "
                + "pool it partitions, got " + eta + " against " + whole);
        // mu = 0.61 + min(1, 1.51) + 0.90 = 2.51, over peak * min(1, 1.51/3)
        assertEquals(2.51 / 1.51, eta, 1e-9);
        assertEquals(3.0 / 1.51, whole, 1e-9);
    }
}
