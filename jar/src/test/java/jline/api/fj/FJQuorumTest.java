/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.fj;

import jline.api.fj.FJ_quorum.StepCDF;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the k-of-n quorum join kernel against exact identities.
 *
 * <p>There is no external oracle for quorum joins: lqns rejects quorum input unless built with
 * GSL, and its GSL-gated quorum code no longer compiles. These tests therefore check the kernel
 * against closed forms that hold independently of any implementation.</p>
 */
public class FJQuorumTest {

    private static final double TOL = 1e-9;

    /** Deterministic branch, i.e. a single point mass at d. */
    private static StepCDF det(double d) {
        return new StepCDF(new double[]{d}, new double[]{1.0});
    }

    /** Exponential branch sampled on a fine grid, for comparison against closed forms. */
    private static StepCDF exponential(double rate, double horizon, int points) {
        double[] t = new double[points];
        double[] A = new double[points];
        for (int i = 0; i < points; i++) {
            t[i] = horizon * (i + 1) / points;
            A[i] = 1.0 - Math.exp(-rate * t[i]);
        }
        return new StepCDF(t, A);
    }

    /**
     * With k = n the join waits for every branch, so its CDF is the product of the branch CDFs.
     */
    @Test
    public void testKEqualsNIsProductOfCDFs() {
        StepCDF[] b = {det(1.0), det(2.0), det(3.0)};
        StepCDF join = FJ_quorum.quorumKofN(b, 3);
        // The maximum of 1, 2, 3 is 3, so no mass accumulates before t = 3.
        assertEquals(0.0, join.at(2.999), TOL);
        assertEquals(1.0, join.at(3.0), TOL);
        assertEquals(3.0, join.mean(), TOL);
        assertEquals(0.0, join.variance(), TOL);
    }

    /**
     * With k = 1 the join fires on the first completion, so its value is the minimum.
     */
    @Test
    public void testKEqualsOneIsMinimum() {
        StepCDF[] b = {det(1.0), det(2.0), det(3.0)};
        StepCDF join = FJ_quorum.quorumKofN(b, 1);
        assertEquals(0.0, join.at(0.999), TOL);
        assertEquals(1.0, join.at(1.0), TOL);
        assertEquals(1.0, join.mean(), TOL);
    }

    /**
     * For deterministic branches the k-th order statistic is the k-th smallest deadline.
     */
    @Test
    public void testDeterministicBranchesGiveKthSmallest() {
        StepCDF[] b = {det(5.0), det(1.0), det(3.0), det(9.0)};
        double[] expected = {1.0, 3.0, 5.0, 9.0};
        for (int k = 1; k <= 4; k++) {
            StepCDF join = FJ_quorum.quorumKofN(b, k);
            assertEquals(expected[k - 1], join.mean(), TOL, "k=" + k);
            assertEquals(0.0, join.variance(), TOL, "k=" + k);
        }
    }

    /**
     * For identically distributed branches the kernel must reproduce the closed binomial form of
     * the order-statistic CDF already in FJ_order_stat, which is derived independently.
     */
    @Test
    public void testAgreesWithIidOrderStatisticCDF() {
        int n = 4;
        double rate = 0.7;
        double horizon = 20.0;
        StepCDF[] b = new StepCDF[n];
        for (int i = 0; i < n; i++) {
            b[i] = exponential(rate, horizon, 4000);
        }
        for (int k = 1; k <= n; k++) {
            StepCDF join = FJ_quorum.quorumKofN(b, k);
            for (double x : new double[]{0.5, 1.0, 2.0, 4.0, 8.0}) {
                double marginal = 1.0 - Math.exp(-rate * x);
                double expected = FJ_order_stat.fj_order_stat_cdf(marginal, k, n);
                assertEquals(expected, join.at(x), 1e-3, "k=" + k + " t=" + x);
            }
        }
    }

    /**
     * For n i.i.d. exponentials of rate mu the k-th order statistic has mean
     * sum_{i=n-k+1..n} 1/(i*mu), by the memorylessness of the successive gaps.
     */
    @Test
    public void testIidExponentialOrderStatisticMean() {
        int n = 3;
        double rate = 2.0;
        StepCDF[] b = new StepCDF[n];
        for (int i = 0; i < n; i++) {
            b[i] = exponential(rate, 40.0, 40000);
        }
        for (int k = 1; k <= n; k++) {
            double expected = 0.0;
            for (int i = n - k + 1; i <= n; i++) {
                expected += 1.0 / (i * rate);
            }
            StepCDF join = FJ_quorum.quorumKofN(b, k);
            assertEquals(expected, join.mean(), 2e-3, "k=" + k);
        }
    }

    /**
     * The completion time must be monotone in the quorum: requiring more branches cannot make the
     * join finish sooner.
     */
    @Test
    public void testMeanIsMonotoneInQuorum() {
        double[] means = {1.0, 2.0, 3.0, 4.0, 5.0};
        double[] vars = {0.5, 1.0, 2.0, 1.5, 0.25};
        double prev = -1.0;
        for (int k = 1; k <= means.length; k++) {
            double m = FJ_quorum.quorumMoments(means, vars, k)[0];
            assertTrue(m >= prev - TOL, "k=" + k + " mean " + m + " < previous " + prev);
            prev = m;
        }
    }

    /**
     * A single branch reduces to that branch regardless of the quorum.
     */
    @Test
    public void testSingleBranch() {
        double[] moments = FJ_quorum.quorumMoments(new double[]{3.0}, new double[]{4.0}, 1);
        assertEquals(3.0, moments[0], 1e-6);
        assertEquals(4.0, moments[1], 1e-6);
    }

    /**
     * The three-point fit matches the first two moments of the branch it represents.
     */
    @Test
    public void testThreePointFitMatchesMoments() {
        double[][] cases = {{1.0, 0.25}, {2.0, 4.0}, {0.5, 0.1}, {10.0, 100.0}};
        for (double[] c : cases) {
            StepCDF f = FJ_quorum.threePointFit(c[0], c[1]);
            assertEquals(c[0], f.mean(), 1e-9, "mean of " + c[0] + "," + c[1]);
            assertEquals(c[1], f.variance(), 1e-9, "var of " + c[0] + "," + c[1]);
        }
    }

    /**
     * A deterministic branch, i.e. zero variance, degenerates to a single point.
     */
    @Test
    public void testThreePointFitDeterministic() {
        StepCDF f = FJ_quorum.threePointFit(2.5, 0.0);
        assertEquals(1, f.t.length);
        assertEquals(2.5, f.mean(), TOL);
        assertEquals(0.0, f.variance(), TOL);
    }

    /**
     * The result must remain a valid, monotone CDF in [0, 1] for a branch count large enough
     * that the alternating inclusion-exclusion sum would lose all significant digits.
     */
    @Test
    public void testLargeBranchCountIsStable() {
        int n = 40;
        StepCDF[] b = new StepCDF[n];
        for (int i = 0; i < n; i++) {
            b[i] = FJ_quorum.threePointFit(1.0 + 0.1 * i, 0.5 + 0.05 * i);
        }
        for (int k : new int[]{1, n / 2, n}) {
            StepCDF join = FJ_quorum.quorumKofN(b, k);
            double prev = -1.0;
            for (int i = 0; i < join.t.length; i++) {
                double v = join.A[i];
                assertTrue(v >= -1e-12 && v <= 1.0 + 1e-12, "k=" + k + " CDF value out of [0,1]: " + v);
                assertTrue(v >= prev - 1e-12, "k=" + k + " CDF is not monotone at index " + i);
                prev = v;
            }
            assertTrue(join.mean() > 0.0, "k=" + k + " mean must be positive");
        }
    }

    /**
     * A branch count above the documented maximum must fail immediately rather than hang.
     */
    @Test
    public void testBranchCountGuard() {
        int n = FJ_quorum.MAX_BRANCHES + 1;
        StepCDF[] b = new StepCDF[n];
        for (int i = 0; i < n; i++) {
            b[i] = det(1.0);
        }
        try {
            FJ_quorum.quorumKofN(b, 1);
            throw new AssertionError("expected the branch count guard to reject n=" + n);
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("branches"), expected.getMessage());
        }
    }

    /**
     * A quorum join over heterogeneous branches must lie between the minimum and the maximum.
     */
    @Test
    public void testBoundedByMinAndMax() {
        double[] means = {1.0, 4.0, 9.0};
        double[] vars = {0.1, 2.0, 5.0};
        double min = FJ_quorum.quorumMoments(means, vars, 1)[0];
        double mid = FJ_quorum.quorumMoments(means, vars, 2)[0];
        double max = FJ_quorum.quorumMoments(means, vars, 3)[0];
        assertTrue(min <= mid + TOL);
        assertTrue(mid <= max + TOL);
        assertTrue(min > 0.0);
    }
}
