/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Test;

/**
 * Tests for the simulation output-analysis package.
 *
 * <p>The expected values are the ones the MATLAB twins produce and, where an
 * external reference exists, the ones scipy produces: the Shapiro-Wilk statistic
 * and p-value agree with {@code scipy.stats.shapiro} to 5e-10 and 1.5e-7.
 */
public class SimApiTest {
    private static final double TOL = 1e-9;

    /** A reproducible standard exponential stream. */
    private static double[] exponentialSample(int n, long seed) {
        Random rng = new Random(seed);
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            y[i] = -Math.log(1.0 - rng.nextDouble());
        }
        return y;
    }

    @Test
    public void testNormalAndTQuantiles() {
        assertEquals(0.9750021048517795, SimDist.normcdf(1.96), 1e-12);
        assertEquals(1.959963984540054, SimDist.norminv(0.975), 1e-10);
        assertEquals(2.2281388519649385, SimDist.tinv(0.975, 10), 1e-9);
        assertEquals(2.045229642132703, SimDist.tinv(0.975, 29), 1e-9);
        assertEquals(-2.093024054408309, SimDist.tinv(0.025, 19), 1e-9);
    }

    @Test
    public void testShapiroWilkAgainstScipy() {
        double[] x12 = {2.1, 3.4, 1.9, 5.6, 3.3, 4.8, 2.2, 6.1, 3.9, 4.4, 2.8, 5.1};
        HypothesisTestResult r12 = Sim_shapirowilk.sim_shapirowilk(x12);
        assertEquals(0.9500369248, r12.statistic, 1e-9);
        assertEquals(0.6375246422, r12.pvalue, 1e-8);
        assertFalse(r12.reject);

        // n = 3 uses the exact null law rather than Royston's transform
        double[] x3 = {1, 2, 5};
        HypothesisTestResult r3 = Sim_shapirowilk.sim_shapirowilk(x3);
        assertEquals(0.9230769231, r3.statistic, 1e-9);
        assertEquals(0.4632628749, r3.pvalue, 1e-8);
        assertTrue(Double.isNaN(r3.zscore));

        // 4 <= n <= 11 uses the middle branch of the transform
        double[] x7 = {1, 2, 3, 4, 5, 6, 20};
        HypothesisTestResult r7 = Sim_shapirowilk.sim_shapirowilk(x7);
        assertEquals(0.7086833563, r7.statistic, 1e-9);
        assertEquals(0.0046079816, r7.pvalue, 1e-9);
        assertTrue(r7.reject);
    }

    @Test
    public void testShapiroWilkWeightsAreAntisymmetric() {
        for (int n : new int[] {3, 4, 5, 6, 12, 32, 101}) {
            double[] a = Sim_shapirowilk.weights(n);
            for (int i = 0; i < n; i++) {
                assertEquals(-a[n - 1 - i], a[i], 1e-12, "n=" + n + " i=" + i);
            }
            if (n % 2 == 1) {
                assertEquals(0.0, a[n / 2], 1e-12, "odd n has a zero middle weight");
            }
        }
    }

    @Test
    public void testVonNeumann() {
        double[] x = {1, 5, 2, 8, 3, 9, 4, 7, 2, 6, 3, 8};
        HypothesisTestResult r = Sim_vonneumann.sim_vonneumann(x);
        assertEquals(2.8285714286, r.statistic, 1e-9);
        assertEquals(1.5666355475, r.zscore, 1e-8);
        assertEquals(0.1171999039, r.pvalue, 1e-8);
        assertFalse(r.reject);
    }

    @Test
    public void testVonNeumannRejectsARandomWalk() {
        Random rng = new Random(7);
        double[] walk = new double[60];
        double acc = 0.0;
        for (int i = 0; i < walk.length; i++) {
            acc += rng.nextGaussian();
            walk[i] = acc;
        }
        HypothesisTestResult r = Sim_vonneumann.sim_vonneumann(walk);
        assertTrue(r.reject, "a random walk is not a random sequence");
        assertTrue(r.statistic < 2.0, "positive correlation pushes the ratio below 2");
    }

    @Test
    public void testStsQuantileAreasRecoverTheVarianceParameter() {
        // for i.i.d. Exp(1), sigma_p^2 = p(1-p)/f(y_p)^2 = p/(1-p) exactly
        int b = 32;
        int m = 4000;
        double p = 0.9;
        double truth = p / (1.0 - p);
        double sumAp = 0.0;
        double sumNp = 0.0;
        int reps = 40;
        for (int r = 0; r < reps; r++) {
            double[] y = exponentialSample(b * m, 1000 + r);
            StsQuantileStats s = Sim_sts_quantile_areas.sim_sts_quantile_areas(y, b, m, p);
            sumAp += s.Ap;
            sumNp += s.Np;
            assertEquals(b, s.areas.length);
            assertEquals(b, s.bqe.length);
            assertEquals(b * m, s.n);
        }
        assertEquals(truth, sumAp / reps, 0.15 * truth, "A_p estimates sigma_p^2");
        assertEquals(truth, sumNp / reps, 0.15 * truth, "N_p estimates sigma_p^2");
    }

    @Test
    public void testStsQuantileAreasSingleBatchHasNoBetweenBatchDf() {
        double[] y = exponentialSample(2000, 5);
        StsQuantileStats s = Sim_sts_quantile_areas.sim_sts_quantile_areas(y, 1, 2000, 0.5);
        assertEquals(1, s.areas.length);
        assertTrue(Double.isFinite(s.Ap), "A_p is still defined for one batch");
        assertTrue(Double.isNaN(s.Np), "N_p needs two batches");
        assertTrue(Double.isNaN(s.Vp), "V_p needs two batches");
    }

    @Test
    public void testStsQuantileAreasFullSampleQuantileIsAnOrderStatistic() {
        double[] y = exponentialSample(6000, 31);
        double p = 0.75;
        StsQuantileStats s = Sim_sts_quantile_areas.sim_sts_quantile_areas(y, 6, 1000, p);
        double[] sorted = y.clone();
        java.util.Arrays.sort(sorted);
        assertEquals(sorted[(int) Math.ceil(6000 * p) - 1], s.quantile, 0.0,
                "the full-sample estimate is the ceil(np)-th order statistic");
    }

    @Test
    public void testFquestBracketsTheExactQuantileOnIidData() {
        // y_p = -log(1-p) for Exp(1)
        for (double p : new double[] {0.5, 0.9, 0.99}) {
            double truth = -Math.log(1.0 - p);
            int covered = 0;
            int reps = 12;
            for (int r = 0; r < reps; r++) {
                double[] y = exponentialSample(200000, 9000 + r);
                QuantileCIResult res = Sim_fquest.sim_fquest(y, p, 0.05);
                assertEquals(truth, res.estimate, 0.05 * truth,
                        "the estimate is close to the truth");
                assertTrue(res.upper > res.lower, "the interval is nondegenerate");
                assertTrue(res.truncated > 0, "the first batch is deleted");
                assertTrue(res.n > 0, "some observations are used");
                if (res.lower <= truth && truth <= res.upper) {
                    covered++;
                }
            }
            assertTrue(covered >= reps - 3, "p=" + p + " covered " + covered + "/" + reps);
        }
    }

    @Test
    public void testFquestDeliveredHalfWidthIsTheCombinedEstimator() {
        double[] y = exponentialSample(60000, 424242);
        QuantileCIResult res = Sim_fquest.sim_fquest(y, 0.9, 0.05);
        assertTrue(res.n <= y.length);
        assertEquals(res.halfwidth, (res.upper - res.lower) / 2.0, TOL);
        assertTrue(Double.isFinite(res.Vp));
        if (!res.heuristic) {
            double expected = SimDist.tinv(0.975, 2.0 * res.b - 1)
                    * Math.sqrt(res.Vp / res.n);
            assertEquals(expected, res.halfwidth, 1e-12,
                    "the delivered half-width is t*sqrt(Vp/n)");
            double vp = (res.b * res.Ap + (res.b - 1) * res.Np) / (2.0 * res.b - 1);
            assertEquals(vp, res.Vp, 1e-12, "Vp combines Ap and Np on 2b-1 degrees of freedom");
        }
    }

    @Test
    public void testFirquestBracketsTheExactQuantileOnIidData() {
        double p = 0.9;
        double truth = -Math.log(1.0 - p);
        int R = 5;
        int covered = 0;
        int reps = 10;
        for (int rep = 0; rep < reps; rep++) {
            double[][] y = new double[R][];
            for (int r = 0; r < R; r++) {
                y[r] = exponentialSample(40000, 77000 + rep * R + r);
            }
            QuantileCIResult res = Sim_firquest.sim_firquest(y, p, 0.05);
            assertEquals(R, res.R);
            assertEquals(R * res.b * res.m, res.n, "the pooled sample is R*b*m");
            assertTrue(res.upper > res.lower, "the interval is nondegenerate");
            if (res.lower <= truth && truth <= res.upper) {
                covered++;
            }
        }
        assertTrue(covered >= reps - 3, "covered " + covered + "/" + reps);
    }

    @Test
    public void testFirquestTruncatesEveryReplication() {
        int R = 4;
        double[][] y = new double[R][];
        for (int r = 0; r < R; r++) {
            y[r] = exponentialSample(30000, 555 + r);
        }
        QuantileCIResult res = Sim_firquest.sim_firquest(y, 0.5, 0.05);
        assertTrue(res.truncated > 0, "the longest warmup batch is deleted everywhere");
        assertTrue(res.b * res.m <= 30000 - res.truncated,
                "each replication contributes b*m of its truncated length");
    }

    @Test
    public void testFirquestDefaultBatchCounts() {
        assertArrayEqualsInt(new int[] {14, 11, 8, 5}, Sim_firquest.defaultBatchCounts(2));
        assertArrayEqualsInt(new int[] {10, 8, 6, 4}, Sim_firquest.defaultBatchCounts(3));
        assertArrayEqualsInt(new int[] {6, 5, 4, 3}, Sim_firquest.defaultBatchCounts(4));
        assertArrayEqualsInt(new int[] {5, 4, 3, 2}, Sim_firquest.defaultBatchCounts(5));
        assertArrayEqualsInt(new int[] {5, 4, 3, 2}, Sim_firquest.defaultBatchCounts(9));
        assertArrayEqualsInt(new int[] {4, 3, 2, 1}, Sim_firquest.defaultBatchCounts(10));
        assertArrayEqualsInt(new int[] {4, 3, 2, 1}, Sim_firquest.defaultBatchCounts(16));
        assertArrayEqualsInt(new int[] {3, 2, 1}, Sim_firquest.defaultBatchCounts(17));
        assertArrayEqualsInt(new int[] {2, 1}, Sim_firquest.defaultBatchCounts(23));
        assertArrayEqualsInt(new int[] {1}, Sim_firquest.defaultBatchCounts(33));
        assertArrayEqualsInt(new int[] {1}, Sim_firquest.defaultBatchCounts(100));
    }

    @Test
    public void testInputValidation() {
        assertThrows(IllegalArgumentException.class,
                () -> Sim_fquest.sim_fquest(exponentialSample(1000, 1), 1.0));
        assertThrows(IllegalArgumentException.class,
                () -> Sim_firquest.sim_firquest(
                        new double[][] {exponentialSample(1000, 1)}, 0.5));
        assertThrows(IllegalArgumentException.class,
                () -> Sim_vonneumann.sim_vonneumann(new double[] {1, 1, 1, 1, 1}));
        assertThrows(IllegalArgumentException.class,
                () -> Sim_shapirowilk.sim_shapirowilk(new double[] {1, 2}));
        assertThrows(IllegalArgumentException.class,
                () -> Sim_sts_quantile_areas.sim_sts_quantile_areas(
                        new double[100], 3, 10, 0.5));
    }

    private static void assertArrayEqualsInt(int[] expected, int[] actual) {
        assertEquals(expected.length, actual.length, "length");
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], actual[i], "index " + i);
        }
    }
}
