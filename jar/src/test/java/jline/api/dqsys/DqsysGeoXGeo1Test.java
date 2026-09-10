/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.dqsys;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Analytic validation of the discrete-time Geo^X/Geo/1 batch-arrival formulas.
 *
 * Assertions are identities, not stored numbers:
 * - the degenerate batch (beta = 1) reproduces Dqsys_geogeo1 exactly;
 * - Little's law and the in-service decomposition hold in both conventions;
 * - the generating function normalizes and its derivative at 1 reproduces the
 *   closed-form mean;
 * - the closed form agrees with a brute-force power iteration of the
 *   batch-arrival chain over a grid of (a, beta, s).
 */
public class DqsysGeoXGeo1Test {

    private static final double TOL = 1e-9;

    /** (a, beta, s) with a/beta < s, so the queue is stable. */
    private static final double[][] GRID = {
            {0.1, 0.5, 0.9},   // lambda = 0.2
            {0.2, 0.8, 0.6},   // lambda = 0.25
            {0.15, 0.4, 0.8},  // lambda = 0.375
            {0.3, 1.0, 0.7},   // degenerate batch, lambda = 0.3
            {0.05, 0.2, 0.5},  // heavy batches, lambda = 0.25
            {0.4, 0.9, 0.9}    // lambda = 0.4444
    };

    // ------------------------------------------------------------------
    // Reduction to the single-arrival case
    // ------------------------------------------------------------------

    @Test
    public void degenerateBatchReproducesGeoGeo1() {
        for (GeoGeo1Convention conv : GeoGeo1Convention.values()) {
            for (double[] pt : new double[][]{{0.2, 0.5}, {0.1, 0.4}, {0.45, 0.5}, {0.25, 1.0}}) {
                double a = pt[0];
                double s = pt[1];
                GeoGeo1Result single = Dqsys_geogeo1.dqsys_geogeo1(a, s, conv);
                GeoXGeo1Result batch = Dqsys_geoxgeo1.dqsys_geoxgeo1(a, 1.0, s, conv);
                String at = " (a=" + a + ", s=" + s + ", " + conv + ")";
                assertEquals(single.getMeanQueueLength(), batch.getMeanQueueLength(), TOL,
                        "E[N] must collapse to Geo/Geo/1" + at);
                assertEquals(single.getMeanSojournTime(), batch.getMeanSojournTime(), TOL,
                        "E[T] must collapse to Geo/Geo/1" + at);
                assertEquals(single.getMeanWaitingTime(), batch.getMeanWaitingTime(), TOL,
                        "E[W] must collapse to Geo/Geo/1" + at);
                assertEquals(single.getUtilization(), batch.getUtilization(), TOL,
                        "U must collapse to Geo/Geo/1" + at);
                assertEquals(single.getThroughput(), batch.getThroughput(), TOL,
                        "X must collapse to Geo/Geo/1" + at);
            }
        }
    }

    // ------------------------------------------------------------------
    // Operational laws
    // ------------------------------------------------------------------

    @Test
    public void littlesLawHolds() {
        for (GeoGeo1Convention conv : GeoGeo1Convention.values()) {
            for (double[] pt : GRID) {
                GeoXGeo1Result r = analyze(pt, conv);
                assertEquals(r.getMeanQueueLength(), r.getThroughput() * r.getMeanSojournTime(),
                        1e-9, "E[N] = lambda E[T], " + conv);
                assertEquals(r.getMeanWaitingQueue(), r.getThroughput() * r.getMeanWaitingTime(),
                        1e-9, "E[Nw] = lambda E[W], " + conv);
            }
        }
    }

    @Test
    public void sojournSplitsIntoWaitingPlusService() {
        for (GeoGeo1Convention conv : GeoGeo1Convention.values()) {
            for (double[] pt : GRID) {
                GeoXGeo1Result r = analyze(pt, conv);
                assertEquals(r.getMeanSojournTime(),
                        r.getMeanWaitingTime() + r.getMeanServiceTime(), 1e-9,
                        "E[T] = E[W] + E[S], " + conv);
                assertEquals(r.getMeanQueueLength() - r.getMeanWaitingQueue(),
                        r.getThroughput() * r.getMeanServiceTime(), 1e-9,
                        "in-service population = lambda E[S], " + conv);
            }
        }
    }

    @Test
    public void utilizationAndBoundaryEmptyProbAreLoadComplements() {
        for (double[] pt : GRID) {
            GeoXGeo1Result r = analyze(pt, GeoGeo1Convention.LAS_DA);
            double lambda = pt[0] / pt[1];
            assertEquals(lambda, r.getArrivalRate(), TOL, "lambda = a E[X]");
            assertEquals(lambda / pt[2], r.getUtilization(), TOL, "U = lambda/s");
            assertEquals(1.0 - r.getUtilization(), r.getBoundaryEmptyProb(), TOL,
                    "p0 = 1 - lambda/s");
        }
    }

    @Test
    public void conventionsDifferByExactlyOneSlot() {
        for (double[] pt : GRID) {
            GeoXGeo1Result las = analyze(pt, GeoGeo1Convention.LAS_DA);
            GeoXGeo1Result eas = analyze(pt, GeoGeo1Convention.EAS);
            assertEquals(1.0, las.getMeanSojournTime() - eas.getMeanSojournTime(), 1e-9,
                    "sojourn differs by one slot");
            assertEquals(las.getArrivalRate(),
                    las.getMeanQueueLength() - eas.getMeanQueueLength(), 1e-9,
                    "content differs by the departures of the slot");
            assertEquals(las.getMeanWaitingTime(), eas.getMeanWaitingTime(), 1e-9,
                    "queueing delay is epoch independent");
        }
    }

    @Test
    public void batchingIncreasesCongestionAtEqualLoad() {
        // Hold lambda = a/beta fixed and make the batches larger: the extra
        // within-batch queueing must show up as a longer sojourn.
        double s = 0.8;
        double lambda = 0.4;
        double previous = 0.0;
        for (double beta : new double[]{1.0, 0.8, 0.5, 0.25}) {
            double a = lambda * beta;
            GeoXGeo1Result r = Dqsys_geoxgeo1.dqsys_geoxgeo1(a, beta, s);
            assertEquals(lambda, r.getArrivalRate(), 1e-12, "load held fixed");
            assertTrue(r.getMeanSojournTime() > previous,
                    "larger batches must not reduce the sojourn, beta=" + beta);
            previous = r.getMeanSojournTime();
        }
    }

    // ------------------------------------------------------------------
    // Generating function
    // ------------------------------------------------------------------

    @Test
    public void pgfNormalizesAndReproducesTheMean() {
        for (double[] pt : GRID) {
            double a = pt[0];
            double beta = pt[1];
            GeoXGeo1Result r = analyze(pt, GeoGeo1Convention.LAS_DA);
            assertEquals(1.0, r.pgf(1.0, 1.0), TOL, "P(1) = 1");

            // Second-order backward difference at z=1, using P(1)=1 exactly.
            // Both the numerator and the denominator of P vanish at z=1, so
            // evaluating within rounding distance of 1 loses every significant
            // digit; stay a finite step away and take the higher-order stencil.
            double h = 1e-4;
            double p1 = r.pgf(1.0 - h, slotPgf(a, beta, 1.0 - h));
            double p2 = r.pgf(1.0 - 2.0 * h, slotPgf(a, beta, 1.0 - 2.0 * h));
            double slope = (3.0 - 4.0 * p1 + p2) / (2.0 * h);
            assertEquals(r.getMeanQueueLength(), slope, 1e-4 * Math.max(1.0, r.getMeanQueueLength()),
                    "P'(1) must equal E[N], a=" + a + " beta=" + beta);
        }
    }

    // ------------------------------------------------------------------
    // Brute-force cross-check of the batch-arrival chain
    // ------------------------------------------------------------------

    @Test
    public void closedFormMatchesBruteForceChain() {
        for (double[] pt : GRID) {
            double a = pt[0];
            double beta = pt[1];
            double s = pt[2];
            double[] pi = solveBatchChain(a, beta, s);

            double mass = 0.0;
            double mean = 0.0;
            for (int n = 0; n < pi.length; n++) {
                mass += pi[n];
                mean += n * pi[n];
            }
            GeoXGeo1Result r = Dqsys_geoxgeo1.dqsys_geoxgeo1(a, beta, s, GeoGeo1Convention.LAS_DA);
            String at = " (a=" + a + ", beta=" + beta + ", s=" + s + ")";
            assertEquals(1.0, mass, 1e-8, "chain mass" + at);
            assertEquals(r.getBoundaryEmptyProb(), pi[0], 1e-5, "p0" + at);
            assertEquals(r.getMeanQueueLength(), mean, 1e-4 * Math.max(1.0, mean), "E[N]" + at);
        }
    }

    // ------------------------------------------------------------------
    // Argument validation
    // ------------------------------------------------------------------

    @Test
    public void invalidParametersThrow() {
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geoxgeo1.dqsys_geoxgeo1(0.5, 0.5, 0.9), "lambda = 1.0 > s");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geoxgeo1.dqsys_geoxgeo1(0.0, 0.5, 0.9), "a must be positive");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geoxgeo1.dqsys_geoxgeo1(0.1, 1.5, 0.9), "beta must not exceed one");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geoxgeo1.dqsys_geoxgeo1(0.1, 0.5, 1.5), "s must not exceed one");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geoxgeo1.dqsys_geoxgeo1(0.1, 0.5, 0.9, null), "convention is required");
        // A batch that arrives carries at least one job.
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geoxgeo1.dqsys_geoxgeo1_moments(0.1, 0.5, 0.0, 0.9,
                        GeoGeo1Convention.LAS_DA), "E[X] < 1");
        // E[X^2] >= E[X]^2 rules this pair out.
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geoxgeo1.dqsys_geoxgeo1_moments(0.1, 3.0, 0.0, 0.9,
                        GeoGeo1Convention.LAS_DA), "E[X(X-1)] below E[X]^2-E[X]");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geoxgeo1.dqsys_geoxgeo1(0.1, 0.5, 0.9).pgf(1.5, 1.0), "z outside (0,1]");
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    private static GeoXGeo1Result analyze(double[] pt, GeoGeo1Convention conv) {
        return Dqsys_geoxgeo1.dqsys_geoxgeo1(pt[0], pt[1], pt[2], conv);
    }

    /** pgf of the number of jobs arriving in one slot, A(z) = 1-a+a X(z). */
    private static double slotPgf(double a, double beta, double z) {
        double batch = beta * z / (1.0 - (1.0 - beta) * z);
        return 1.0 - a + a * batch;
    }

    /** P(batch size = k) for the geometric batch on {1,2,...}. */
    private static double batchPmf(double beta, int k) {
        return beta * Math.pow(1.0 - beta, k - 1);
    }

    /**
     * Power-iterates the truncated slot-boundary chain
     * {@code X(t+1) = X(t) - D(t) + A(t)} with the departure resolved first.
     */
    private static double[] solveBatchChain(double a, double beta, double s) {
        final int N = 3000;
        final int MAXBATCH = 600;
        double[] batch = new double[MAXBATCH + 1];
        for (int k = 1; k <= MAXBATCH; k++) {
            batch[k] = batchPmf(beta, k);
        }
        // Fold the truncated batch tail onto the last supported size so the
        // arrival law stays a probability distribution.
        double tail = 1.0;
        for (int k = 1; k <= MAXBATCH; k++) tail -= batch[k];
        batch[MAXBATCH] += tail;

        double[] pi = new double[N + 1];
        pi[0] = 1.0;
        double[] next = new double[N + 1];
        for (int it = 0; it < 200000; it++) {
            java.util.Arrays.fill(next, 0.0);
            for (int n = 0; n <= N; n++) {
                double mass = pi[n];
                if (mass == 0.0) continue;
                // Departure first, then the batch arrival.
                double downProb = (n >= 1) ? s : 0.0;
                for (int d = 0; d <= 1; d++) {
                    double pd = (d == 1) ? downProb : 1.0 - downProb;
                    if (pd == 0.0) continue;
                    int base = n - d;
                    // No batch this slot.
                    next[Math.min(base, N)] += mass * pd * (1.0 - a);
                    for (int k = 1; k <= MAXBATCH; k++) {
                        double pk = batch[k];
                        if (pk == 0.0) continue;
                        next[Math.min(base + k, N)] += mass * pd * a * pk;
                    }
                }
            }
            double delta = 0.0;
            for (int n = 0; n <= N; n++) delta += Math.abs(next[n] - pi[n]);
            System.arraycopy(next, 0, pi, 0, N + 1);
            if (delta < 1e-14) break;
        }
        return pi;
    }
}
