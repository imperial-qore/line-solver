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
 * Analytic validation of the discrete-time Geo/Geo/1 formulas.
 *
 * Every assertion is an exact identity rather than a stored number:
 * - the stationary pmf normalizes and reproduces the closed-form mean;
 * - Little's law and the in-service decomposition hold in both conventions;
 * - the two conventions differ by exactly one slot of sojourn and agree on the
 *   waiting time;
 * - the closed forms agree with a brute-force power-iteration solve of the
 *   underlying birth-death chain over a grid of (a,s).
 */
public class DqsysGeoGeo1Test {

    private static final double TOL = 1e-9;
    private static final double NUM_TOL = 1e-7; // truncated series / power iteration

    /** Parameter grid; every pair satisfies a < s so the queue is stable. */
    private static final double[][] GRID = {
            {0.1, 0.9}, {0.2, 0.5}, {0.3, 0.4}, {0.05, 0.95},
            {0.4, 0.8}, {0.45, 0.5}, {0.25, 1.0}, {0.6, 0.75}
    };

    // ------------------------------------------------------------------
    // Distributional identities
    // ------------------------------------------------------------------

    @Test
    public void pmfNormalizesInBothConventions() {
        for (GeoGeo1Convention conv : GeoGeo1Convention.values()) {
            for (double[] pt : GRID) {
                GeoGeo1Result r = Dqsys_geogeo1.dqsys_geogeo1(pt[0], pt[1], conv);
                double mass = 0.0;
                for (int n = 0; n <= 20000; n++) {
                    mass += r.pmf(n);
                }
                assertEquals(1.0, mass, NUM_TOL,
                        "pmf must normalize, " + conv + " a=" + pt[0] + " s=" + pt[1]);
            }
        }
    }

    @Test
    public void pmfMeanMatchesClosedForm() {
        for (GeoGeo1Convention conv : GeoGeo1Convention.values()) {
            for (double[] pt : GRID) {
                GeoGeo1Result r = Dqsys_geogeo1.dqsys_geogeo1(pt[0], pt[1], conv);
                double mean = 0.0;
                for (int n = 1; n <= 20000; n++) {
                    mean += n * r.pmf(n);
                }
                assertEquals(r.getMeanQueueLength(), mean, NUM_TOL,
                        "sum n*pi_n must equal E[N], " + conv + " a=" + pt[0] + " s=" + pt[1]);
            }
        }
    }

    @Test
    public void emptyProbMatchesPmfAtZero() {
        for (GeoGeo1Convention conv : GeoGeo1Convention.values()) {
            for (double[] pt : GRID) {
                GeoGeo1Result r = Dqsys_geogeo1.dqsys_geogeo1(pt[0], pt[1], conv);
                assertEquals(r.getEmptyProb(), r.pmf(0), TOL, "pi_0 accessor vs pmf(0)");
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
                GeoGeo1Result r = analyze(pt, conv);
                assertEquals(r.getMeanQueueLength(), r.getThroughput() * r.getMeanSojournTime(),
                        TOL, "E[N] = a E[T], " + conv);
                assertEquals(r.getMeanWaitingQueue(), r.getThroughput() * r.getMeanWaitingTime(),
                        TOL, "E[Nw] = a E[W], " + conv);
            }
        }
    }

    @Test
    public void sojournSplitsIntoWaitingPlusService() {
        for (GeoGeo1Convention conv : GeoGeo1Convention.values()) {
            for (double[] pt : GRID) {
                GeoGeo1Result r = analyze(pt, conv);
                assertEquals(r.getMeanSojournTime(),
                        r.getMeanWaitingTime() + r.getMeanServiceTime(), TOL,
                        "E[T] = E[W] + E[S], " + conv);
                assertEquals(r.getMeanQueueLength() - r.getMeanWaitingQueue(),
                        r.getThroughput() * r.getMeanServiceTime(), TOL,
                        "in-service population = a E[S], " + conv);
            }
        }
    }

    @Test
    public void utilizationIsLoadInBothConventions() {
        for (GeoGeo1Convention conv : GeoGeo1Convention.values()) {
            for (double[] pt : GRID) {
                GeoGeo1Result r = analyze(pt, conv);
                assertEquals(pt[0] / pt[1], r.getUtilization(), TOL, "U = a/s, " + conv);
            }
        }
    }

    // ------------------------------------------------------------------
    // Relation between the two conventions
    // ------------------------------------------------------------------

    @Test
    public void conventionsDifferByExactlyOneSlot() {
        for (double[] pt : GRID) {
            GeoGeo1Result las = Dqsys_geogeo1.dqsys_geogeo1(pt[0], pt[1], GeoGeo1Convention.LAS_DA);
            GeoGeo1Result eas = Dqsys_geogeo1.dqsys_geogeo1(pt[0], pt[1], GeoGeo1Convention.EAS);
            assertEquals(1.0, las.getMeanSojournTime() - eas.getMeanSojournTime(), TOL,
                    "LAS-DA sojourn exceeds EAS sojourn by one slot");
            assertEquals(pt[0], las.getMeanQueueLength() - eas.getMeanQueueLength(), TOL,
                    "content differs by the arrivals of the current slot");
            assertEquals(las.getMeanWaitingTime(), eas.getMeanWaitingTime(), TOL,
                    "queueing delay is convention independent");
        }
    }

    @Test
    public void lasDaIsTheDefaultConvention() {
        GeoGeo1Result def = Dqsys_geogeo1.dqsys_geogeo1(0.2, 0.5);
        GeoGeo1Result las = Dqsys_geogeo1.dqsys_geogeo1(0.2, 0.5, GeoGeo1Convention.LAS_DA);
        assertEquals(GeoGeo1Convention.LAS_DA, def.getConvention());
        assertEquals(las.getMeanSojournTime(), def.getMeanSojournTime(), TOL);
    }

    // ------------------------------------------------------------------
    // Degenerate limits
    // ------------------------------------------------------------------

    @Test
    public void vanishingLoadLeavesOnlyTheServiceTime() {
        double s = 0.4;
        GeoGeo1Result las = Dqsys_geogeo1.dqsys_geogeo1(1e-9, s, GeoGeo1Convention.LAS_DA);
        GeoGeo1Result eas = Dqsys_geogeo1.dqsys_geogeo1(1e-9, s, GeoGeo1Convention.EAS);
        assertEquals(1.0 / s, las.getMeanSojournTime(), 1e-6, "a -> 0 gives E[T] -> 1/s");
        assertEquals((1.0 - s) / s, eas.getMeanSojournTime(), 1e-6, "a -> 0 gives E[T] -> (1-s)/s");
        assertEquals(0.0, las.getMeanWaitingTime(), 1e-6, "no queueing at vanishing load");
    }

    @Test
    public void deterministicUnitServiceNeverQueues() {
        // s = 1 makes the service exactly one slot, so under LAS-DA every job
        // departs one slot after it arrives, whatever the arrival rate.
        for (double a : new double[]{0.05, 0.3, 0.7, 0.99}) {
            GeoGeo1Result r = Dqsys_geogeo1.dqsys_geogeo1(a, 1.0, GeoGeo1Convention.LAS_DA);
            assertEquals(1.0, r.getMeanSojournTime(), TOL, "a=" + a);
            assertEquals(0.0, r.getMeanWaitingTime(), TOL, "a=" + a);
            assertEquals(0.0, r.getRatio(), TOL, "tail vanishes when s=1");
        }
    }

    // ------------------------------------------------------------------
    // Brute-force cross-check of the underlying chains
    // ------------------------------------------------------------------

    @Test
    public void lasDaMatchesBruteForceChain() {
        for (double[] pt : GRID) {
            double a = pt[0];
            double s = pt[1];
            // Departure resolved before arrival within a slot: from n >= 1 the
            // content falls with probability s(1-a) and rises with a(1-s); the
            // empty state rises with probability a because no job can depart.
            double[] pi = solveBirthDeath(a, s * (1.0 - a), a * (1.0 - s));
            GeoGeo1Result r = Dqsys_geogeo1.dqsys_geogeo1(a, s, GeoGeo1Convention.LAS_DA);
            assertChainAgrees(pi, r, "LAS-DA a=" + a + " s=" + s);
        }
    }

    @Test
    public void easMatchesBruteForceChain() {
        for (double[] pt : GRID) {
            double a = pt[0];
            double s = pt[1];
            // A job arriving into an empty system may depart in the same slot,
            // so the empty state rises with probability a(1-s) as well.
            double[] pi = solveBirthDeath(a * (1.0 - s), s * (1.0 - a), a * (1.0 - s));
            GeoGeo1Result r = Dqsys_geogeo1.dqsys_geogeo1(a, s, GeoGeo1Convention.EAS);
            assertChainAgrees(pi, r, "EAS a=" + a + " s=" + s);
        }
    }

    // ------------------------------------------------------------------
    // Argument validation
    // ------------------------------------------------------------------

    @Test
    public void unstableAndOutOfRangeParametersThrow() {
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geogeo1.dqsys_geogeo1(0.5, 0.5), "a = s is not stable");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geogeo1.dqsys_geogeo1(0.7, 0.5), "a > s is not stable");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geogeo1.dqsys_geogeo1(0.0, 0.5), "a must be positive");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geogeo1.dqsys_geogeo1(0.2, 1.5), "s must not exceed one");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geogeo1.dqsys_geogeo1(1.2, 1.0), "a must not exceed one");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geogeo1.dqsys_geogeo1(0.2, 0.5, null), "convention is required");
        assertThrows(IllegalArgumentException.class,
                () -> Dqsys_geogeo1.dqsys_geogeo1(0.2, 0.5).pmf(-1), "queue length must be non-negative");
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    private static GeoGeo1Result analyze(double[] pt, GeoGeo1Convention conv) {
        return Dqsys_geogeo1.dqsys_geogeo1(pt[0], pt[1], conv);
    }

    /**
     * Power-iterates a truncated discrete birth-death chain to stationarity.
     *
     * @param up0  probability of moving 0 -&gt; 1
     * @param down probability of moving n -&gt; n-1 for n &gt;= 1
     * @param up   probability of moving n -&gt; n+1 for n &gt;= 1
     * @return the stationary distribution over the truncated state space
     */
    private static double[] solveBirthDeath(double up0, double down, double up) {
        final int N = 4000;
        double[] pi = new double[N + 1];
        pi[0] = 1.0;
        double[] next = new double[N + 1];
        for (int it = 0; it < 400000; it++) {
            java.util.Arrays.fill(next, 0.0);
            next[0] += pi[0] * (1.0 - up0);
            next[1] += pi[0] * up0;
            for (int n = 1; n <= N; n++) {
                double stay = 1.0 - down - (n < N ? up : 0.0);
                next[n - 1] += pi[n] * down;
                next[n] += pi[n] * stay;
                if (n < N) next[n + 1] += pi[n] * up;
            }
            double delta = 0.0;
            for (int n = 0; n <= N; n++) {
                delta += Math.abs(next[n] - pi[n]);
            }
            System.arraycopy(next, 0, pi, 0, N + 1);
            if (delta < 1e-14) break;
        }
        return pi;
    }

    private static void assertChainAgrees(double[] pi, GeoGeo1Result r, String label) {
        double mass = 0.0;
        double mean = 0.0;
        for (int n = 0; n < pi.length; n++) {
            mass += pi[n];
            mean += n * pi[n];
        }
        assertEquals(1.0, mass, 1e-9, "chain mass, " + label);
        assertEquals(r.getEmptyProb(), pi[0], 1e-6, "pi_0, " + label);
        assertEquals(r.getMeanQueueLength(), mean, 1e-5, "E[N], " + label);
        for (int n = 0; n <= 10; n++) {
            assertEquals(r.pmf(n), pi[n], 1e-6, "pi_" + n + ", " + label);
        }
        assertTrue(r.getRatio() < 1.0, "tail ratio must be sub-unit, " + label);
    }
}
