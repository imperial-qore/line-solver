/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.trace;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the marked-trace statistics APIs (jline.api.trace) on a
 * deterministic two-class trace with hand-computable statistics:
 *
 * T = {1,2,1,2,1,2}, L = {0,1,0,1,0,1}
 * - class 0 interarrivals: {1,1,1} -> mean 1, second moment 1
 * - class 1 interarrivals: {2,2,2} -> mean 2, second moment 4
 * - class probabilities: {1/2, 1/2}
 * - label transitions (5 pairs): 0->1 three times, 1->0 twice
 */
public class TraceApiTest {

    private static final double TOL = 1e-12;

    private static final double[] T = {1, 2, 1, 2, 1, 2};
    private static final int[] L = {0, 1, 0, 1, 0, 1};

    @Test
    public void perClassMeansMatchHandComputation() {
        Matrix mean = Mtrace_mean.mtrace_mean(T, 2, L);
        assertEquals(1.0, mean.get(0), TOL, "class 0 mean interarrival");
        assertEquals(2.0, mean.get(1), TOL, "class 1 mean interarrival");
    }

    @Test
    public void classProbabilitiesMatchHandComputation() {
        Matrix pc = Mtrace_pc.mtrace_pc(T, L);
        assertEquals(0.5, pc.get(0), TOL, "class 0 probability");
        assertEquals(0.5, pc.get(1), TOL, "class 1 probability");
    }

    @Test
    public void transitionCrossMomentsMatchHandComputation() {
        // mtrace_moment_simple returns MC(i,j) = mean of T[t]^k conditioned on
        // the label transition i->j. In the alternating trace every 0->1
        // transition carries T=2 and every 1->0 transition carries T=1; the
        // diagonal transitions never occur (NaN by definition).
        Matrix m1 = Mtrace_moment_simple.mtrace_moment_simple(T, L, 1);
        assertEquals(2.0, m1.get(0, 1), TOL, "0->1 first cross moment");
        assertEquals(1.0, m1.get(1, 0), TOL, "1->0 first cross moment");
        assertTrue(Double.isNaN(m1.get(0, 0)), "0->0 transition never occurs");
        assertTrue(Double.isNaN(m1.get(1, 1)), "1->1 transition never occurs");
        Matrix m2 = Mtrace_moment_simple.mtrace_moment_simple(T, L, 2);
        assertEquals(4.0, m2.get(0, 1), TOL, "0->1 second cross moment");
        assertEquals(1.0, m2.get(1, 0), TOL, "1->0 second cross moment");
    }

    @Test
    public void labelTransitionFrequenciesMatchHandComputation() {
        Matrix sigma = Mtrace_sigma.mtrace_sigma(T, L);
        // 5 consecutive pairs: (0,1) x3, (1,0) x2, no same-label pairs
        assertEquals(0.0, sigma.get(0, 0), TOL, "0->0 frequency");
        assertEquals(3.0 / 5.0, sigma.get(0, 1), TOL, "0->1 frequency");
        assertEquals(2.0 / 5.0, sigma.get(1, 0), TOL, "1->0 frequency");
        assertEquals(0.0, sigma.get(1, 1), TOL, "1->1 frequency");
    }

    @Test
    public void mergedTracesInterleaveDeterministically() {
        // Trace 1 events at 2,4 (interarrivals {2,2}); trace 2 events at 1,3
        double[] t1 = {2, 2};
        double[] t2 = {1, 2};
        jline.util.Pair<double[], int[]> merged = Mtrace_merge.mtrace_merge(t1, t2);
        assertNotNull(merged, "mtrace_merge returned null");
        double[] mergedIat = merged.getLeft();
        int[] labels = merged.getRight();
        assertEquals(4, labels.length, "merged trace must have 4 events");
        // Reconstruct absolute times and verify the interleaving 1,2,3,4
        double tAbs = 0;
        double[] expected = {1, 2, 3, 4};
        for (int i = 0; i < mergedIat.length; i++) {
            tAbs += mergedIat[i];
            assertEquals(expected[i], tAbs, 1e-9, "merged event time " + i);
        }
    }

    @Test
    public void summaryIsConsistentWithMoments() {
        Mtrace_summary.MtraceSummary summary = Mtrace_summary.mtrace_summary(T, L);
        assertNotNull(summary, "mtrace_summary returned null");
        // Aggregate moments: M1 = (1+2)/2 = 1.5, M2 = (1+4)/2 = 2.5
        assertEquals(1.5, summary.M[0], TOL, "aggregate first moment");
        assertEquals(2.5, summary.M[1], TOL, "aggregate second moment");
        assertEquals(0.5, summary.Pc.get(0), TOL, "summary class 0 probability");
        assertEquals(0.5, summary.Pc.get(1), TOL, "summary class 1 probability");
    }

    @Test
    public void iat2countsProducesConsistentCountProcess() {
        // With scale 3, windows of length 3 each contain exactly one class-0
        // and one class-1 arrival (period is 3 time units)
        Matrix counts = Mtrace_iat2counts.mtrace_iat2counts(T, L, 3.0);
        assertNotNull(counts, "mtrace_iat2counts returned null");
        double total = 0;
        for (int i = 0; i < counts.getNumRows(); i++) {
            for (int j = 0; j < counts.getNumCols(); j++) {
                double v = counts.get(i, j);
                assertTrue(v >= 0, "counts must be nonnegative");
                total += v;
            }
        }
        assertTrue(total > 0, "count process must register the arrivals");
    }
}
