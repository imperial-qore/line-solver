/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.wf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.api.wf.Wf_pattern_updater.ServiceParameters;
import jline.api.wf.Wf_pattern_updater.UpdatedWorkflow;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The phase-type algebra behind the workflow collapse. Every convolution has a
 * closed-form mean its result must match, and none of those means depends on
 * the representation the formula happens to build:
 * sequence, the means add; parallel, the maximum of two independent
 * exponentials has mean 1/l1 + 1/l2 - 1/(l1+l2); loop, a geometric number of
 * repetitions has mean m/(1-p); branch, the means mix with the probabilities.
 * The mean of an (alpha, T) law is alpha (-T)^-1 e, computed here by a solve,
 * so the test never reuses the assembly it is checking.
 */
public class WfPatternUpdaterTest {

    private static final double TOL = 1e-10;

    /** An exponential of rate lambda as an (alpha, T) pair. */
    private static ServiceParameters expo(double lambda) {
        Matrix alpha = Matrix.zeros(1, 1);
        alpha.set(0, 0, 1.0);
        Matrix T = Matrix.zeros(1, 1);
        T.set(0, 0, -lambda);
        return new ServiceParameters(alpha, T);
    }

    /** Erlang-k of total mean k/rate. */
    private static ServiceParameters erlang(double rate, int k) {
        Matrix alpha = Matrix.zeros(1, k);
        alpha.set(0, 0, 1.0);
        Matrix T = Matrix.zeros(k, k);
        for (int i = 0; i < k; i++) {
            T.set(i, i, -rate);
            if (i + 1 < k) T.set(i, i + 1, rate);
        }
        return new ServiceParameters(alpha, T);
    }

    /** alpha (-T)^-1 e, by Gauss-Jordan rather than by the assembly under test. */
    private static double phMean(ServiceParameters p) {
        Matrix Tm = p.getT();
        int n = Tm.getNumRows();
        double[][] a = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) a[i][j] = -Tm.get(i, j);
            a[i][n] = 1.0;
        }
        for (int c = 0; c < n; c++) {
            int piv = c;
            for (int r = c + 1; r < n; r++) {
                if (Math.abs(a[r][c]) > Math.abs(a[piv][c])) piv = r;
            }
            double[] t = a[c];
            a[c] = a[piv];
            a[piv] = t;
            for (int r = 0; r < n; r++) {
                if (r == c) continue;
                double f = a[r][c] / a[c][c];
                for (int k = c; k <= n; k++) a[r][k] -= f * a[c][k];
            }
        }
        double mean = 0.0;
        Matrix alpha = p.getAlpha();
        int idx = 0;
        for (int i = 0; i < alpha.getNumRows(); i++) {
            for (int j = 0; j < alpha.getNumCols(); j++) {
                mean += alpha.get(i, j) * (a[idx][n] / a[idx][idx]);
                idx++;
            }
        }
        return mean;
    }

    private static List<ServiceParameters> list(ServiceParameters... ps) {
        List<ServiceParameters> out = new ArrayList<ServiceParameters>();
        for (ServiceParameters p : ps) out.add(p);
        return out;
    }

    private static List<Double> probs(double... vs) {
        List<Double> out = new ArrayList<Double>();
        for (double v : vs) out.add(v);
        return out;
    }

    private static List<Integer> ints(int... vs) {
        List<Integer> out = new ArrayList<Integer>();
        for (int v : vs) out.add(v);
        return out;
    }

    @Test
    public void convolveSequenceAddsTheMeans() {
        ServiceParameters c = Wf_pattern_updater.convolveSequence(
                list(expo(2.0), erlang(4.0, 3), expo(5.0)));
        assertEquals(5, c.getAlpha().getNumRows() * c.getAlpha().getNumCols());
        assertEquals(0.5 + 0.75 + 0.2, phMean(c), TOL);

        assertEquals(1, Wf_pattern_updater.convolveSequence(list()).getAlpha().getNumCols());
        assertEquals(1.0 / 3.0, phMean(Wf_pattern_updater.convolveSequence(list(expo(3.0)))), TOL);
    }

    @Test
    public void convolveParallelGivesTheMaximumOfTwoExponentials() {
        double l1 = 2.0;
        double l2 = 3.0;
        ServiceParameters c = Wf_pattern_updater.convolveParallel(list(expo(l1), expo(l2)));
        assertEquals(3, c.getAlpha().getNumRows() * c.getAlpha().getNumCols());
        assertEquals(1.0 / l1 + 1.0 / l2 - 1.0 / (l1 + l2), phMean(c), TOL);
    }

    @Test
    public void convolveParallelDominatesEachBranchAndIsSymmetricInTheMean() {
        double m1 = phMean(Wf_pattern_updater.convolveParallel(list(erlang(4.0, 2), expo(1.0))));
        double m2 = phMean(Wf_pattern_updater.convolveParallel(list(expo(1.0), erlang(4.0, 2))));
        assertEquals(m1, m2, TOL);
        assertTrue(m1 > 1.0, "the maximum exceeds the slower branch's own mean");
    }

    @Test
    public void convolveLoopInflatesTheMeanByTheGeometricFactor() {
        assertEquals(0.5 / (1.0 - 0.4),
                phMean(Wf_pattern_updater.convolveLoop(erlang(4.0, 2), 0.4)), TOL);
        // outside (0,1) the law is returned unchanged, as in the reference
        assertEquals(0.5, phMean(Wf_pattern_updater.convolveLoop(erlang(4.0, 2), 0.0)), TOL);
        assertEquals(0.5, phMean(Wf_pattern_updater.convolveLoop(erlang(4.0, 2), 1.0)), TOL);
    }

    @Test
    public void convolveBranchesMixesTheMeansByTheBranchProbabilities() {
        List<ServiceParameters> ps = list(expo(2.0), erlang(4.0, 4));
        ServiceParameters c = Wf_pattern_updater.convolveBranches(ps, probs(0.25, 0.75));
        assertEquals(5, c.getAlpha().getNumRows() * c.getAlpha().getNumCols());
        assertEquals(0.25 * 0.5 + 0.75 * 1.0, phMean(c), TOL);

        // unnormalized probabilities are renormalized, so the answer is the same
        assertEquals(phMean(c), phMean(Wf_pattern_updater.convolveBranches(ps, probs(1.0, 3.0))), TOL);
        // a zero total falls back to a uniform choice
        assertEquals(0.75, phMean(Wf_pattern_updater.convolveBranches(ps, probs(0.0, 0.0))), TOL);
    }

    @Test
    public void removeMatrixRowsDropsEveryListedRow() {
        Matrix m = Matrix.zeros(5, 1);
        for (int i = 0; i < 5; i++) m.set(i, 0, i);
        Matrix r = Wf_pattern_updater.removeMatrixRows(m, ints(1, 3));
        assertEquals(3, r.getNumRows());
        assertEquals(0.0, r.get(0, 0), 0.0);
        assertEquals(2.0, r.get(1, 0), 0.0);
        assertEquals(4.0, r.get(2, 0), 0.0);
    }

    @Test
    public void updatePatternsCollapsesAChainAndConvolvesItsServiceLaws() {
        // 1 -> 2 -> 3 -> 4 with 2,3,4 service nodes: the sequence 2-3-4 collapses
        // onto node 2, whose law becomes the convolution of the three.
        Matrix link = Matrix.zeros(3, 3);
        link.set(0, 0, 1); link.set(0, 1, 2); link.set(0, 2, 1.0);
        link.set(1, 0, 2); link.set(1, 1, 3); link.set(1, 2, 1.0);
        link.set(2, 0, 3); link.set(2, 1, 4); link.set(2, 2, 1.0);

        Map<Integer, ServiceParameters> params = new HashMap<Integer, ServiceParameters>();
        params.put(2, expo(2.0));
        params.put(3, expo(4.0));
        params.put(4, expo(5.0));

        UpdatedWorkflow w = Wf_pattern_updater.updatePatterns(
                link, ints(2, 3, 4), ints(), ints(), ints(), params);

        assertTrue(w.getLinkMatrix().getNumRows() < link.getNumRows(), "the chain's edges are gone");
        assertTrue(w.getServiceParameters().containsKey(2));
        assertFalse(w.getServiceParameters().containsKey(3), "node 3 is absorbed");
        assertFalse(w.getServiceParameters().containsKey(4), "node 4 is absorbed");
        assertEquals(0.95, phMean(w.getServiceParameters().get(2)), TOL);
    }

    @Test
    public void updatePatternsIsTheIdentityWhenNothingIsDetected() {
        Matrix link = Matrix.zeros(1, 3);
        link.set(0, 0, 1); link.set(0, 1, 2); link.set(0, 2, 1.0);
        Map<Integer, ServiceParameters> params = new HashMap<Integer, ServiceParameters>();
        params.put(2, expo(1.0));

        UpdatedWorkflow w = Wf_pattern_updater.updatePatterns(
                link, ints(2), ints(), ints(), ints(), params);
        assertEquals(1, w.getLinkMatrix().getNumRows());
        assertEquals(1, w.getServiceParameters().size());
        assertEquals(1.0, phMean(w.getServiceParameters().get(2)), TOL);
    }

    @Test
    public void validateSequenceAcceptsAConnectedChainAndRejectsAGap() {
        Matrix link = Matrix.zeros(2, 3);
        link.set(0, 0, 0); link.set(0, 1, 1); link.set(0, 2, 1.0);
        link.set(1, 0, 1); link.set(1, 1, 2); link.set(1, 2, 1.0);
        assertTrue(Wf_sequence_detector.validateSequence(ints(0, 1, 2), link));
        assertFalse(Wf_sequence_detector.validateSequence(ints(0, 2), link));
    }
}
