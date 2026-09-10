/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.cache;

import java.util.List;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the cost-capped cache normalizing constant against the worked
 * example of Casale-Gast, IEEE/ACM Trans. Networking 29(2), 2021, Sec. IX
 * (Model 7: tree cache, n = 10, m = (2,1,1,2), sizes 1 and 2, caps (2,1,2,4)).
 */
public class CacheCostCapTest {

    private static final double TOL = 1e-4;

    private static Matrix paperGamma() {
        // Fig. 2 tree: p(1)=p(2)=0, p(3)=p(4)=2. Stream 1 requests items 1..5 at
        // rate 0.9 with c0 = 0.025 and c = c' = 0; stream 2 requests items 6..10
        // at rate 1.0 with c0 = 1.0 and c = c' = 0.5.
        Matrix gamma = new Matrix(10, 4);
        for (int i = 0; i < 5; i++) {
            gamma.set(i, 0, 0.9 * (1 - 0.025));
            gamma.set(i, 1, 0.9 * 0.025);
            gamma.set(i, 2, 0.0);
            gamma.set(i, 3, 0.0);
        }
        for (int i = 5; i < 10; i++) {
            gamma.set(i, 0, 0.0);
            gamma.set(i, 1, 1.0);
            gamma.set(i, 2, 0.5);
            gamma.set(i, 3, 0.5);
        }
        return gamma;
    }

    private static Matrix rowVector(double... v) {
        Matrix out = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            out.set(0, i, v[i]);
        }
        return out;
    }

    @Test
    public void constrainedNormalizingConstantMatchesThePaper() {
        Matrix gamma = paperGamma();
        Matrix m = rowVector(2, 1, 1, 2);
        Matrix sigma = rowVector(1, 1, 1, 1, 1, 2, 2, 2, 2, 2);
        Matrix k = rowVector(2, 1, 2, 4);

        assertEquals(238.7982, Cache_erec.cache_erec(gamma, m).value(), TOL);
        assertEquals(7.7963, Cache_erec.cache_erec(gamma, m, sigma, k).value(), TOL);
    }

    @Test
    public void unitSizesWithSaturatedCapsReduceToTheUnconstrainedConstant() {
        Matrix gamma = paperGamma();
        Matrix m = rowVector(2, 1, 1, 2);
        Matrix unit = rowVector(1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
        assertEquals(Cache_erec.cache_erec(gamma, m).value(),
                Cache_erec.cache_erec(gamma, m, unit, m).value(), 1e-9);
    }

    @Test
    public void constrainedMarginalsFillEveryListExactly() {
        Matrix gamma = paperGamma();
        Matrix m = rowVector(2, 1, 1, 2);
        Matrix sigma = rowVector(1, 1, 1, 1, 1, 2, 2, 2, 2, 2);
        Matrix k = rowVector(2, 1, 2, 4);
        Matrix pij = Cache_prob_erec.cache_prob_erec(gamma, m, sigma, k);
        for (int j = 0; j < 4; j++) {
            double colSum = 0.0;
            for (int i = 0; i < 10; i++) {
                colSum += pij.get(i, j + 1);
            }
            assertEquals(m.get(j), colSum, 1e-9);
        }
        // items of size 2 cannot enter list 2, whose cap is 1
        for (int i = 5; i < 10; i++) {
            assertEquals(0.0, pij.get(i, 2), 1e-12);
        }
        // the mean per-list cost saturates every cap in this instance
        Matrix K = Cache_cost.cache_cost(gamma, m, sigma, k, pij);
        for (int j = 0; j < 4; j++) {
            assertEquals(k.get(j), K.get(j), 1e-9);
        }
    }

    @Test
    public void samplingEstimatorAgreesWithTheExactConstrainedConstant() {
        // small linear cache where every size-feasible state is reachable
        Matrix gamma = new Matrix(6, 2);
        for (int i = 0; i < 6; i++) {
            gamma.set(i, 0, 1.0 / 3.0);
            gamma.set(i, 1, 1.0 / 9.0);
        }
        Matrix m = rowVector(1, 1);
        Matrix sigma = rowVector(1, 1, 1, 2, 2, 2);
        Matrix k = rowVector(2, 1);
        double exact = Cache_erec.cache_erec(gamma, m, sigma, k).value();
        double sampled = Cache_is.cache_is(gamma, m, 200000, sigma, k).E;
        assertEquals(exact, sampled, 0.05 * exact);
    }

    @Test
    public void pathCheckFlagsTheBlockedPromotionOfThePaperInstance() {
        Matrix gamma = paperGamma();
        Matrix sigma = rowVector(1, 1, 1, 1, 1, 2, 2, 2, 2, 2);
        Matrix k = rowVector(2, 1, 2, 4);
        int[] parent = {-1, -1, 1, 1}; // p(1)=p(2)=miss, p(3)=p(4)=list 2
        List<Cache_cost_pathcheck.BlockedPair> viol =
                Cache_cost_pathcheck.cache_cost_pathcheck(gamma, sigma, k, parent);
        // size-2 items may sit in lists 3 and 4 but can never traverse list 2
        assertFalse(viol.isEmpty());
        for (Cache_cost_pathcheck.BlockedPair p : viol) {
            assertTrue(p.item >= 5);
            assertTrue(p.list == 2 || p.list == 3);
            assertEquals(1, p.blockingList);
        }
        // the same cache without cost caps on the intermediate list is clean
        Matrix kOpen = rowVector(2, 2, 2, 4);
        assertTrue(Cache_cost_pathcheck.cache_cost_pathcheck(gamma, sigma, kOpen, parent).isEmpty());
    }
}
