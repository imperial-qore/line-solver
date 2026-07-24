/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.cache;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the cache hit-probability APIs (jline.api.cache) against
 * probability invariants and the MATLAB cache_prob_erec reference.
 */
public class CacheApiTest {

    private static final double TOL = 1e-9;

    private static Matrix accessRates() {
        // 3 items, 2 cache lists; gamma(i,j) is the access rate of item i to
        // list j (rows = items, cols = lists)
        Matrix gamma = new Matrix(3, 2);
        double[][] g = {{2.0, 0.5}, {1.0, 1.5}, {0.8, 0.3}};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 2; j++) {
                gamma.set(i, j, g[i][j]);
            }
        }
        return gamma;
    }

    private static Matrix listSizes() {
        Matrix m = new Matrix(1, 2);
        m.set(0, 0, 1);
        m.set(0, 1, 1);
        return m;
    }

    @Test
    public void cacheProbErecMatchesMatlabReferenceAndSumsToOne() {
        Matrix prob = Cache_prob_erec.cache_prob_erec(accessRates(), listSizes());
        // Reference values from MATLAB cache_prob_erec on the same input
        double[][] expected = {
                {0.25, 0.6, 0.15},
                {1.0 / 6.0, 2.0 / 15.0, 0.7},
                {7.0 / 12.0, 4.0 / 15.0, 0.15}
        };
        for (int i = 0; i < 3; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < 3; j++) {
                assertEquals(expected[i][j], prob.get(i, j), 1e-6,
                        "cache_prob_erec[" + i + "][" + j + "]");
                assertTrue(prob.get(i, j) >= -TOL, "cache probabilities must be nonnegative");
                rowSum += prob.get(i, j);
            }
            // Per item: miss probability + sum of per-list hit probabilities = 1
            assertEquals(1.0, rowSum, 1e-9, "cache row " + i + " probabilities must sum to 1");
        }
    }

    @Test
    public void cacheErecNormalizingConstantIsPositive() {
        Matrix eMat = Cache_erec.cache_erec(accessRates(), listSizes());
        double e = eMat.get(0);
        assertTrue(e > 0 && Double.isFinite(e),
                "cache rational normalizing constant must be positive and finite");
    }
}
