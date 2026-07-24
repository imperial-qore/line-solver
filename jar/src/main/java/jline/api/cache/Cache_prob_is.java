/**
 * @file Cache Probability Computation via Importance Sampling
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import java.util.Random;

import org.apache.commons.math3.util.FastMath;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class Cache_prob_is {
    private Cache_prob_is() {}

    /**
     * Computes cache hit probabilities using Monte Carlo importance sampling.
     */
    public static Matrix cache_prob_is(Matrix gamma, Matrix m, int samples) {
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        int mt = (int) m.elementSum();

        Matrix prob = new Matrix(n, h + 1);

        if (n == 0 || mt == 0) {
            for (int i = 0; i < n; i++) {
                prob.set(i, 0, 1.0);
            }
            return prob;
        }

        if (n < mt) {
            InputOutput.line_warning("cache_prob_is",
                    "Number of items (%d) less than cache capacity (%d).", n, mt);
            for (int i = 0; i < n; i++) {
                prob.set(i, 0, 1.0);
            }
            return prob;
        }

        if (n == mt) {
            return Cache_prob_erec.cache_prob_erec(gamma, m);
        }

        Matrix logGamma = new Matrix(n, h);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < h; j++) {
                logGamma.set(i, j, FastMath.log(gamma.get(i, j) + 1e-300));
            }
        }

        double logMFact = 0.0;
        for (int j = 0; j < h; j++) {
            logMFact += Cache_is.factln((int) m.get(j));
        }

        double logCombinations = Cache_is.factln(n) - Cache_is.factln(mt) - Cache_is.factln(n - mt);
        double logMultinomial = Cache_is.factln(mt) - logMFact;

        Matrix itemLevelWeight = new Matrix(n, h);
        double totalWeight = 0.0;

        Random random = new Random();

        for (int s = 0; s < samples; s++) {
            int[] selected = Cache_is.sampleWithoutReplacement(n, mt, random);

            int[][] assignment = Cache_is.assignItemsToLevels(m, selected, random);

            double logStateProb = logMFact;
            for (int j = 0; j < h; j++) {
                for (int item : assignment[j]) {
                    logStateProb += logGamma.get(item, j);
                }
            }

            double logProposal = -logCombinations - logMultinomial;

            double logIsWeight = logStateProb - logProposal;
            double isWeight = FastMath.exp(logIsWeight - 50);

            totalWeight += isWeight;
            for (int j = 0; j < h; j++) {
                for (int item : assignment[j]) {
                    itemLevelWeight.set(item, j, itemLevelWeight.get(item, j) + isWeight);
                }
            }
        }

        if (totalWeight > 0) {
            for (int i = 0; i < n; i++) {
                double hitSum = 0.0;
                for (int j = 0; j < h; j++) {
                    double hitProb = itemLevelWeight.get(i, j) / totalWeight;
                    prob.set(i, j + 1, hitProb);
                    hitSum += hitProb;
                }
                prob.set(i, 0, FastMath.max(0.0, 1.0 - hitSum));
            }
        } else {
            for (int i = 0; i < n; i++) {
                prob.set(i, 0, 1.0);
            }
        }

        return prob;
    }

    public static Matrix cache_prob_is(Matrix gamma, Matrix m) {
        return cache_prob_is(gamma, m, 100000);
    }
}
