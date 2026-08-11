/**
 * @file Cache Probability Analysis via Ray Integration
 *
 * Computes cache state probabilities using ray integration methods for
 * partial differential equation solutions. Provides PDE-based approximations
 * for complex cache systems where traditional methods become intractable.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.Collections;
import java.util.stream.IntStream;

public final class Cache_prob_spm {
    private Cache_prob_spm() {}

    /**
     * Computes the cache state probabilities using the ray method.
     * This method calculates the probabilities of the cache being in different states based on the access factors
     * and capacity, utilizing the logarithm of the partition function obtained through the ray method.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m     Matrix representing the cache capacity vector.
     * @return Matrix containing the computed cache state probabilities for each item and level.
     */
    public static Matrix cache_prob_spm(Matrix gamma, Matrix m) {
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        double lE = Cache_spm.cache_spm(gamma, m).lZ;
        Matrix prob = new Matrix(n, h + 1);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < h; j++) {
                Matrix subGamma = gamma.copy();
                subGamma.removeRows(Collections.singleton(i));
                double lEi = Cache_spm.cache_spm(subGamma, Matrix.oner(m, j)).lZ;
                double value = m.get(j) * gamma.get(i, j) * FastMath.exp(lEi - lE);
                prob.set(i, j + 1, value);
            }
            final int finalI = i;
            double rowSum = IntStream.range(1, h + 1).mapToDouble(j -> prob.get(finalI, j)).sum();
            prob.set(i, 0, FastMath.abs(1 - rowSum));
        }
        return prob;
    }
}
