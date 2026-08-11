/**
 * @file Exact Recursive Cache Probability Computation
 *
 * Computes exact cache state probabilities using recursive methods based on
 * exact recursive (EREC) algorithms. Provides precise analysis for small to
 * medium-sized cache systems where computational complexity is manageable.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.Collections;
import java.util.stream.IntStream;

public final class Cache_prob_erec {
    private Cache_prob_erec() {}

    /**
     * Computes the cache state probabilities using an exact recursive method.
     * This method calculates the probabilities of the cache being in different states based on the cache
     * access factors and capacity.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m Matrix representing the cache capacity vector.
     * @return matrix containing the computed cache state probabilities for each item and level.
     */
    public static Matrix cache_prob_erec(Matrix gamma, Matrix m) {
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        Matrix E = Cache_erec.cache_erec(gamma, m);
        Matrix prob = new Matrix(n, h + 1);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < h; j++) {
                Matrix subGamma = gamma.copy();
                subGamma.removeRows(Collections.singleton(i));
                Matrix Ei = Cache_erec.cache_erec(subGamma, Matrix.oner(m, j));
                double value = m.get(j) * gamma.get(i, j) * Ei.value() / E.value();
                prob.set(i, j + 1, value);
            }
            final int finalI = i;
            double rowSum = IntStream.range(1, h + 1).mapToDouble(j -> prob.get(finalI, j)).sum();
            prob.set(i, 0, FastMath.abs(1 - rowSum));
        }
        return prob;
    }
}
