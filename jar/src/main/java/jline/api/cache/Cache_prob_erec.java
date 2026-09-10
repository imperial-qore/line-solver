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
        return cache_prob_erec(gamma, m, null, null);
    }

    /**
     * Computes the cache state probabilities under per-list storage cost caps,
     * pi_ij = m_j gamma_ij E_i(m-1_j, k-sigma_i 1_j) / E(m,k). See Casale-Gast,
     * IEEE/ACM Trans. Networking 29(2), 2021, Sec. IX.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m Matrix representing the cache capacity vector.
     * @param sigma Item storage costs (sizes); null or empty for none.
     * @param k Per-list storage cost caps; null or empty for none.
     * @return matrix containing the computed cache state probabilities for each item and level.
     */
    public static Matrix cache_prob_erec(Matrix gamma, Matrix m, Matrix sigma, Matrix k) {
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        boolean capped = sigma != null && k != null && !sigma.isEmpty() && !k.isEmpty();
        Matrix E = Cache_erec.cache_erec(gamma, m, sigma, k);
        Matrix prob = new Matrix(n, h + 1);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < h; j++) {
                Matrix subGamma = gamma.copy();
                subGamma.removeRows(Collections.singleton(i));
                double value;
                if (capped) {
                    Matrix kij = k.copy();
                    kij.set(j, kij.get(j) - sigma.get(i));
                    if (kij.get(j) < 0) {
                        value = 0.0;
                    } else {
                        Matrix Ei = Cache_erec.cache_erec(subGamma, Matrix.oner(m, j), dropEntry(sigma, i), kij);
                        value = m.get(j) * gamma.get(i, j) * Ei.value() / E.value();
                    }
                } else {
                    Matrix Ei = Cache_erec.cache_erec(subGamma, Matrix.oner(m, j));
                    value = m.get(j) * gamma.get(i, j) * Ei.value() / E.value();
                }
                prob.set(i, j + 1, value);
            }
            final int finalI = i;
            double rowSum = IntStream.range(1, h + 1).mapToDouble(j -> prob.get(finalI, j)).sum();
            prob.set(i, 0, FastMath.abs(1 - rowSum));
        }
        return prob;
    }

    /**
     * Returns the vector with the entry at the given position removed.
     */
    static Matrix dropEntry(Matrix v, int i) {
        int len = v.length();
        Matrix out = new Matrix(1, len - 1);
        int c = 0;
        for (int t = 0; t < len; t++) {
            if (t != i) {
                out.set(0, c, v.get(t));
                c++;
            }
        }
        return out;
    }
}
