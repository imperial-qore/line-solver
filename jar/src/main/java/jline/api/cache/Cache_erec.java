/**
 * @file Exact Recursive Cache Analysis
 *
 * Implements exact recursive (EREC) algorithms for cache system analysis.
 * Provides numerically exact solutions for cache performance metrics in
 * systems where computational complexity allows for precise evaluation.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.util.matrix.Matrix;

public final class Cache_erec {
    private Cache_erec() {}

    /**
     * Computes the cache miss rate using an exact recursive method.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m     Matrix representing the cache capacity vector.
     * @return Matrix containing the computed cache miss rates.
     */
    public static Matrix cache_erec(Matrix gamma, Matrix m) {
        return cache_erec_aux(gamma, m, gamma.getNumRows());
    }

    /**
     * Auxiliary method for computing the cache miss rate using an exact recursive method.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m     Matrix representing the cache capacity vector.
     * @param k     Integer representing the current number of rows in the recursive step.
     * @return Matrix containing the computed cache miss rates.
     */
    public static Matrix cache_erec_aux(Matrix gamma, Matrix m, int k) {
        int h = m.getNumElements();

        if (m.elementSum() == 0.0) {
            return Matrix.singleton(1.0);
        }

        if (m.elementSum() > k || m.elementMin() < 0) {
            return Matrix.singleton(0.0);
        }

        if (k == 1 && m.elementSum() == 1.0) {
            int j = (int) m.find().value(); // Find the index of the non-zero element in m
            return Matrix.singleton(gamma.get(0, j));
        }

        Matrix E = cache_erec_aux(gamma, m, k - 1);
        for (int j = 0; j < h; j++) {
            if (m.get(j) > 0) {
                Matrix onerM = Matrix.oner(m, j);
                Matrix term = cache_erec_aux(gamma, onerM, k - 1).scale(gamma.get(k - 1, j) * m.get(j));
                E = E.add(term);
            }
        }
        return E;
    }
}
