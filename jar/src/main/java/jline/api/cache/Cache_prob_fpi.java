/**
 * @file Cache Probability Analysis via Fixed Point Iteration
 *
 * Computes cache state probabilities using fixed-point iteration algorithms.
 * Provides iterative solution methods for cache systems where direct
 * analytical solutions are not feasible due to system complexity.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Cache_prob_fpi {
    private Cache_prob_fpi() {}

    /**
     * Estimate asymptotic values of the cache state probabilities at steady-state.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m     Matrix representing the cache capacity vector.
     * @return A matrix containing the estimated cache state probabilities at equilibrium.
     */
    public static Matrix cache_prob_fpi(Matrix gamma, Matrix m) {
        // FPI method
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();

        Ret.cacheXiFp result = Cache_xi_fp.cache_xi_fp(gamma, m, null);
        Matrix xi = result.xi;
        Matrix prob = new Matrix(n, h + 1);
        for (int i = 0; i < n; i++) {
            Matrix mul = Matrix.extractRows(gamma, i, i + 1, null).mult(xi.columnMajorOrder(), null);
            prob.set(i, 0, 1 / (1 + mul.get(0)));
            for (int j = 1; j < h + 1; j++) {
                prob.set(i, j, mul.get(0) / (1 + mul.get(0)));
            }
        }
        return prob;
    }
}
