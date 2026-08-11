/**
 * @file Load-dependent scaling vector position shifting utility
 *
 * Provides utility function for shifting load-dependent scaling vectors by one position,
 * commonly used in recursive normalizing constant computations. Supports the manipulation
 * of state-dependent service rate vectors during iterative solution algorithms.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import jline.util.matrix.Matrix;

public final class Pfqn_mushift {
    private Pfqn_mushift() {}

    /**
     * Shifts a load-dependent scaling vector by one position
     *
     * @param mu - load-dependent scalings
     * @return normalizing constant and its logarithm
     */
    public static Matrix pfqn_mushift(Matrix mu, int k) {
        int M = mu.getNumRows();
        int N = mu.getNumCols();
        Matrix mushift = new Matrix(M, N - 1);
        Matrix.extract(mu, 0, M, 0, N - 1, mushift, 0, 0);
        for (int j = 0; j < N - 1; j++) {
            mushift.set(k, j, mu.get(k, j + 1));
        }
        return mushift;
    }
}
