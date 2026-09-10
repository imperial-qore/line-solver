/**
 * @file Normalizing Constant Approximation (NCA) for single-class networks
 *
 * Implements the Normalizing Constant Approximation method for single-class closed
 * queueing networks with load-dependent service rates. Uses iterative computation
 * with service rate functions to approximate the normalizing constant.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import jline.util.matrix.Matrix;

public final class Pfqn_nca {
    private Pfqn_nca() {}

    public static Matrix pfqn_nca(Matrix L, int N, Matrix alpha) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (R > 1 && M > 1) {
            throw new RuntimeException("pfqn_nca supports only single class queuing networks; L must be a column vector.");
        }
        Matrix X = Matrix.ones(M, N);
        Matrix pi = new Matrix(M, 1 + N);
        pi.zero();
        for (int j = 0; j < N; j++) {
            X.set(0, j, alpha.get(0, j) / L.get(0));
        }
        for (int i = 1; i < M; i++) {
            pi.set(i, 0, 1.0);
            for (int k = 0; k < N; k++) {
                for (int j = 0; j < k; j++) {
                    pi.set(i, 1 + j, pi.get(i - 1, j) * L.get(i) / alpha.get(i, j));
                }
                pi.set(i, 1, pi.get(i, 1) / X.get(i - 1, k));
                X.set(i, k, 1 / pi.sumSubMatrix(i, i, 1, k + 1));
                for (int j = -1; j < k; j++) {
                    pi.set(i, j + 1, X.get(i, k) * pi.get(i, j + 1));
                }
            }
        }
        return X;
    }
}
