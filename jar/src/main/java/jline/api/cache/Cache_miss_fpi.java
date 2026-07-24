/**
 * @file Cache Miss Analysis via Fixed Point Iteration
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Cache_miss_fpi {
    private Cache_miss_fpi() {}

    public static CacheMissFpiResult cache_miss_fpi(Matrix gamma, Matrix m) {
        return cache_miss_fpi(gamma, m, null);
    }

    /**
     * Compute cache miss rates using Fixed Point Iteration (FPI) method.
     */
    public static CacheMissFpiResult cache_miss_fpi(Matrix gamma, Matrix m, MatrixCell lambda) {
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();

        Matrix xi = Cache_xi_fp.cache_xi_fp(gamma, m, null).xi;

        Matrix pi0 = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            double gammaXiSum = 0.0;
            for (int j = 0; j < h; j++) {
                gammaXiSum += gamma.get(i, j) * xi.get(0, j);
            }
            pi0.set(i, 0, 1.0 / (1.0 + gammaXiSum));
        }

        if (lambda == null) {
            double M = pi0.elementSum();
            return new CacheMissFpiResult(M, null, null, pi0);
        }

        int u = lambda.size();

        Matrix MI = new Matrix(n, 1);
        for (int i = 0; i < n; i++) {
            double lambdaSum = 0.0;
            for (int v = 0; v < u; v++) {
                MatrixCell userMatrix = (MatrixCell) (Object) lambda.get(v);
                lambdaSum += userMatrix.get(i).get(0, 0);
            }
            double gammaXiSum = 0.0;
            for (int j = 0; j < h; j++) {
                gammaXiSum += gamma.get(i, j) * xi.get(0, j);
            }
            MI.set(i, 0, lambdaSum / (1.0 + gammaXiSum));
        }

        Matrix MU = new Matrix(u, 1);
        for (int v = 0; v < u; v++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) {
                MatrixCell userMatrix = (MatrixCell) (Object) lambda.get(v);
                sum += userMatrix.get(i).get(0, 0) * pi0.get(i, 0);
            }
            MU.set(v, 0, sum);
        }

        double M = MI.elementSum();
        return new CacheMissFpiResult(M, MU, MI, pi0);
    }
}
