/**
 * @file General Cache Miss Rate Analysis
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.util.matrix.Matrix;

public final class Cache_miss {
    private Cache_miss() {}

    /**
     * Result of cache_miss analysis.
     */
    public static final class CacheMissResult {
        public final double globalMissRate;
        public final Matrix perUserMissRate;
        public final Matrix perItemMissRate;
        public final Matrix perItemMissProb;

        public CacheMissResult(double globalMissRate, Matrix perUserMissRate, Matrix perItemMissRate, Matrix perItemMissProb) {
            this.globalMissRate = globalMissRate;
            this.perUserMissRate = perUserMissRate;
            this.perItemMissRate = perItemMissRate;
            this.perItemMissProb = perItemMissProb;
        }

        public CacheMissResult(double globalMissRate) {
            this(globalMissRate, null, null, null);
        }
    }

    public static CacheMissResult cache_miss(Matrix gamma, Matrix m, Matrix lambda) {
        Matrix ma = m.copy();
        ma.set(0, ma.get(0) + 1.0);

        double globalMissRate = Cache_erec.cache_erec(gamma, ma).get(0)
                / Cache_erec.cache_erec(gamma, m).get(0);

        if (lambda.isEmpty()) {
            return new CacheMissResult(globalMissRate);
        }

        int u = lambda.getNumRows();
        int n = lambda.getNumCols();

        Matrix pi0 = Matrix.zeros(1, n);
        double denominatorValue = Cache_erec.cache_erec(gamma, m).get(0);

        for (int k = 0; k < n; k++) {
            Matrix gammaWithoutK = Matrix.zeros(gamma.getNumRows(), gamma.getNumCols() - 1);
            int colIndex = 0;
            for (int j = 0; j < gamma.getNumCols(); j++) {
                if (j != k) {
                    for (int i = 0; i < gamma.getNumRows(); i++) {
                        gammaWithoutK.set(i, colIndex, gamma.get(i, j));
                    }
                    colIndex++;
                }
            }

            double numeratorValue = Cache_erec.cache_erec(gammaWithoutK, m).get(0);
            pi0.set(0, k, numeratorValue / denominatorValue);
        }

        Matrix MU = Matrix.zeros(u, 1);
        for (int v = 0; v < u; v++) {
            for (int k = 0; k < n; k++) {
                MU.set(v, 0, MU.get(v, 0) + lambda.get(v, k) * pi0.get(0, k));
            }
        }

        Matrix MI = Matrix.zeros(n, 1);
        for (int k = 0; k < n; k++) {
            double sum = 0.0;
            for (int v = 0; v < u; v++) {
                sum += lambda.get(v, k);
            }
            MI.set(k, 0, sum * pi0.get(0, k));
        }

        return new CacheMissResult(globalMissRate, MU, MI, pi0);
    }

    public static CacheMissResult cache_miss(Matrix gamma, Matrix m) {
        return cache_miss(gamma, m, new Matrix(gamma.getNumRows(), gamma.getNumCols()));
    }
}
