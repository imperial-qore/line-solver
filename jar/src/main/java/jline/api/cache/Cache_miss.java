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
        return cache_miss(gamma, m, lambda, null, null);
    }

    /**
     * Computes cache miss metrics, optionally under per-list storage cost caps.
     *
     * @param gamma Cache access factors (n x h).
     * @param m Cache capacity vector (1 x h).
     * @param lambda Arrival rates per user per item (u x n).
     * @param sigma Item storage costs (sizes); null or empty for none.
     * @param cap Per-list storage cost caps; null or empty for none.
     * @return the miss metrics.
     */
    public static CacheMissResult cache_miss(Matrix gamma, Matrix m, Matrix lambda, Matrix sigma, Matrix cap) {
        Matrix ma = m.copy();
        ma.set(0, ma.get(0) + 1.0);

        double denominatorValue = Cache_erec.cache_erec(gamma, m, sigma, cap).get(0);
        double globalMissRate = Cache_erec.cache_erec(gamma, ma, sigma, cap).get(0) / denominatorValue;

        if (lambda.isEmpty()) {
            return new CacheMissResult(globalMissRate);
        }

        int u = lambda.getNumRows();
        int n = lambda.getNumCols();

        Matrix pi0 = Matrix.zeros(1, n);

        for (int k = 0; k < n; k++) {
            // E_k is the constant of the model WITHOUT item k, i.e. gamma row k dropped
            Matrix gammaWithoutK = Matrix.zeros(gamma.getNumRows() - 1, gamma.getNumCols());
            int rowIndex = 0;
            for (int i = 0; i < gamma.getNumRows(); i++) {
                if (i != k) {
                    for (int j = 0; j < gamma.getNumCols(); j++) {
                        gammaWithoutK.set(rowIndex, j, gamma.get(i, j));
                    }
                    rowIndex++;
                }
            }
            Matrix sigmaWithoutK = (sigma == null || sigma.isEmpty()) ? sigma : Cache_prob_erec.dropEntry(sigma, k);

            double numeratorValue = Cache_erec.cache_erec(gammaWithoutK, m, sigmaWithoutK, cap).get(0);
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
