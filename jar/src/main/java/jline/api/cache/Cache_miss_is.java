/**
 * @file Cache Miss Analysis via Importance Sampling
 *
 * Computes cache miss rates using Monte Carlo importance sampling.
 * Provides global, per-user, and per-item miss rate estimates.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Cache_miss_is {
    private Cache_miss_is() {}

    /**
     * Computes cache miss rates using Monte Carlo importance sampling.
     *
     * @param gamma Matrix representing the cache access factors (n x h).
     * @param m Matrix representing the cache capacity vector (1 x h).
     * @param lambda MatrixCell representing the request rates for different users or items.
     * @param samples Number of Monte Carlo samples.
     * @return Ret.cacheMissSpm containing miss rate metrics (M, MU, MI, pi0, lE).
     */
    public static Ret.cacheMissSpm cache_miss_is(Matrix gamma, Matrix m, MatrixCell lambda, int samples) {
        return cache_miss_is(gamma, m, lambda, samples, null, null);
    }

    /**
     * Computes cache miss rates by importance sampling, optionally under
     * per-list storage cost caps.
     *
     * @param gamma Cache access factors (n x h).
     * @param m Cache capacity vector (1 x h).
     * @param lambda Request rates per user per item.
     * @param samples Number of Monte Carlo samples.
     * @param sigma Item storage costs (sizes); null or empty for none.
     * @param cap Per-list storage cost caps; null or empty for none.
     * @return miss rate metrics (M, MU, MI, pi0, lE).
     */
    public static Ret.cacheMissSpm cache_miss_is(Matrix gamma, Matrix m, MatrixCell lambda, int samples, Matrix sigma, Matrix cap) {
        int n = gamma.getNumRows();

        // Compute normalizing constant via importance sampling
        Ret.cacheIs isResult = Cache_is.cache_is(gamma, m, samples, sigma, cap);
        double lE = isResult.lE;

        // Compute hit probabilities via importance sampling
        Matrix pij = Cache_prob_is.cache_prob_is(gamma, m, samples, sigma, cap);

        // Extract miss probabilities (first column)
        double[] pi0 = new double[n];
        for (int i = 0; i < n; i++) {
            pi0[i] = pij.get(i, 0);
        }

        int u = lambda.size();

        // Per-user miss rate
        double[] MU = new double[u];
        for (int v = 0; v < u; v++) {
            for (int k = 0; k < n; k++) {
                MU[v] += lambda.get(v).get(k, 0) * pi0[k];
            }
        }

        // Per-item miss rate
        double[] MI = new double[n];
        for (int k = 0; k < n; k++) {
            MI[k] = lambda.cellsum(k, 0) * pi0[k];
        }

        // Global miss rate (sum of per-item miss rates)
        double M = 0.0;
        for (int k = 0; k < n; k++) {
            M += MI[k];
        }

        return new Ret.cacheMissSpm(M, MU, MI, pi0, lE);
    }

    public static Ret.cacheMissSpm cache_miss_is(Matrix gamma, Matrix m, MatrixCell lambda) {
        return cache_miss_is(gamma, m, lambda, 100000);
    }
}
