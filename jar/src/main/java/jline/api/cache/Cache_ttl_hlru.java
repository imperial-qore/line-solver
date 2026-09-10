/**
 * @file TTL (characteristic-time) approximation for h-LRU / LRU(m) caches
 *
 * Steady-state list occupancy probabilities of the list-based h-LRU (LRU(m))
 * replacement policy under the characteristic-time approximation (Gast and
 * Van Houdt, SIGMETRICS 2015). For h=1 this reduces exactly to the Che
 * approximation for LRU.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.util.matrix.Matrix;

public final class Cache_ttl_hlru {
    private Cache_ttl_hlru() {}

    /**
     * Steady-state list occupancy probabilities for an h-LRU cache.
     *
     * @param lambda per-class request rates; lambda[v] is (n x h+1) with the
     *               same per-item rate replicated across list slices (the MVA
     *               cache analyzer layout); rates are aggregated over classes
     * @param m      (1 x h) list capacities
     * @return (n x h+1) probabilities; column 0 = not cached, column 1+l = in
     *         list l
     */
    public static Matrix cache_ttl_hlru(Matrix[] lambda, Matrix m) {
        int n = lambda[0].getNumRows();
        int h = m.length();
        int rateCol = Math.min(1, lambda[0].getNumCols() - 1);
        Matrix lam = new Matrix(n, 1);
        for (int v = 0; v < lambda.length; v++) {
            for (int k = 0; k < n; k++) {
                lam.set(k, 0, lam.get(k, 0) + lambda[v].get(k, rateCol));
            }
        }

        Matrix t = Cache_t_hlru.cache_t_hlru(lam, m);
        double[] tv = new double[h];
        for (int l = 0; l < h; l++) {
            tv[l] = t.get(0, l);
        }

        Matrix pij = new Matrix(n, h + 1);
        for (int k = 0; k < n; k++) {
            double[] w = Cache_t_hlru.levelWeights(lam.get(k, 0), tv, h);
            double sum = 0.0;
            for (int s = 0; s <= h; s++) {
                sum += w[s];
            }
            for (int s = 0; s <= h; s++) {
                pij.set(k, s, w[s] / sum);
            }
        }
        return pij;
    }
}
