/**
 * @file TTL approximation of LRU(m) caches under MAP request streams
 *
 * Front-end of the Gast-Van Houdt (Performance Evaluation 2017) TTL
 * approximation for LRU(m) with per-item MAP request processes. The
 * characteristic times are computed by Cache_t_lrum_map; this class exposes
 * the per-item request-weighted hit/miss probabilities and the aggregate
 * hit rate. Intended for caches whose items have genuinely distinct or
 * correlated request processes (e.g. marked MAP arrivals); when the items
 * are i.i.d. marks of a common stream the request sequence is IRM and the
 * Poisson-based TTL approximations already apply.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Cache_ttl_lrum_map {
    private Cache_ttl_lrum_map() {}

    /**
     * Per-item request-weighted probabilities under the LRU(m)-MAP TTL
     * approximation.
     *
     * @param D0Matrix per-item MAP hidden-transition cells; D0Matrix[k].get(l)
     *                 is the D0 matrix of item k while in list l (l=0..h-1)
     * @param D1Matrix per-item MAP arrival cells, same layout
     * @param m        cache capacity vector (1 x h)
     * @return (n x h+1) matrix; column 0 is the probability that a request
     *         for item k misses, column l the probability that it hits in
     *         list l.
     */
    public static Matrix cache_ttl_lrum_map_pij(MatrixCell[] D0Matrix, MatrixCell[] D1Matrix, Matrix m) {
        int n = D0Matrix.length;
        int h = D0Matrix[0].size();

        Matrix tMatrix = Cache_t_lrum_map.cache_t_lrum_map(D0Matrix, D1Matrix, m);
        double[] t = new double[(int) tMatrix.length()];
        for (int ti = 0; ti < t.length; ti++) {
            t[ti] = tMatrix.get(ti);
        }

        Matrix pij = new Matrix(n, h + 1);
        for (int k = 0; k < n; k++) {
            Cache_t_lrum_map.LevelStats st = Cache_t_lrum_map.levelStats(t, D0Matrix[k], D1Matrix[k], h);
            double hitsum = 0.0;
            for (int l = 0; l < h; l++) {
                pij.set(k, l + 1, st.hitfrac[l]);
                hitsum += st.hitfrac[l];
            }
            pij.set(k, 0, Math.max(0.0, 1.0 - hitsum));
        }
        return pij;
    }

    /**
     * Aggregate request-weighted hit rate under the LRU(m)-MAP TTL
     * approximation.
     */
    public static double cache_ttl_lrum_map(MatrixCell[] D0Matrix, MatrixCell[] D1Matrix, Matrix m) {
        int n = D0Matrix.length;
        int d = D0Matrix[0].get(0).getNumRows();
        Matrix pij = cache_ttl_lrum_map_pij(D0Matrix, D1Matrix, m);

        double hits = 0.0;
        double lamsum = 0.0;
        for (int k = 0; k < n; k++) {
            double lam = Cache_t_lrum_map.mapRate(D0Matrix[k].get(0), D1Matrix[k].get(0), d);
            double hitk = 0.0;
            for (int l = 1; l <= pij.getNumCols() - 1; l++) {
                hitk += pij.get(k, l);
            }
            hits += lam * hitk;
            lamsum += lam;
        }
        return lamsum > 0 ? hits / lamsum : 0.0;
    }
}
