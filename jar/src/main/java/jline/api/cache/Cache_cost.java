/**
 * @file Mean per-list storage cost of a cache with item sizes
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class Cache_cost {
    private Cache_cost() {}

    /**
     * Computes K_j = sum_i sigma_i pi_ij, the expected storage cost of the items
     * resident in each list at steady state, as defined in Casale-Gast, IEEE/ACM
     * Trans. Networking 29(2), 2021, Sec. IX.
     *
     * @param gamma Cache access factors (n x h).
     * @param m Cache capacity vector (1 x h).
     * @param sigma Item storage costs (sizes), one per item.
     * @param k Per-list storage cost caps; null or empty for none.
     * @return the mean storage cost of each list (1 x h).
     */
    public static Matrix cache_cost(Matrix gamma, Matrix m, Matrix sigma, Matrix k) {
        return cache_cost(gamma, m, sigma, k, null);
    }

    /**
     * Computes the mean per-list storage cost from a precomputed occupancy matrix.
     *
     * @param gamma Cache access factors (n x h).
     * @param m Cache capacity vector (1 x h).
     * @param sigma Item storage costs (sizes), one per item.
     * @param k Per-list storage cost caps; null or empty for none.
     * @param pij Occupancy matrix (n x (h+1)) with column 0 the miss probability; null to recompute.
     * @return the mean storage cost of each list (1 x h).
     */
    public static Matrix cache_cost(Matrix gamma, Matrix m, Matrix sigma, Matrix k, Matrix pij) {
        int n = gamma.getNumRows();
        int h = gamma.getNumCols();
        if (sigma == null || sigma.length() != n) {
            InputOutput.line_error(InputOutput.mfilename(new Object()),
                    "The item size vector must have one entry per item.");
        }
        Matrix occupancy = pij;
        if (occupancy == null) {
            occupancy = Cache_prob_erec.cache_prob_erec(gamma, m, sigma, k);
        }
        Matrix K = new Matrix(1, h);
        for (int j = 0; j < h; j++) {
            double c = 0.0;
            for (int i = 0; i < n; i++) {
                c += sigma.get(i) * occupancy.get(i, j + 1);
            }
            K.set(0, j, c);
        }
        return K;
    }
}
