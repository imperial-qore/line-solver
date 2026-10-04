/**
 * @file Marked MAP counting-process power moments
 *
 * Computes per-class power moments of the counting process of a Marked MAP by
 * reducing each class to its marginal MAP and applying {@link Map_count_moment}.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_count_moment {
    private Mmap_count_moment() {}

    /**
     * Per-class power moments of counts, at resolution t, of a Marked MAP.
     *
     * <p>For each class k the marginal MAP is {@code {D0 + sum_{j!=k} D1_j, D1_k}}
     * (all non-class-k transitions are folded into the hidden generator), and its
     * count moments are obtained from {@link Map_count_moment#map_count_moment}.</p>
     *
     * @param MMAP   marked MAP: {D0, D1, D1_1, ..., D1_K}
     * @param t      resolution
     * @param orders orders of the moments to compute
     * @return matrix of size (orders.length x K); column k holds the moments for class k
     */
    public static Matrix mmap_count_moment(MatrixCell MMAP, double t, int[] orders) {
        int K = MMAP.size() - 2;
        Matrix D0 = MMAP.get(0);
        int n = D0.getNumRows();
        Matrix out = new Matrix(orders.length, K);

        for (int k = 0; k < K; k++) {
            // Marginal generator D0' = D0 + sum_{j != k} D1_j
            Matrix D0marg = D0.copy();
            for (int j = 0; j < K; j++) {
                if (j != k) {
                    D0marg = D0marg.add(MMAP.get(2 + j));
                }
            }
            Matrix D1marg = MMAP.get(2 + k);
            MatrixCell marginal = new MatrixCell(D0marg, D1marg);
            double[] m = Map_count_moment.map_count_moment(marginal, t, orders);
            for (int i = 0; i < orders.length; i++) {
                out.set(i, k, m[i]);
            }
        }
        return out;
    }
}
