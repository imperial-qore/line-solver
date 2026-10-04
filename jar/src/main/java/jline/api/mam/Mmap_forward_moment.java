/**
 * @file Marked Markovian Arrival Process forward moment analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.CombinatoricsUtils;

public final class Mmap_forward_moment {
    private Mmap_forward_moment() {}

    /**
     * Computes the forward moments of an MMAP for specified orders with normalization.
     *
     * @param MMAP   the MMAP from which to compute the forward moments
     * @param ORDERS a matrix specifying the orders of the moments to be computed
     * @param NORM   normalization flag
     * @return a matrix where each row corresponds to a type and each column corresponds to a moment order
     */
    public static Matrix mmap_forward_moment(MatrixCell MMAP, Matrix ORDERS, int NORM) {
        int C = MMAP.size() - 2;
        int K = ORDERS.length();

        Matrix MOMENTS = new Matrix(C, K, C * K);
        Matrix pie = Map_pie.map_pie(MMAP.get(0), MMAP.get(1));

        Matrix neg_D0 = MMAP.get(0).copy();
        neg_D0.scaleEq(-1.0);
        Matrix M = neg_D0.inv();

        for (int a = 0; a < C; a++) {
            double pa;
            if (NORM == 1) {
                pa = pie.mult(M.mult(MMAP.get(2 + a))).elementSum();
            } else {
                pa = 1.0;
            }
            for (int h = 0; h < ORDERS.length(); h++) {
                int k = (int) ORDERS.get(h);
                double fk = (double) CombinatoricsUtils.factorial(k);
                MOMENTS.set(a, h, fk / pa * pie.mult(M.mult(MMAP.get(2 + a))).mult(Matrix.pow(M, k)).elementSum());
            }
        }
        return MOMENTS;
    }

    public static Matrix mmap_forward_moment(MatrixCell MMAP, Matrix ORDERS) {
        return mmap_forward_moment(MMAP, ORDERS, 1);
    }
}
