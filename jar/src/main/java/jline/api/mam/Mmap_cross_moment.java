/**
 * @file Marked Markovian Arrival Process cross-moment analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

public final class Mmap_cross_moment {
    private Mmap_cross_moment() {}

    /**
     * Computes the k-th cross-moment matrix for a given MMAP.
     *
     * @param mmap the MMAP for which to compute the cross-moments
     * @param k    the order of the moment to compute
     * @return a matrix representing the k-th cross-moment for each pair of types
     */
    public static Matrix mmap_cross_moment(MatrixCell mmap, int k) {
        int C = mmap.size() - 2;

        MatrixCell TG = new MatrixCell();
        Matrix MC = new Matrix(C, C, (int) FastMath.pow((double) C, 2));

        for (int i = 0; i < C; i++) {
            Matrix a = mmap.get(0).inv().mult(mmap.get(2 + i));
            a.scaleEq(-1.0);
            TG.set(i, Map_pie.map_pie(mmap.get(0), mmap.get(1)).mult(a).sumCols());
        }

        for (int i = 0; i < C; i++) {
            Matrix a = mmap.get(0).inv().mult(mmap.get(2 + i));
            a.scaleEq(-1.0);
            Matrix start = Map_pie.map_pie(mmap.get(0), mmap.get(1)).mult(a).mult(TG.get(i).inv());
            for (int j = 0; j < C; j++) {
                Matrix neg_mmap0 = mmap.get(0).copy();
                neg_mmap0.scaleEq(-1.0);
                MC.set(i, j,
                        CombinatoricsUtils.factorial(k)
                                * start.mult(Matrix.pow(neg_mmap0.inv(), k + 1).mult(mmap.get(2 + j))).elementSum());
                MC.set(i, j, MC.get(i, j) / start.mult(neg_mmap0.mult(mmap.get(2 + j).inv())).elementSum());
            }
        }
        return MC;
    }
}
