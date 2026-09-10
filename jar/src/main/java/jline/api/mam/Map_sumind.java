/**
 * @file Markovian Arrival Process independent summation operations
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_sumind {
    private Map_sumind() {}

    /**
     * Computes the Markovian Arrival Process (MAP) representing the sum of `n` independent MAPs.
     *
     * @param MAPs an array of MatrixCells, each containing transition matrices D0 and D1
     * @return a MatrixCell containing the transition matrices D0 and D1 of the resulting summed MAP
     */
    public static MatrixCell map_sumind(MatrixCell[] MAPs) {
        int n = MAPs.length;
        Matrix order = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            order.set(i, (double) MAPs[i].get(0).getNumRows());
        }
        int total = (int) order.elementSum();
        Matrix D0 = Matrix.zeros(total, total);
        Matrix D1 = Matrix.zeros(total, total);
        int curpos = 0;
        for (int i = 0; i < n; i++) {
            int oi = (int) order.get(i);
            D0.setSliceEq(curpos, curpos + oi, curpos, curpos + oi, MAPs[i].get(0));
            if (i < (n - 1)) {
                int oi1 = (int) order.get(i + 1);
                D0.setSliceEq(curpos,
                        curpos + oi,
                        curpos + oi,
                        curpos + oi + oi1,
                        MAPs[i].get(1).mult(Matrix.ones(oi, 1)).mult(Map_pie.map_pie(MAPs[i + 1])));
            } else {
                D1.setSliceEq(curpos,
                        curpos + oi,
                        0,
                        oi,
                        MAPs[i].get(1).mult(Matrix.ones(oi, 1)).mult(Map_pie.map_pie(MAPs[i])));
            }
            curpos += oi;
        }
        return new MatrixCell(D0, D1);
    }
}
