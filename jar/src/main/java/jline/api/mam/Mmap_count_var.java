/**
 * @file Marked Markovian Arrival Process counting process variance analysis
 *
 * Computes variance vectors for counting processes of each marked class in MMAP.
 * Essential for analyzing variability and dispersion in multiclass arrival systems.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.Maths;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_count_var {
    private Mmap_count_var() {}

    /**
     * Computes the variance of the count vector of events of different types in a MMAP over a time period.
     *
     * @param MMAP the MatrixCell containing the transition matrices of the MMAP, with D0, D1, ..., Dc representing different types of events
     * @param t    the time period over which to compute the variances of counts
     * @return a Matrix containing the variance vector of counts over time t
     */
    public static Matrix mmap_count_var(MatrixCell MMAP, double t) {
        Matrix D0 = MMAP.get(0);
        Matrix D1 = MMAP.get(1);
        Matrix theta = Map_piq.map_piq(D0, D1);
        int C = MMAP.size() - 2;
        Matrix vt = new Matrix(1, C);
        Matrix Q = Map_infgen.map_infgen(D0, D1);
        int n = D0.getNumRows();
        Matrix I = Matrix.eye(n);
        Matrix e = Matrix.ones(n, 1);
        Matrix tmp = Matrix.ones(n, 1).mult(theta).add(-1.0, Q).inv();
        Matrix lk = new Matrix(1, C);
        MatrixCell ck = new MatrixCell(C);
        MatrixCell dk = new MatrixCell(C);
        Matrix llk = new Matrix(1, C);
        for (int c = 0; c < C; c++) {
            Matrix Dc = MMAP.get(2 + c);
            lk.set(c, theta.mult(Dc.mult(e)).value());
            ck.set(c, theta.mult(Dc.mult(tmp)));
            dk.set(c, tmp.mult(Dc.mult(e)));
            llk.set(c, theta.mult(Dc.mult(e)).value());
        }
        for (int c = 0; c < C; c++) {
            Matrix Dc = MMAP.get(2 + c);
            double vc = llk.get(c) - 2 * lk.get(c) * lk.get(c);
            vc += ck.get(c).scale(2.0).mult(Dc.mult(e)).value();
            vc *= t;
            vc -= 2 * ck.get(c).mult(I.sub(Maths.matrixExp(Q.scale(t)))).mult(dk.get(c)).value();
            vt.set(c, vc);
        }

        return vt;
    }
}
