/**
 * @file Marked Markovian Arrival Process counting covariance analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_count_mcov {
    private Mmap_count_mcov() {}

    /**
     * Computes the count covariance between each pair of classes at a given time scale.
     */
    public static Matrix mmap_count_mcov(MatrixCell MMAP, double t) {
        int m = MMAP.size() - 2;

        Matrix mV = Mmap_count_var.mmap_count_var(MMAP, t);
        Matrix S = Matrix.zeros(m, m);

        for (int i = 0; i < m; i++) {
            S.set(i, i, mV.get(0, i));
        }

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                if (i != j) {
                    MatrixCell mmap2 = new MatrixCell(4);
                    mmap2.set(0, MMAP.get(0));
                    mmap2.set(1, MMAP.get(1));
                    mmap2.set(2, MMAP.get(2 + i).add(MMAP.get(2 + j)));
                    mmap2.set(3, mmap2.get(1).sub(mmap2.get(2)));

                    Matrix pV = Mmap_count_var.mmap_count_var(mmap2, t);

                    double covariance = 0.5 * (pV.get(0, 0) - mV.get(0, i) - mV.get(0, j));
                    S.set(i, j, covariance);
                }
            }
        }

        return S;
    }

    /**
     * Array-based overload.
     */
    public static Matrix mmap_count_mcov(Matrix[] mmap, double t) {
        MatrixCell mmapCell = new MatrixCell(mmap.length);
        for (int i = 0; i < mmap.length; i++) {
            mmapCell.set(i, mmap[i]);
        }
        return mmap_count_mcov(mmapCell, t);
    }
}
