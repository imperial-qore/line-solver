/**
 * @file Marked Markovian Arrival Process n-fold interarrival sum
 *
 * Computes the MMAP whose interarrival time is the sum of n interarrival times
 * of a base MMAP (n-fold convolution). This is the operation MATLAB and the
 * native-Python line_solver name mmap_sum; the superposition of independent
 * streams is a distinct operation provided by {@link Mmap_super}.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_sum {
    private Mmap_sum() {}

    /**
     * Returns the MMAP representing the sum of n identical copies of the given
     * MMAP: its interarrival time is distributed as the sum of n interarrival
     * times of the base process, and the class marking is inherited from the
     * n-th (last) of the n base arrivals.
     *
     * <p>Mirrors MATLAB mmap_sum(MMAP,n) and native-Python mmap_sum(mmap,n).
     * The MMAP is a {@link MatrixCell} laid out as
     * {D0, D1(aggregate), C_1, ..., C_K}.
     *
     * @param mmap the base MMAP (K classes, order m)
     * @param n    the number of interarrival times to sum (n &gt;= 1)
     * @return the n-fold summed MMAP (K classes, order n*m)
     */
    public static MatrixCell mmap_sum(MatrixCell mmap, int n) {
        if (n < 1) {
            throw new IllegalArgumentException("n must be at least 1");
        }
        int K = mmap.size() - 2;            // number of classes
        int m = mmap.get(0).getNumRows();   // order of each copy
        int total = n * m;

        Matrix d0 = new Matrix(total, total);
        // Output cells 0..K of D1: index 0 is the aggregate D1, 1..K are the
        // per-class matrices (matching the base MMAP cells 1..K+1)
        Matrix[] d1 = new Matrix[K + 1];
        for (int k = 0; k <= K; k++) {
            d1[k] = new Matrix(total, total);
        }

        int curpos = 0;
        for (int i = 0; i < n; i++) {
            // Diagonal block: the copy's own D0 (no arrival emitted)
            d0.insertSubMatrix(curpos, curpos, curpos + m, curpos + m, mmap.get(0));
            if (i < n - 1) {
                // Hidden interarrival: on a base arrival, advance to the next
                // copy without emitting an arrival in the summed process
                d0.insertSubMatrix(curpos, curpos + m, curpos + m, curpos + 2 * m, mmap.get(1));
            } else {
                // The n-th base arrival emits the summed arrival, carrying the
                // base marking, and restarts the chain at the first copy
                for (int k = 1; k <= K + 1; k++) {
                    d1[k - 1].insertSubMatrix(curpos, 0, curpos + m, m, mmap.get(k));
                }
            }
            curpos += m;
        }

        Matrix[] cells = new Matrix[K + 2];
        cells[0] = d0;
        for (int k = 0; k <= K; k++) {
            cells[k + 1] = d1[k];
        }
        return new MatrixCell(cells);
    }
}
