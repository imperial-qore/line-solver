/**
 * @file Marked Markovian Arrival Process counting process mean analysis
 *
 * Computes mean count vectors for each marked class in MMAP counting processes.
 * Essential for determining expected arrival rates per class in multiclass systems.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_count_mean {
    private Mmap_count_mean() {}

    /**
     * Computes the mean count vector of events of different types in a Markovian Arrival
     * Process with marked arrivals (MMAP) over a time period.
     *
     * @param MMAP the MatrixCell containing the transition matrices of the MMAP, with
     *             D0, D1, ..., Dc representing different types of events
     * @param t    the time period over which to compute the mean counts
     * @return a Matrix containing the mean count vector of events over time t
     */
    public static Matrix mmap_count_mean(MatrixCell MMAP, double t) {
        Matrix theta = Map_piq.map_piq(MMAP.get(0), MMAP.get(1));
        int C = MMAP.size() - 2;
        Matrix et = new Matrix(1, C);
        for (int c = 0; c < C; c++) {
            et.set(c, theta.mult(MMAP.get(2 + c)).elementSum() * t);
        }
        return et;
    }
}
