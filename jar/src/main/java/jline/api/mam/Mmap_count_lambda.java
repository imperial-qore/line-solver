/**
 * @file Marked Markovian Arrival Process class-specific arrival rates
 *
 * Computes arrival rate vectors for each marked class in MMAP processes.
 * Fundamental for multiclass queueing analysis and capacity planning.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_count_lambda {
    private Mmap_count_lambda() {}

    /**
     * Computes the arrival rate vector of the counting process for the given Marked MAP (MMAP).
     *
     * This method calculates the arrival rates for each job class in the MMAP. The arrival rate for each class is determined
     * by multiplying the stationary distribution of the underlying Markov chain (computed by map_piq) with the corresponding
     * marking matrix and summing the elements. The result is a vector of arrival rates, one for each job class.
     *
     * @param mmap the Marked MAP represented by a MatrixCell
     * @return a matrix (vector) with the arrival rate for each job class
     */
    public static Matrix mmap_count_lambda(MatrixCell mmap) {
        int K = mmap.size() - 2;

        // symbolic map haven't been implemented
        Matrix lk = new Matrix(K, 1, K);

        Matrix theta = Map_piq.map_piq(mmap.get(0), mmap.get(1));

        for (int k = 0; k < K; k++) {
            lk.set(k, 0, theta.mult(mmap.get(2 + k)).elementSum());
        }

        return lk.transpose();
    }
}
