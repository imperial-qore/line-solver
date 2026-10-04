package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * MMAP arrival rate computation algorithms.
 *
 * Provides methods for computing arrival rates of Marked Markovian Arrival Processes (MMAP).
 * MMAPs extend MAPs by allowing multiple types of arrivals, each with its own marking,
 * making them essential for multi-class queueing system modeling.
 *
 * @since LINE 3.0
 */
public final class Mmap_lambda {
    private Mmap_lambda() {}

    /**
     * Alias for mmap_count_lambda.
     *
     * This method provides an alternative name for `mmap_count_lambda`, offering the same functionality.
     *
     * @param MMAP the Marked MAP represented by a MatrixCell
     * @return a matrix (vector) with the arrival rate for each job class
     */
    public static Matrix mmap_lambda(MatrixCell MMAP) {
        return Mmap_count_lambda.mmap_count_lambda(MMAP);
    }
}
