/**
 * @file Marked Markovian Arrival Process exponential distribution construction
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_exponential {
    private Mmap_exponential() {}

    /**
     * Fits an order-n MMAP based on the given arrival rates for each job class.
     *
     * @param lambda a 1xK matrix describing the arrival rates for each job class
     * @param n      the number of states in the MMAP
     * @return a MatrixCell representing the MMAP
     */
    public static MatrixCell mmap_exponential(Matrix lambda, int n) {
        int K = lambda.length();
        MatrixCell MMAP = new MatrixCell(2 + K);
        MMAP.set(0, new Matrix(n, n, n ^ 2));
        MMAP.set(1, new Matrix(n, n, n ^ 2));
        for (int k = 0; k < K; k++) {
            double a = lambda.get(0, k);
            Matrix m = new Matrix(n, n, n ^ 2);
            for (int i = 0; i < n; i++) {
                m.set(i, n - 1 - i, a);
            }
            MMAP.set(2 + k, m);
            MMAP.set(1, MMAP.get(1).add(1.0, m));
        }
        return Mmap_normalize.mmap_normalize(MMAP);
    }

    /**
     * Fits a single-state MMAP based on the given arrival rates for each job class.
     *
     * @param lambda a 1xK matrix describing the arrival rates for each job class
     * @return a MatrixCell representing the MMAP
     */
    public static MatrixCell mmap_exponential(Matrix lambda) {
        return mmap_exponential(lambda, 1);
    }
}
