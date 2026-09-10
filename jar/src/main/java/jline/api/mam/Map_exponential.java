/**
 * @file Markovian Arrival Process exponential distribution construction
 *
 * Creates MAP representations of exponential inter-arrival time distributions with specified
 * means. The simplest MAP form with single-state Poisson arrival processes.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_exponential {
    private Map_exponential() {}

    /**
     * Creates a Markovian Arrival Process (MAP) with an exponential inter-arrival time distribution.
     *
     * @param mean the desired mean of the exponential inter-arrival times
     * @return a MatrixCell containing the MAP transition matrices for the exponential distribution
     */
    public static MatrixCell map_exponential(double mean) {
        double mu = 1.0 / mean;
        MatrixCell MAP = new MatrixCell();
        Matrix D0 = new Matrix(1, 1, 1);
        D0.set(0, 0, -mu);
        Matrix D1 = new Matrix(1, 1, 1);
        D1.set(0, 0, mu);
        MAP.set(0, D0);
        MAP.set(1, D1);
        return MAP;
    }
}
