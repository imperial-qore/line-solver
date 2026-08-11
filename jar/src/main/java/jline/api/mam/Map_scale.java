/**
 * @file Markovian Arrival Process temporal scaling operations
 *
 * Rescales MAP inter-arrival time distributions to achieve specified mean values.
 * Essential for parameter adjustment and model calibration in queueing analysis.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_scale {
    private Map_scale() {}

    /**
     * Rescales the mean inter-arrival time of a Markovian Arrival Process (MAP) to a specified new mean.
     *
     * <p>This method adjusts the transition matrices D0 and D1 of the MAP to achieve a desired mean
     * inter-arrival time while preserving the structure and relative rates of transitions. The MAP is
     * represented by two matrices: D0 and D1, where D0 is the hidden transition matrix and D1 is the
     * visible transition matrix.
     *
     * @param D0      the hidden transition matrix of the MAP
     * @param D1      the visible transition matrix of the MAP
     * @param newMean the desired new mean inter-arrival time
     * @return a MatrixCell containing the scaled MAP transition matrices
     */
    public static MatrixCell map_scale(Matrix D0, Matrix D1, double newMean) {
        double ratio = Map_mean.map_mean(D0, D1) / newMean;
        Matrix D0s = D0.copy();
        Matrix D1s = D1.copy();
        D0s.scaleEq(ratio);
        D1s.scaleEq(ratio);
        return Map_normalize.map_normalize(D0s, D1s);
    }

    /**
     * Rescales the mean inter-arrival time of a MAP stored in a MatrixCell that contains the MAP's
     * transition matrices.
     *
     * @param MAP     a MatrixCell containing the transition matrices D0 and D1 of the MAP
     * @param newMean the desired new mean inter-arrival time
     * @return a MatrixCell containing the scaled MAP transition matrices
     */
    public static MatrixCell map_scale(MatrixCell MAP, double newMean) {
        return map_scale(MAP.get(0), MAP.get(1), newMean);
    }
}
