/**
 * @file Markovian Arrival Process squared coefficient of variation analysis
 *
 * Computes SCV of MAP inter-arrival times as normalized dispersion measure.
 * Fundamental metric for characterizing variability and burstiness in arrival processes.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_scv {
    private Map_scv() {}

    /**
     * Computes the squared coefficient of variation (SCV) of the inter-arrival times of a Markovian Arrival Process (MAP).
     *
     * The SCV is a normalized measure of the dispersion of the inter-arrival time distribution. It is calculated as the
     * variance of the inter-arrival times divided by the square of the mean inter-arrival time. The MAP is represented by
     * two matrices: D0 and D1, where D0 is the hidden transition matrix and D1 is the visible transition matrix.
     *
     * @param D0 the hidden transition matrix of the MAP
     * @param D1 the visible transition matrix of the MAP
     * @return the squared coefficient of variation (SCV) of the inter-arrival times
     */
    public static double map_scv(Matrix D0, Matrix D1) {
        double e1 = Map_moment.map_moment(D0, D1, 1);
        double e2 = Map_moment.map_moment(D0, D1, 2);

        double var = e2 - e1 * e1;
        double scv = var / e1 / e1;
        return scv;
    }

    /**
     * Computes the squared coefficient of variation (SCV) of the inter-arrival times of a MAP
     * stored in a MatrixCell that contains the MAP's transition matrices.
     *
     * @param MAP a MatrixCell containing the transition matrices D0 and D1 of the MAP
     * @return the squared coefficient of variation (SCV) of the inter-arrival times
     */
    public static double map_scv(MatrixCell MAP) {
        return map_scv(MAP.get(0), MAP.get(1));
    }
}
