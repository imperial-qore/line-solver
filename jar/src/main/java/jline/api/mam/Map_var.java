package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

/**
 * MAP variance computation algorithms.
 */
public final class Map_var {
    private Map_var() {}

    /**
     * Computes the variance of the inter-arrival times for a MAP.
     *
     * @param D0 the hidden transition matrix of the MAP
     * @param D1 the visible transition matrix of the MAP
     * @return the variance of the inter-arrival times
     */
    public static double map_var(Matrix D0, Matrix D1) {
        return Map_moment.map_moment(D0, D1, 2) - FastMath.pow(Map_mean.map_mean(D0, D1), 2);
    }

    /**
     * Computes the variance of the inter-arrival times for a MAP using a MatrixCell.
     *
     * @param MAP the MatrixCell representing the MAP
     * @return the variance of the inter-arrival times
     */
    public static double map_var(MatrixCell MAP) {
        return map_var(MAP.get(0), MAP.get(1));
    }
}
