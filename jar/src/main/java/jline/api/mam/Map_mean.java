package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * MAP mean inter-arrival time computation algorithms.
 *
 * @since LINE 3.0
 */
public final class Map_mean {
    private Map_mean() {}

    /**
     * Computes the mean inter-arrival time of a Markovian Arrival Process (MAP).
     *
     * @param D0 The hidden transition matrix of the MAP.
     * @param D1 The visible transition matrix of the MAP.
     * @return The mean inter-arrival time of the MAP.
     */
    public static double map_mean(Matrix D0, Matrix D1) {
        return 1.0 / Map_lambda.map_lambda(D0, D1);
    }

    /**
     * Computes the mean inter-arrival time of a MAP using matrices stored in a MatrixCell.
     *
     * @param MAP The Markovian Arrival Process stored in a MatrixCell.
     * @return The mean inter-arrival time of the MAP.
     */
    public static double map_mean(MatrixCell MAP) {
        return map_mean(MAP.get(0), MAP.get(1));
    }
}
