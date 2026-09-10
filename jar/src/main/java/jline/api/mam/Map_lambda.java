package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * MAP arrival rate computation algorithms.
 */
public final class Map_lambda {
    private Map_lambda() {}

    /**
     * Computes the arrival rate (lambda) of a Markovian Arrival Process (MAP).
     *
     * @param D0 The hidden transition matrix of the MAP.
     * @param D1 The visible transition matrix of the MAP.
     * @return The arrival rate (lambda) of the MAP.
     */
    public static double map_lambda(Matrix D0, Matrix D1) {
        Matrix e = Matrix.ones(D0.getNumRows(), 1);
        Matrix lambda = Map_piq.map_piq(D0, D1); // piq
        lambda = lambda.mult(D1);
        lambda = lambda.mult(e);
        return lambda.toDouble();
    }

    /**
     * Computes the arrival rate (lambda) of a MAP using matrices stored in a MatrixCell.
     *
     * @param MAP The Markovian Arrival Process stored in a MatrixCell.
     * @return The arrival rate (lambda) of the MAP.
     */
    public static double map_lambda(MatrixCell MAP) {
        return map_lambda(MAP.get(0), MAP.get(1));
    }
}
