/**
 * @file Markovian Arrival Process Erlang-k distribution fitting
 *
 * Constructs MAP representations of Erlang-k processes with specified means and phases.
 * Used for modeling low-variability arrival processes with coefficient of variation less than 1.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_erlang {
    private Map_erlang() {}

    /**
     * Fits an Erlang-k process as a Markovian Arrival Process (MAP).
     *
     * The Erlang-k process is characterized by k phases, each with an exponential distribution, resulting in a distribution
     * with a mean and a coefficient of variation smaller than 1. This method constructs a MAP that approximates the behavior
     * of an Erlang-k process with a specified mean and shape parameter k.
     *
     * @param mean the desired mean of the Erlang-k process
     * @param k    the shape parameter, representing the number of phases in the Erlang-k process
     * @return a MatrixCell containing the transition matrices D0 and D1 of the fitted MAP
     */
    public static MatrixCell map_erlang(double mean, int k) {
        double mu = (double) k / mean;
        MatrixCell MAP = new MatrixCell();
        Matrix D0 = new Matrix(k, k, 2 * k - 1);
        D0.set(0, 0, -mu);
        for (int i = 0; i < k - 1; i++) {
            D0.set(i, i + 1, mu);
        }
        Matrix D1 = new Matrix(k, k, 1);
        D1.set(k - 1, 0, mu);
        MAP.set(0, D0);
        MAP.set(1, D1);
        MAP = Map_normalize.map_normalize(MAP);
        return MAP;
    }
}
