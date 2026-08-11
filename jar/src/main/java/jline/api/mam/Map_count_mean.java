/**
 * @file Markovian Arrival Process counting process mean analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.MatrixCell;

public final class Map_count_mean {
    private Map_count_mean() {}

    /**
     * Computes the mean of the counting process over multiple specified interval lengths for a given MAP.
     *
     * @param MAP The Markovian Arrival Process stored in a MatrixCell.
     * @param t   An array of interval lengths over which to compute the mean.
     * @return An array of doubles, where each element represents the mean of the counting process.
     */
    public static double[] map_count_mean(MatrixCell MAP, double[] t) {
        double[] ret = new double[t.length];
        double lambda = Map_lambda.map_lambda(MAP);
        for (int i = 0; i < t.length; i++) {
            ret[i] = lambda * t[i];
        }
        return ret;
    }

    /**
     * Computes the mean of the counting process over a specified interval length for a given MAP.
     *
     * @param MAP The Markovian Arrival Process stored in a MatrixCell.
     * @param t   The length of the interval over which to compute the mean.
     * @return The mean of the counting process over the interval `t`.
     */
    public static double map_count_mean(MatrixCell MAP, double t) {
        return Map_lambda.map_lambda(MAP) * t;
    }
}
