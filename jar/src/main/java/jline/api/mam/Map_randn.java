/**
 * @file Markovian Arrival Process random generation with noise
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public final class Map_randn {
    private Map_randn() {}

    /**
     * Generates a random MAP with 2 states using N(1, 2^2) magnitudes.
     *
     * @return a MatrixCell representing the random MAP transition matrices
     */
    public static MatrixCell map_randn() {
        return map_randn(2, 1.0, 2.0);
    }

    /**
     * Generates a random MAP with K states using normal distribution.
     *
     * @param K     the number of states
     * @param mu    the mean of the normal distribution
     * @param sigma the standard deviation of the normal distribution
     * @return a MatrixCell representing the random MAP transition matrices
     */
    public static MatrixCell map_randn(int K, double mu, double sigma) {
        Random random = new Random();
        // Randomize D0
        List<List<Double>> D0 = new ArrayList<List<Double>>();
        for (int i = 0; i < K; i++) {
            List<Double> row = new ArrayList<Double>();
            for (int j = 0; j < K; j++) {
                row.add(Math.abs(mu + sigma * random.nextGaussian()));
            }
            D0.add(row);
        }
        // Randomize D1
        List<List<Double>> D1 = new ArrayList<List<Double>>();
        for (int i = 0; i < K; i++) {
            List<Double> row = new ArrayList<Double>();
            for (int j = 0; j < K; j++) {
                row.add(Math.abs(mu + sigma * random.nextGaussian()));
            }
            D1.add(row);
        }
        return Map_normalize.map_normalize(Matrix.fromRows(D0), Matrix.fromRows(D1));
    }
}
