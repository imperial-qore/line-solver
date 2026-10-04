/**
 * @file Markovian Arrival Process random generation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public final class Map_rand {
    private Map_rand() {}

    /**
     * Generates a random Markovian Arrival Process (MAP) with 2 states.
     *
     * @return a MatrixCell representing the random MAP transition matrices
     */
    public static MatrixCell map_rand() {
        return map_rand(2);
    }

    /**
     * Generates a random Markovian Arrival Process (MAP) with K states.
     *
     * @param K the number of states
     * @return a MatrixCell representing the random MAP transition matrices
     */
    public static MatrixCell map_rand(int K) {
        Random random = new Random();
        // Randomize D0
        List<List<Double>> D0 = new ArrayList<List<Double>>();
        for (int i = 0; i < K; i++) {
            List<Double> row = new ArrayList<Double>();
            for (int j = 0; j < K; j++) {
                row.add(random.nextDouble());
            }
            D0.add(row);
        }
        // Randomize D1
        List<List<Double>> D1 = new ArrayList<List<Double>>();
        for (int i = 0; i < K; i++) {
            List<Double> row = new ArrayList<Double>();
            for (int j = 0; j < K; j++) {
                row.add(random.nextDouble());
            }
            D1.add(row);
        }
        return Map_normalize.map_normalize(Matrix.fromRows(D0), Matrix.fromRows(D1));
    }
}
