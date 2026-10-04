/**
 * @file Markov Modulated Poisson Process random generation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.Random;

public final class Mmpp_rand {
    private Mmpp_rand() {}

    /**
     * Generates a random Markov Modulated Poisson Process (MMPP) with K states.
     *
     * @param K the number of states
     * @return a MatrixCell representing the random MMPP transition matrices
     */
    public static MatrixCell mmpp_rand(int K) {
        Random random = new Random();

        // Generate random D0 matrix
        Matrix D0 = new Matrix(K, K);
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                D0.set(i, j, random.nextDouble());
            }
        }

        // Generate random D1 matrix (diagonal for MMPP)
        Matrix D1 = new Matrix(K, K);
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                if (i == j) {
                    D1.set(i, j, random.nextDouble());
                } else {
                    D1.set(i, j, 0.0);
                }
            }
        }

        return Map_normalize.map_normalize(D0, D1);
    }

    /**
     * Generates a random Markov Modulated Poisson Process (MMPP) with 2 states.
     *
     * @return a MatrixCell representing the random MMPP transition matrices
     */
    public static MatrixCell mmpp_rand() {
        return mmpp_rand(2);
    }
}
