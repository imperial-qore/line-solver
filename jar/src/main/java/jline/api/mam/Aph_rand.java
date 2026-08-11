package jline.api.mam;

import java.util.Random;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Aph_rand {
    private Aph_rand() {}

    /**
     * Generates a random Acyclic Phase-type (APH) distribution with K phases.
     *
     * @param K the number of phases
     * @return a MatrixCell representing the random APH distribution
     */
    public static MatrixCell aph_rand(int K) {
        Random random = new Random();

        // Generate random D1 matrix (full matrix)
        Matrix D1 = new Matrix(K, K);
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                D1.set(i, j, random.nextDouble());
            }
        }

        // Generate random D0 matrix with acyclic structure
        Matrix D0 = new Matrix(K, K);
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                if (j < i) {
                    D0.set(i, j, 0.0);
                } else {
                    D0.set(i, j, random.nextDouble());
                }
            }
        }

        // Create temporary MAP
        MatrixCell tempMAP = new MatrixCell();
        tempMAP.set(0, D0);
        tempMAP.set(1, D1);

        // Apply renewal process and normalize
        MatrixCell renewalMAP = Map_renewal.map_renewal(tempMAP);
        return Map_normalize.map_normalize(renewalMAP);
    }

    /**
     * Generates a random Acyclic Phase-type (APH) distribution with 2 phases.
     */
    public static MatrixCell aph_rand() {
        return aph_rand(2);
    }
}
