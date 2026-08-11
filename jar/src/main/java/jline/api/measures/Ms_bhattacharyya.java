/**
 * @file Bhattacharyya distance for probability distributions
 *
 * Computes the Bhattacharyya distance measuring the similarity between probability
 * distributions based on the Bhattacharyya coefficient. Widely used in pattern
 * recognition, classification, and feature selection for comparing statistical models.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_bhattacharyya {
    private Ms_bhattacharyya() {}

    /**
     * Bhattacharyya distance between two probability distributions.
     * Part of the fidelity family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Bhattacharyya distance
     */
    public static double ms_bhattacharyya(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            sum += Math.sqrt(P.get(i) * Q.get(i));
        }

        return -Math.log(sum);
    }
}
