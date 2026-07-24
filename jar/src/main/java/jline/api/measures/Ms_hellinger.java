/**
 * @file Hellinger distance for probability distributions
 *
 * Computes the Hellinger distance measuring dissimilarity between probability distributions
 * based on the square root of their densities. Bounded between 0 and sqrt(2), it's commonly
 * used in statistics for comparing continuous and discrete probability distributions.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_hellinger {
    private Ms_hellinger() {}

    /**
     * Hellinger distance between two probability distributions.
     * Part of the fidelity family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Hellinger distance
     */
    public static double ms_hellinger(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double diff = Math.sqrt(P.get(i)) - Math.sqrt(Q.get(i));
            sum += diff * diff;
        }

        return Math.sqrt(2.0 * sum);
    }
}
