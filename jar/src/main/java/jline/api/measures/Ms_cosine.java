/**
 * @file Cosine distance for probability distributions
 *
 * Computes cosine distance (1 - cosine similarity) measuring the angle between two
 * vectors in high-dimensional space. Particularly useful for text analysis and
 * information retrieval where magnitude is less important than direction.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_cosine {
    private Ms_cosine() {}

    /**
     * Cosine distance between two probability distributions.
     * Part of the inner product family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Cosine distance
     */
    public static double ms_cosine(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double dotProduct = 0.0;
        double sumP2 = 0.0;
        double sumQ2 = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            dotProduct += pi * qi;
            sumP2 += pi * pi;
            sumQ2 += qi * qi;
        }

        double s = dotProduct / Math.sqrt(sumP2) / Math.sqrt(sumQ2);
        return 1.0 - s;
    }
}
