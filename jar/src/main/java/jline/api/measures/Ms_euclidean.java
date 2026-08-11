/**
 * @file Euclidean distance metric
 *
 * Implements the standard Euclidean distance d(p,q) = sqrt(sum((pi-qi)^2)) between probability
 * distributions or vectors. The most commonly used distance metric providing geometric
 * interpretation of dissimilarity in statistical analysis and machine learning.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_euclidean {
    private Ms_euclidean() {}

    /**
     * Euclidean distance between two probability distributions.
     * Part of the Minkowski family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Euclidean distance
     */
    public static double ms_euclidean(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double diff = P.get(i) - Q.get(i);
            sum += diff * diff;
        }

        return Math.sqrt(sum);
    }
}
