/**
 * Minkowski Distance for Probability Distributions
 *
 * Implements the generalized Minkowski distance (Lp norm) between probability distributions.
 * Special cases include Manhattan distance (p=1) and Euclidean distance (p=2). Widely
 * used in machine learning for measuring dissimilarity between feature vectors.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_minkowski {
    private Ms_minkowski() {}

    /**
     * Minkowski distance between two probability distributions.
     * Generalization of Euclidean (p=2) and Manhattan (p=1) distances.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @param p order parameter (p &gt;= 1)
     * @return Minkowski distance of order p
     */
    public static double ms_minkowski(Matrix P, Matrix Q, double p) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }
        if (p < 1.0) {
            throw new IllegalArgumentException("Order parameter p must be >= 1");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            sum += Math.pow(Math.abs(P.get(i) - Q.get(i)), p);
        }

        return Math.pow(sum, 1.0 / p);
    }
}
