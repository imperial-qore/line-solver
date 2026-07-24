package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_squaredeuclidean {
    private Ms_squaredeuclidean() {}

    /**
     * Squared Euclidean distance between two probability distributions.
     * Part of the squared L2 or chi-squared family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Squared Euclidean distance
     */
    public static double ms_squaredeuclidean(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double diff = P.get(i) - Q.get(i);
            sum += diff * diff;
        }

        return sum;
    }
}
