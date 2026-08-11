package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_divergence {
    private Ms_divergence() {}

    /**
     * Divergence distance between two probability distributions.
     * Part of the squared L2 or chi-squared family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Divergence distance
     */
    public static double ms_divergence(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double sumPQ = pi + qi;
            if (sumPQ > 0) {
                double diff = pi - qi;
                sum += (diff * diff) / (sumPQ * sumPQ);
            }
        }

        return 2.0 * sum;
    }
}
