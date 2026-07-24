package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_neymanchisquared {
    private Ms_neymanchisquared() {}

    /**
     * Neyman chi-squared distance between two probability distributions.
     * Part of the squared L2 or chi-squared family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Neyman chi-squared distance
     */
    public static double ms_neymanchisquared(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            if (pi > 0) {
                double diff = pi - qi;
                sum += (diff * diff) / pi;
            }
        }

        return sum;
    }
}
