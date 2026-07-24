package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_squaredchisquared {
    private Ms_squaredchisquared() {}

    /**
     * Squared chi-squared distance between two probability distributions.
     * Part of the squared L2 or chi-squared family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Squared chi-squared distance
     */
    public static double ms_squaredchisquared(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double denom = pi + qi;
            if (denom > 0) {
                double diff = pi - qi;
                sum += (diff * diff) / denom;
            }
        }

        return sum;
    }
}
