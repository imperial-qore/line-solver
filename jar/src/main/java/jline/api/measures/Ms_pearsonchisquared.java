package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_pearsonchisquared {
    private Ms_pearsonchisquared() {}

    /**
     * Pearson chi-squared distance between two probability distributions.
     * Part of the squared L2 or chi-squared family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Pearson chi-squared distance
     */
    public static double ms_pearsonchisquared(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            if (qi > 0) {
                double diff = pi - qi;
                sum += (diff * diff) / qi;
            }
        }

        return sum;
    }
}
