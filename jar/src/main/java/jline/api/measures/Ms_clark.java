package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_clark {
    private Ms_clark() {}

    /**
     * Clark distance between two probability distributions.
     * Part of the squared L2 or chi-squared family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Clark distance
     */
    public static double ms_clark(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double denom = pi + qi;
            if (denom > 0) {
                double term = Math.abs(pi - qi) / denom;
                sum += term * term;
            }
        }

        return Math.sqrt(sum);
    }
}
