package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_kdivergence {
    private Ms_kdivergence() {}

    /**
     * K-divergence between two probability distributions.
     * Part of Shannon's entropy family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return K-divergence
     */
    public static double ms_kdivergence(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double denom = pi + qi;
            if (pi > 0 && denom > 0) {
                sum += pi * Math.log(2.0 * pi / denom);
            }
        }

        return sum;
    }
}
