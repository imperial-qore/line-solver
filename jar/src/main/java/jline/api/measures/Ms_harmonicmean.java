package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_harmonicmean {
    private Ms_harmonicmean() {}

    /**
     * Harmonic mean distance between two probability distributions.
     * Part of the inner product family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Harmonic mean distance
     */
    public static double ms_harmonicmean(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double denom = pi + qi;
            if (denom > 0) {
                sum += pi * qi / denom;
            }
        }

        return 1.0 - 2.0 * sum;
    }
}
