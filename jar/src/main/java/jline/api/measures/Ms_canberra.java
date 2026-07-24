package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_canberra {
    private Ms_canberra() {}

    /**
     * Canberra distance between two probability distributions.
     * Part of the L1 family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Canberra distance
     */
    public static double ms_canberra(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double denom = pi + qi;
            if (denom > 0) {
                sum += Math.abs(pi - qi) / denom;
            }
        }

        return sum;
    }
}
