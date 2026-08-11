package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_gower {
    private Ms_gower() {}

    /**
     * Gower distance between two probability distributions.
     * Part of the L1 family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Gower distance
     */
    public static double ms_gower(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            sum += Math.abs(P.get(i) - Q.get(i));
        }

        return sum / P.length();
    }
}
