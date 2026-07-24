package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_lorentzian {
    private Ms_lorentzian() {}

    /**
     * Lorentzian distance between two probability distributions.
     * Part of the L1 family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Lorentzian distance
     */
    public static double ms_lorentzian(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            sum += Math.log(1.0 + Math.abs(P.get(i) - Q.get(i)));
        }

        return sum;
    }
}
