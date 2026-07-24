package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_jeffreys {
    private Ms_jeffreys() {}

    /**
     * Jeffreys divergence between two probability distributions.
     * Part of Shannon's entropy family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Jeffreys divergence
     */
    public static double ms_jeffreys(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            if (pi > 0 && qi > 0) {
                sum += (pi - qi) * Math.log(pi / qi);
            }
        }

        return sum;
    }
}
