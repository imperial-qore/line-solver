package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_soergel {
    private Ms_soergel() {}

    /**
     * Soergel distance between two probability distributions.
     * Part of the L1 family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Soergel distance
     */
    public static double ms_soergel(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sumDiff = 0.0;
        double sumMax = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            sumDiff += Math.abs(pi - qi);
            sumMax += Math.max(pi, qi);
        }

        return (sumMax > 0) ? sumDiff / sumMax : 0.0;
    }
}
