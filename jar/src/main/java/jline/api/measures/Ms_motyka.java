package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_motyka {
    private Ms_motyka() {}

    /**
     * Motyka distance between two probability distributions.
     * Part of the intersection family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Motyka distance
     */
    public static double ms_motyka(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sumMax = 0.0;
        double sumTotal = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            sumMax += Math.max(pi, qi);
            sumTotal += pi + qi;
        }

        return (sumTotal > 0) ? sumMax / sumTotal : 0.0;
    }
}
