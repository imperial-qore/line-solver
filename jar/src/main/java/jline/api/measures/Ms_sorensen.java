package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_sorensen {
    private Ms_sorensen() {}

    /**
     * Sorensen distance between two probability distributions.
     * Part of the L1 family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Sorensen distance
     */
    public static double ms_sorensen(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sumDiff = 0.0;
        double sumTotal = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            sumDiff += Math.abs(pi - qi);
            sumTotal += pi + qi;
        }

        return (sumTotal > 0) ? sumDiff / sumTotal : 0.0;
    }
}
