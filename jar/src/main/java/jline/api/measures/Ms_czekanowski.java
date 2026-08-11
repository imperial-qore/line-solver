package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_czekanowski {
    private Ms_czekanowski() {}

    /**
     * Czekanowski distance between two probability distributions.
     * Part of the intersection family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Czekanowski distance
     */
    public static double ms_czekanowski(Matrix P, Matrix Q) {
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
