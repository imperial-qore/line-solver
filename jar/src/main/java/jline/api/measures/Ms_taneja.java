package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_taneja {
    private Ms_taneja() {}

    /**
     * Taneja distance between two probability distributions.
     * Part of the combination family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Taneja distance
     */
    public static double ms_taneja(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double m = (pi + qi) / 2.0;
            double sqrtPQ = Math.sqrt(pi * qi);
            if (m > 0 && sqrtPQ > 0) {
                sum += m * Math.log(m / sqrtPQ);
            }
        }

        return sum;
    }
}
