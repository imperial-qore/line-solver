package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_ruzicka {
    private Ms_ruzicka() {}

    /**
     * Ruzicka distance between two probability distributions.
     * Part of the intersection family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Ruzicka distance
     */
    public static double ms_ruzicka(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sumMin = 0.0;
        double sumMax = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            sumMin += Math.min(pi, qi);
            sumMax += Math.max(pi, qi);
        }

        return (sumMax > 0) ? 1.0 - sumMin / sumMax : 0.0;
    }
}
