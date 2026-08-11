package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_kulczynskis {
    private Ms_kulczynskis() {}

    /**
     * Kulczynski s distance between two probability distributions.
     * Part of the intersection family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Kulczynski s distance
     */
    public static double ms_kulczynskis(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sumMin = 0.0;
        double sumDiff = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            sumMin += Math.min(pi, qi);
            sumDiff += Math.abs(pi - qi);
        }

        return (sumMin > 0) ? sumDiff / sumMin : 0.0;
    }
}
