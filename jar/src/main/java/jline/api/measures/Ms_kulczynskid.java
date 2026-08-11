package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_kulczynskid {
    private Ms_kulczynskid() {}

    /**
     * Kulczynski d distance between two probability distributions.
     * Part of the L1 family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Kulczynski d distance
     */
    public static double ms_kulczynskid(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sumDiff = 0.0;
        double sumMin = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            sumDiff += Math.abs(pi - qi);
            sumMin += Math.min(pi, qi);
        }

        return (sumMin > 0) ? sumDiff / sumMin : 0.0;
    }
}
