package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_tanimoto {
    private Ms_tanimoto() {}

    /**
     * Tanimoto distance between two probability distributions.
     * Part of the intersection family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Tanimoto distance
     */
    public static double ms_tanimoto(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sumDiff = 0.0;
        double sumMax = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double maxVal = Math.max(pi, qi);
            double minVal = Math.min(pi, qi);
            sumDiff += maxVal - minVal;
            sumMax += maxVal;
        }

        return sumMax > 0 ? sumDiff / sumMax : 0.0;
    }
}
