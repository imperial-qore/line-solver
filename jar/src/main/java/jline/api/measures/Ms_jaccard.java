package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_jaccard {
    private Ms_jaccard() {}

    /**
     * Jaccard distance between two probability distributions.
     * Part of the inner product family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Jaccard distance
     */
    public static double ms_jaccard(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double dotProduct = 0.0;
        double sumP2 = 0.0;
        double sumQ2 = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            dotProduct += pi * qi;
            sumP2 += pi * pi;
            sumQ2 += qi * qi;
        }

        double s = dotProduct / (sumP2 + sumQ2 - dotProduct);
        return 1.0 - s;
    }
}
