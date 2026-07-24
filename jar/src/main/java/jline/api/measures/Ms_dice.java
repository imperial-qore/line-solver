package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_dice {
    private Ms_dice() {}

    /**
     * Dice distance between two probability distributions.
     * Part of the inner product family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Dice distance
     */
    public static double ms_dice(Matrix P, Matrix Q) {
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

        double s = 2.0 * dotProduct / (sumP2 + sumQ2);
        return 1.0 - s;
    }
}
