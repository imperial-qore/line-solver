package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_squaredchord {
    private Ms_squaredchord() {}

    /**
     * Squared chord distance between two probability distributions.
     * Part of the fidelity family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Squared chord distance
     */
    public static double ms_squaredchord(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double diff = Math.sqrt(P.get(i)) - Math.sqrt(Q.get(i));
            sum += diff * diff;
        }

        return sum;
    }
}
