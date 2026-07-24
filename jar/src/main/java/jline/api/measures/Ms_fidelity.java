package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_fidelity {
    private Ms_fidelity() {}

    /**
     * Fidelity distance between two probability distributions.
     * Part of the fidelity family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Fidelity distance
     */
    public static double ms_fidelity(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            sum += Math.sqrt(P.get(i) * Q.get(i));
        }

        return 1.0 - sum;
    }
}
