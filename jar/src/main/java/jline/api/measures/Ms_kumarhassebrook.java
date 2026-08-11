package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_kumarhassebrook {
    private Ms_kumarhassebrook() {}

    /**
     * Kumar-Hassebrook distance between two probability distributions.
     * Part of the inner product family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Kumar-Hassebrook distance
     */
    public static double ms_kumarhassebrook(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double dotProduct = 0.0;
        double sumDiff2 = 0.0;

        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            dotProduct += pi * qi;
            double diff = Math.abs(pi - qi);
            sumDiff2 += diff * diff;
        }

        double s = dotProduct / (sumDiff2 + dotProduct);
        return 1.0 - s;
    }
}
