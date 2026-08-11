package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_kumarjohnson {
    private Ms_kumarjohnson() {}

    /**
     * Kumar-Johnson distance between two probability distributions.
     * Part of the combination family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Kumar-Johnson distance
     */
    public static double ms_kumarjohnson(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double denom = 2.0 * Math.pow(pi * qi, 3.0 / 2.0);
            if (denom > 0) {
                double diff = pi * pi - qi * qi;
                sum += (diff * diff) / denom;
            }
        }

        return sum;
    }
}
