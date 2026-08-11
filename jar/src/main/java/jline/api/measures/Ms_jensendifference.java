package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_jensendifference {
    private Ms_jensendifference() {}

    /**
     * Jensen difference divergence between two probability distributions.
     * Part of Shannon's entropy family.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Jensen difference
     */
    public static double ms_jensendifference(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            double pi = P.get(i);
            double qi = Q.get(i);
            double m = (pi + qi) / 2.0;

            double term = 0.0;
            if (pi > 0) term += pi * Math.log(pi);
            if (qi > 0) term += qi * Math.log(qi);
            term /= 2.0;

            if (m > 0) term -= m * Math.log(m);

            sum += term;
        }

        return sum;
    }
}
