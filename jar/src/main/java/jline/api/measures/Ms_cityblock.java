/**
 * @file City block (Manhattan) distance metric
 *
 * Implements the City block distance d(p,q) = sum|p_i - q_i| between probability
 * distributions. Also known as Manhattan distance or L1 norm, commonly used
 * in clustering and nearest neighbor applications for robust distance measurement.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_cityblock {
    private Ms_cityblock() {}

    /**
     * City block (Manhattan) distance between two probability distributions.
     * Part of the Minkowski family. Also known as L1 distance.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return City block distance
     */
    public static double ms_cityblock(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double sum = 0.0;
        for (int i = 0; i < P.length(); i++) {
            sum += Math.abs(P.get(i) - Q.get(i));
        }

        return sum;
    }
}
