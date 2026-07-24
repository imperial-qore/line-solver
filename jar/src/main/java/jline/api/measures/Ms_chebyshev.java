/**
 * Chebyshev Distance for Probability Distributions
 *
 * Implements the Chebyshev distance (L-infinity norm) computing the maximum absolute difference
 * between corresponding elements. Also known as chess-board distance, it represents
 * the limiting case of Minkowski distance as p approaches infinity.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.util.matrix.Matrix;

public final class Ms_chebyshev {
    private Ms_chebyshev() {}

    /**
     * Chebyshev distance between two probability distributions.
     * Part of the Minkowski family. Also known as L-infinity distance.
     *
     * @param P first probability distribution
     * @param Q second probability distribution
     * @return Chebyshev distance
     */
    public static double ms_chebyshev(Matrix P, Matrix Q) {
        if (P.length() != Q.length()) {
            throw new IllegalArgumentException("Distributions must have the same size");
        }

        double maxDist = 0.0;
        for (int i = 0; i < P.length(); i++) {
            maxDist = Math.max(maxDist, Math.abs(P.get(i) - Q.get(i)));
        }

        return maxDist;
    }
}
