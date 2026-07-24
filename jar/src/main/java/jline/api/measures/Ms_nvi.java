/**
 * @file Normalized variation of information for discrete variables
 *
 * Computes normalized variation information (1-I(X,Y)/H(X,Y)) measuring the normalized
 * amount of information lost in going from the joint distribution to individual variables.
 * Used as a distance metric in clustering and classification validation.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import static jline.api.measures.Ms_entropy.ms_entropy;
import static jline.api.measures.Ms_jointentropy.ms_jointentropy;

import jline.util.matrix.Matrix;

public final class Ms_nvi {
    private Ms_nvi() {}

    /**
     * Compute normalized variation information z=(1-I(x,y)/H(x,y)) of two discrete variables x and y.
     *
     * @param x first matrix
     * @param y second matrix of the same length
     * @return normalized variation information z=(1-I(x,y)/H(x,y))
     */
    public static double ms_nvi(Matrix x, Matrix y) {
        if (x.length() != y.length()) {
            throw new IllegalArgumentException("Input matrices must have the same length");
        }

        if (x.isEmpty()) return 0.0;

        double hx = ms_entropy(x);
        double hy = ms_entropy(y);
        double hxy = ms_jointentropy(x, y);

        // Handle edge case where joint entropy is zero
        if (hxy == 0.0) return 0.0;

        // ms_nvi = 2 - (H(x) + H(y)) / H(x,y)
        double ms_nvi = 2.0 - (hx + hy) / hxy;
        return Math.max(0.0, ms_nvi);
    }
}
