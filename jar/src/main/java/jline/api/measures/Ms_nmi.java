/**
 * @file Normalized mutual information for discrete variables
 *
 * Computes normalized mutual information I(X,Y)/sqrt(H(X)*H(Y)) providing a scale-invariant
 * measure of dependence between variables. Commonly used in clustering evaluation and
 * feature selection where normalized measures are preferred over raw mutual information.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.util.matrix.Matrix;

import static jline.api.measures.Ms_entropy.ms_entropy;
import static jline.api.measures.Ms_mutinfo.ms_mutinfo;

public final class Ms_nmi {
    private Ms_nmi() {}

    /**
     * Compute normalized mutual information I(x,y)/sqrt(H(x)*H(y)) of two discrete variables x and y.
     *
     * @param x first matrix
     * @param y second matrix of the same length
     * @return normalized mutual information z=I(x,y)/sqrt(H(x)*H(y))
     */
    public static double ms_nmi(Matrix x, Matrix y) {
        if (x.length() != y.length()) {
            throw new IllegalArgumentException("Input matrices must have the same length");
        }

        if (x.isEmpty()) return 0.0;

        double hx = ms_entropy(x);
        double hy = ms_entropy(y);

        // Handle edge cases where entropy is zero
        if (hx == 0.0 || hy == 0.0) return 0.0;

        double mi = ms_mutinfo(x, y);

        // Normalized mutual information
        double nmi = Math.sqrt((mi / hx) * (mi / hy));
        return Math.max(0.0, nmi);
    }
}
