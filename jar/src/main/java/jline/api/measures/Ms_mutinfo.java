package jline.api.measures;

import static jline.api.measures.Ms_entropy.ms_entropy;
import static jline.api.measures.Ms_jointentropy.ms_jointentropy;

import jline.util.matrix.Matrix;

/**
 * Mutual information for discrete variables.
 *
 * <p>Computes mutual information I(X,Y) = H(X) + H(Y) - H(X,Y) measuring the amount of
 * information obtained about one variable through observing the other. Fundamental
 * measure in feature selection, clustering validation, and dependency analysis.
 *
 * @since LINE 3.0
 */
public final class Ms_mutinfo {
    private Ms_mutinfo() {}

    /**
     * Compute mutual information I(x,y) of two discrete variables x and y.
     *
     * @param x first matrix
     * @param y second matrix of the same length
     * @return mutual information z=I(x,y)
     */
    public static double ms_mutinfo(Matrix x, Matrix y) {
        if (x.length() != y.length()) {
            throw new IllegalArgumentException("Input matrices must have the same length");
        }

        if (x.isEmpty()) return 0.0;

        // I(x,y) = H(x) + H(y) - H(x,y)
        double hx = ms_entropy(x);
        double hy = ms_entropy(y);
        double hxy = ms_jointentropy(x, y);

        double mutualInfo = hx + hy - hxy;
        return Math.max(0.0, mutualInfo);
    }
}
