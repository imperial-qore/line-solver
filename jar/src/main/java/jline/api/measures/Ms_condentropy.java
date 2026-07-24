package jline.api.measures;

import static jline.api.measures.Ms_entropy.ms_entropy;
import static jline.api.measures.Ms_jointentropy.ms_jointentropy;

import jline.util.matrix.Matrix;

public final class Ms_condentropy {
    private Ms_condentropy() {}

    /**
     * Compute conditional entropy z=H(x|y) of two discrete variables x and y.
     *
     * @param x first matrix
     * @param y second matrix of the same length
     * @return conditional entropy z=H(x|y)
     */
    public static double ms_condentropy(Matrix x, Matrix y) {
        if (x.length() != y.length()) {
            throw new IllegalArgumentException("Input matrices must have the same length");
        }

        if (x.isEmpty()) {
            return 0.0;
        }

        // H(x|y) = H(x,y) - H(y)
        double hxy = ms_jointentropy(x, y);
        double hy = ms_entropy(y);

        double condEntropy = hxy - hy;
        return Math.max(0.0, condEntropy);
    }
}
