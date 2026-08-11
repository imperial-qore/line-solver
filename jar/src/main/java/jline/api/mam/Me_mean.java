package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * ME mean computation algorithms.
 *
 * @since LINE 3.0
 */
public final class Me_mean {
    private Me_mean() {}

    /**
     * Computes the mean of a Matrix Exponential (ME) distribution.
     *
     * @param alpha The initial vector of the ME distribution
     * @param A The matrix parameter of the ME distribution
     * @return The mean of the ME distribution
     */
    public static double me_mean(Matrix alpha, Matrix A) {
        int n = A.getNumRows();
        Matrix e = Matrix.ones(n, 1);
        Matrix Ainv = A.inv();
        return -alpha.mult(Ainv).mult(e).get(0, 0);
    }

    /**
     * Computes the mean of a Matrix Exponential (ME) distribution using matrices stored in a MatrixCell.
     *
     * @param ME The Matrix Exponential distribution stored in a MatrixCell
     * @return The mean of the ME distribution
     */
    public static double me_mean(MatrixCell ME) {
        return Map_mean.map_mean(ME.get(0), ME.get(1));
    }
}
