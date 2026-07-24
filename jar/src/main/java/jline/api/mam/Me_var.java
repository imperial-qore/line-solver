package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * ME variance computation algorithms.
 *
 * @since LINE 3.0
 */
public final class Me_var {
    private Me_var() {}

    /**
     * Computes the variance of a Matrix Exponential (ME) distribution.
     *
     * @param alpha The initial vector of the ME distribution
     * @param A The matrix parameter of the ME distribution
     * @return The variance of the ME distribution
     */
    public static double me_var(Matrix alpha, Matrix A) {
        Matrix e = Matrix.ones(A.getNumRows(), 1);
        Matrix Ainv = A.inv();
        Matrix Ainv2 = Ainv.mult(Ainv);

        double m1 = -alpha.mult(Ainv).mult(e).get(0, 0);
        double m2 = 2.0 * alpha.mult(Ainv2).mult(e).get(0, 0);

        return m2 - m1 * m1;
    }

    /**
     * Computes the variance of a Matrix Exponential (ME) distribution using a MatrixCell.
     *
     * @param ME The Matrix Exponential distribution stored in a MatrixCell
     * @return The variance of the ME distribution
     */
    public static double me_var(MatrixCell ME) {
        return Map_var.map_var(ME.get(0), ME.get(1));
    }
}
