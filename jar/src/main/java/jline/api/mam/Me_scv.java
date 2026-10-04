package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * ME squared coefficient of variation (SCV) computation algorithms.
 *
 * Provides methods for computing the squared coefficient of variation of Matrix Exponential (ME) distributions.
 * The SCV is the variance divided by the square of the mean: SCV = Var(X) / E[X]^2
 *
 * @since LINE 3.0
 */
public final class Me_scv {
    private Me_scv() {}

    /**
     * Computes the squared coefficient of variation (SCV) of a Matrix Exponential (ME) distribution.
     *
     * The SCV is calculated as: SCV = variance / mean^2
     *
     * @param alpha The initial vector of the ME distribution
     * @param A The matrix parameter of the ME distribution
     * @return The squared coefficient of variation of the ME distribution
     */
    public static double me_scv(Matrix alpha, Matrix A) {
        double mean = Me_mean.me_mean(alpha, A);
        double variance = Me_var.me_var(alpha, A);
        return variance / (mean * mean);
    }

    /**
     * Computes the squared coefficient of variation (SCV) of an ME distribution using matrices stored in a MatrixCell.
     *
     * This is a convenience method that computes the SCV using the process representation.
     *
     * @param ME The Matrix Exponential distribution stored in a MatrixCell
     * @return The squared coefficient of variation of the ME distribution
     */
    public static double me_scv(MatrixCell ME) {
        // Use map_scv on the process representation
        return Map_scv.map_scv(ME.get(0), ME.get(1));
    }
}
