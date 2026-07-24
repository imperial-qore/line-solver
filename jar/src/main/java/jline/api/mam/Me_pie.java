/**
 * @file ME initial probability computation algorithms.
 *
 * Provides methods for computing the stationary initial probability vector of Matrix Exponential (ME)
 * distributions, primarily for use with RAP (Rational Arrival Process) distributions.
 *
 * <p>For RAP distributions represented as (H0, H1), the initial probability is the stationary
 * distribution of the generator H0 + H1.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Me_pie {
    private Me_pie() {}

    /**
     * Computes the stationary initial probability for an ME/RAP distribution.
     *
     * <p>For RAP distributions with matrices (H0, H1), this computes the stationary distribution
     * of the generator Q = H0 + H1, which satisfies: pi * Q = 0 and pi * e = 1
     *
     * <p>This is equivalent to map_pie for the process representation.
     *
     * @param H0 The H0 matrix (hidden transitions)
     * @param H1 The H1 matrix (visible transitions)
     * @return The stationary initial probability vector
     */
    public static Matrix me_pie(Matrix H0, Matrix H1) {
        // Delegate to map_pie which computes the stationary distribution
        // of the generator D0 + D1
        return Map_pie.map_pie(H0, H1);
    }

    /**
     * Computes the stationary initial probability for an ME/RAP distribution using matrices stored in
     * a MatrixCell.
     *
     * @param ME The Matrix Exponential/RAP distribution stored in a MatrixCell
     * @return The stationary initial probability vector
     */
    public static Matrix me_pie(MatrixCell ME) {
        return Map_pie.map_pie(ME.get(0), ME.get(1));
    }
}
