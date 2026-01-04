package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * ME initial probability computation algorithms.
 *
 * Provides methods for computing the stationary initial probability vector of Matrix Exponential (ME)
 * distributions, primarily for use with RAP (Rational Arrival Process) distributions.
 *
 * For RAP distributions represented as (H0, H1), the initial probability is the stationary
 * distribution of the generator H0 + H1.
 *
 * @since LINE 3.0
 */

/**
 * Computes the stationary initial probability for an ME/RAP distribution.
 *
 * For RAP distributions with matrices (H0, H1), this computes the stationary distribution
 * of the generator Q = H0 + H1, which satisfies: pi * Q = 0 and pi * e = 1
 *
 * This is equivalent to map_pie for the process representation.
 *
 * @param H0 The H0 matrix (hidden transitions)
 * @param H1 The H1 matrix (visible transitions)
 * @return The stationary initial probability vector
 */
fun me_pie(H0: Matrix, H1: Matrix): Matrix {
    // Delegate to map_pie which computes the stationary distribution
    // of the generator D0 + D1
    return map_pie(H0, H1)
}

/**
 * Computes the stationary initial probability for an ME/RAP distribution using matrices stored in a MatrixCell.
 *
 * @param ME The Matrix Exponential/RAP distribution stored in a MatrixCell
 * @return The stationary initial probability vector
 */
fun me_pie(ME: MatrixCell): Matrix {
    return map_pie(ME[0], ME[1])
}

/**
 * ME initial probability computation algorithms documentation marker for Dokka.
 */
@Suppress("unused")
class MePie {
    companion object {
        // Class documentation marker for Dokka
    }
}
