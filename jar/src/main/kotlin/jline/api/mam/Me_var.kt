package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * ME variance computation algorithms.
 *
 * Provides methods for computing the variance of Matrix Exponential (ME) distributions.
 * The variance is calculated from the initial vector alpha and matrix parameter A.
 *
 * @since LINE 3.0
 */

/**
 * Computes the variance of a Matrix Exponential (ME) distribution.
 *
 * For an ME distribution with initial vector alpha and matrix parameter A,
 * the variance is computed using the first two moments:
 * - m1 = -alpha * A^(-1) * e
 * - m2 = 2 * alpha * A^(-2) * e
 * - variance = m2 - m1^2
 *
 * @param alpha The initial vector of the ME distribution
 * @param A The matrix parameter of the ME distribution
 * @return The variance of the ME distribution
 */
fun me_var(alpha: Matrix, A: Matrix): Double {
    val n = A.numRows
    val e = Matrix.ones(n, 1)
    val Ainv = A.inv()
    val Ainv2 = Ainv.mult(Ainv)

    val m1 = -alpha.mult(Ainv).mult(e).get(0, 0)
    val m2 = 2.0 * alpha.mult(Ainv2).mult(e).get(0, 0)

    return m2 - m1 * m1
}

/**
 * Computes the variance of a Matrix Exponential (ME) distribution using matrices stored in a MatrixCell.
 *
 * This is a convenience method that computes the variance using the process representation.
 *
 * @param ME The Matrix Exponential distribution stored in a MatrixCell
 * @return The variance of the ME distribution
 */
fun me_var(ME: MatrixCell): Double {
    // Use map_var on the process representation
    return map_var(ME[0], ME[1])
}

/**
 * ME variance computation algorithms documentation marker for Dokka.
 */
@Suppress("unused")
class MeVar {
    companion object {
        // Class documentation marker for Dokka
    }
}
