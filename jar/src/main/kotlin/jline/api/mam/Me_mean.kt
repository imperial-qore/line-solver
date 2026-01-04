package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * ME mean computation algorithms.
 *
 * Provides methods for computing the mean of Matrix Exponential (ME) distributions.
 * The mean is calculated from the initial vector alpha and matrix parameter A.
 *
 * @since LINE 3.0
 */

/**
 * Computes the mean of a Matrix Exponential (ME) distribution.
 *
 * For an ME distribution with initial vector alpha and matrix parameter A,
 * the mean is computed as: mean = -alpha * A^(-1) * e
 * where e is a vector of ones.
 *
 * @param alpha The initial vector of the ME distribution
 * @param A The matrix parameter of the ME distribution
 * @return The mean of the ME distribution
 */
fun me_mean(alpha: Matrix, A: Matrix): Double {
    val n = A.numRows
    val e = Matrix.ones(n, 1)
    val Ainv = A.inv()
    return -alpha.mult(Ainv).mult(e).get(0, 0)
}

/**
 * Computes the mean of a Matrix Exponential (ME) distribution using matrices stored in a MatrixCell.
 *
 * This is a convenience method that extracts alpha and A from a given ME stored in a MatrixCell
 * and computes the mean.
 *
 * Note: For ME distributions, the MatrixCell format is {D0=A, D1=-A*e*alpha'}.
 * We extract A from D0.
 *
 * @param ME The Matrix Exponential distribution stored in a MatrixCell
 * @return The mean of the ME distribution
 */
fun me_mean(ME: MatrixCell): Double {
    // For ME: D0 = A
    // We need to reconstruct alpha from D1 = -A*e*alpha'
    // This is equivalent to using map_mean on the process representation
    return map_mean(ME[0], ME[1])
}

/**
 * ME mean computation algorithms documentation marker for Dokka.
 */
@Suppress("unused")
class MeMean {
    companion object {
        // Class documentation marker for Dokka
    }
}
