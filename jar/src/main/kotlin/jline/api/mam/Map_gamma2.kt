/**
 * @file Markovian Arrival Process eigenvalue-based correlation analysis
 * 
 * Computes largest non-unit eigenvalue of embedded DTMC for MAP correlation characterization.
 * Used for advanced analysis of temporal dependence patterns in Markovian arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell
import org.apache.commons.math3.linear.EigenDecomposition
import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.util.FastMath

/**
 * Returns the largest non-unit eigenvalue (both real and imaginary parts) of the
 * embedded Discrete-Time Markov Chain (DTMC) of a given Markovian Arrival Process (MAP).
 *
 *
 * The function first computes the embedded DTMC transition matrix from the provided MAP.
 * It then calculates the eigenvalues of this matrix and identifies the largest non-unit
 * eigenvalue based on its magnitude. The eigenvalue's real and imaginary parts are returned.
 *
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell.
 * @return A double array of size 2, where the first element is the real part and the second element
 * is the absolute value of the imaginary part of the largest non-unit eigenvalue.
 */
fun map_gamma2(MAP: MatrixCell): DoubleArray {
    val P = map_embedded(MAP)
    val realMatrix = MatrixUtils.createRealMatrix(P.toArray2D())

    val eigDecomp = EigenDecomposition(realMatrix)
    val realPart = eigDecomp.realEigenvalues
    val imgPart = eigDecomp.imagEigenvalues

    // Find the index of the maximum absolute value
    var maxIndex = realPart.size - 1
    var secondMaxIndex = realPart.size - 1
    var maxAbsValue = -Double.MIN_VALUE
    for (i in realPart.indices) {
        val abs = FastMath.sqrt(realPart[i] * realPart[i] + imgPart[i] * imgPart[i])
        if (abs > maxAbsValue) {
            maxAbsValue = abs
            secondMaxIndex = maxIndex
            maxIndex = i
        }
    }

    // Return the real and imaginary parts of the eigenvalue with the maximum real part
    return doubleArrayOf(realPart[secondMaxIndex], FastMath.abs(imgPart[secondMaxIndex]))
}
/**
 * MAP gamma2 algorithms
 */
@Suppress("unused")
class MapGamma2 {
    companion object {
        // Class documentation marker for Dokka
    }
}