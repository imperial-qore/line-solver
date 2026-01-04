/**
 * @file Markovian Arrival Process complementary CDF derivative analysis
 * 
 * Computes derivatives of MAP complementary cumulative distribution functions at zero.
 * Used for advanced moment analysis and joint queue analysis in MAP/MAP/1 queueing systems.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Compute derivative at 0 of a MAP complementary CDF
 * 
 * Based on: A. Horvath et al. A Joint Moments Based Analysis of Networks of
 * MAP/MAP/1 Queues
 *
 * @param MAP MatrixCell containing D0 and D1 matrices of the MAP
 * @param i Derivative order
 * @return Scalar value representing the derivative
 */
fun map_ccdf_derivative(MAP: MatrixCell, i: Int): Double {
    val D0 = MAP[0]
    val n = D0.numRows
    
    val pie = map_pie(MAP)
    val D0PowerI = Matrix.pow(D0, i)
    val ones = Matrix.ones(n, 1)
    
    return pie.mult(D0PowerI).mult(ones).get(0, 0)
}

/**
 * Compute derivative at 0 of a MAP complementary CDF
 *
 * @param D0 Hidden transition matrix of the MAP
 * @param D1 Visible transition matrix of the MAP
 * @param i Derivative order
 * @return Scalar value representing the derivative
 */
fun map_ccdf_derivative(D0: Matrix, D1: Matrix, i: Int): Double {
    val MAP = MatrixCell(arrayOf(D0, D1))
    return map_ccdf_derivative(MAP, i)
}

/**
 * MAP complementary CDF derivative algorithms
 */
@Suppress("unused")
class MapCcdfDerivativeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}