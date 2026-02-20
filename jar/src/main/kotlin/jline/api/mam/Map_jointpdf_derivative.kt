/**
 * @file Markovian Arrival Process joint PDF derivative analysis
 * 
 * Computes partial derivatives of MAP joint probability density functions at origin.
 * Used for advanced statistical analysis and theoretical characterization of arrival dependencies.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Compute partial derivative at 0 of a MAP's joint PDF
 * 
 * Based on: A. Horvath et al. A Joint Moments Based Analysis of Networks of
 * MAP/MAP/1 Queues
 *
 * @param MAP MatrixCell containing D0 and D1 matrices of the MAP
 * @param iset Index set for the derivative computation
 * @return Scalar value representing the derivative
 */
fun map_jointpdf_derivative(MAP: MatrixCell, iset: IntArray): Double {
    val D0 = MAP[0]
    val D1 = MAP[1]
    val n = D0.numRows
    
    var gamma = map_pie(MAP)
    
    for (j in iset) {
        // gamma = gamma * MAP{1}^j * MAP{2}
        val D0PowerJ = Matrix.pow(D0, j)
        gamma = gamma.mult(D0PowerJ).mult(D1)
    }
    
    val ones = Matrix.ones(n, 1)
    return gamma.mult(ones).get(0, 0)
}

/**
 * Compute partial derivative at 0 of a MAP's joint PDF
 *
 * @param D0 Hidden transition matrix of the MAP
 * @param D1 Visible transition matrix of the MAP  
 * @param iset Index set for the derivative computation
 * @return Scalar value representing the derivative
 */
fun map_jointpdf_derivative(D0: Matrix, D1: Matrix, iset: IntArray): Double {
    val MAP = MatrixCell(arrayOf(D0, D1))
    return map_jointpdf_derivative(MAP, iset)
}
/**
 * MAP jointpdf derivative algorithms
 */
@Suppress("unused")
class MapJointpdfDerivativeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}