/**
 * @file Markovian Arrival Process cumulative distribution function
 * 
 * Computes CDF values for MAP inter-arrival times using CTMC uniformization techniques.
 * Essential for probability analysis and performance evaluation of stochastic arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.api.mc.ctmc_uniformization
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the cumulative distribution function (CDF) of the inter-arrival times of a Markovian Arrival Process (MAP).
 *
 *
 * The MAP is represented by two matrices: D0 and D1. D0 is the hidden transition matrix, representing
 * transitions without an observed event, while D1 is the visible transition matrix, representing transitions
 * with an observed event. The CDF values are calculated for a given set of points, which represent the inter-arrival times.
 *
 * @param D0     the hidden transition matrix of the MAP
 * @param D1     the visible transition matrix of the MAP
 * @param points a matrix containing the points at which to compute the CDF
 * @return a matrix containing the CDF values corresponding to the provided points
 */

fun map_cdf(D0: Matrix, D1: Matrix, points: Matrix): Matrix {
    val CDFVals = Matrix(1, points.length())
    val pie = map_pie(D0, D1)
    val e1 = Matrix.ones(D0.numRows, 1)

    var nanVal = 0.0
    for (t in 0..<points.length()) {
        val output = ctmc_uniformization(pie, D0, points[t]).mult(e1)
        var `val` = 1 - output.value()
        if (java.lang.Double.isNaN(`val`)) {
            `val` = nanVal
        } else { // after it finds the first non-zero, set nanVal to 1.0
            nanVal = 1.0
        }
        CDFVals[0, t] = `val`
    }
    return CDFVals
}

/**
 * Computes the cumulative distribution function (CDF) of the inter-arrival times of a MAP
 * stored in a MatrixCell that contains the MAP's transition matrices.
 *
 * @param MAP    a MatrixCell containing the transition matrices D0 and D1 of the MAP
 * @param points a matrix containing the points at which to compute the CDF
 * @return a matrix containing the CDF values corresponding to the provided points
 */

fun map_cdf(MAP: MatrixCell, points: Matrix): Matrix {
    return map_cdf(MAP[0], MAP[1], points)
}

/**
 * MAP cumulative distribution function algorithms
 */
@Suppress("unused")
class MapCdf {
    companion object {
        // Class documentation marker for Dokka
    }
}