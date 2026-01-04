/**
 * @file Markovian Arrival Process counting variance analysis
 * 
 * Computes variance of event counts in MAP processes over specified time intervals.
 * Essential for analyzing variability and dispersion in arrival counting processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.Maths
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the variance of the counts in a MAP over a time period t.
 *
 *
 * This method calculates the variance of the number of events that occur in the MAP
 * over the specified time interval t. It involves matrix operations that consider
 * the MAP's transition structure and steady-state behavior.
 *
 * @param D0 the hidden transition matrix of the MAP
 * @param D1 the visible transition matrix of the MAP
 * @param t  the time period over which to compute the variance of counts
 * @return the variance of the counts over time t
 */

fun map_varcount(D0: Matrix, D1: Matrix, t: Double): Double {
    val n = D0.numRows
    val Q = map_infgen(D0, D1)
    val I = Matrix.eye(n)
    val e = Matrix.ones(n, 1)
    val piq = map_piq(D0, D1)
    val lambda = 1 / map_mean(D0, D1)
    val IpiqQ = Matrix.ones(n, 1).mult(piq).add(-1.0, Q).inv()
    val PRE = (lambda - 2 * lambda * lambda + 2 * piq.mult(D1).mult(IpiqQ).mult(D1).mult(e).value())
    val POST = piq.mult(D1).mult(I.sub(Maths.matrixExp(Q.scale(t)))).mult(IpiqQ.square()).mult(D1).mult(e).value()
    return PRE * t - 2 * POST
}

/**
 * Computes the variance of the counts in a MAP over multiple time periods.
 *
 *
 * This method calculates the variance for each time period specified in the matrix t.
 *
 * @param D0 the hidden transition matrix of the MAP
 * @param D1 the visible transition matrix of the MAP
 * @param t  a matrix containing the time periods over which to compute the variance of counts
 * @return a matrix containing the variance of counts for each time period in t
 */
fun map_varcount(D0: Matrix, D1: Matrix, t: Matrix): Matrix {
    val result = Matrix(t.numRows, t.numCols)
    for (i in 0..<t.numElements) {
        result[i] = map_varcount(D0, D1, t[i])
    }
    return result
}

/**
 * Computes the variance of the counts in a MAP over a time period t using a MatrixCell.
 *
 *
 * This method is a convenience overload that extracts the D0 and D1 matrices from the provided MatrixCell
 * and calculates the variance of the counts over time t.
 *
 * @param MAP the MatrixCell representing the MAP, containing the D0 and D1 matrices
 * @param t   the time period over which to compute the variance of counts
 * @return the variance of the counts over time t
 */
fun map_varcount(MAP: MatrixCell, t: Double): Double {
    return map_varcount(MAP[0], MAP[1], t)
}

/**
 * Computes the variance of the counts in a MAP over multiple time periods using a MatrixCell.
 *
 *
 * This method is a convenience overload that extracts the D0 and D1 matrices from the provided MatrixCell
 * and calculates the variance of the counts over multiple time periods specified in the matrix t.
 *
 * @param MAP the MatrixCell representing the MAP, containing the D0 and D1 matrices
 * @param t   a matrix containing the time periods over which to compute the variance of counts
 * @return a matrix containing the variance of counts for each time period in t
 */
fun map_varcount(MAP: MatrixCell, t: Matrix): Matrix {
    return map_varcount(MAP[0], MAP[1], t)
}
/**
 * MAP varcount algorithms
 */
@Suppress("unused")
class MapVarcountAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}