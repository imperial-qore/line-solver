/**
 * @file Markovian Arrival Process raw moment computation
 * 
 * Computes raw moments of MAP inter-arrival times using matrix inversion techniques.
 * Essential for statistical characterization and parameter fitting of arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the raw moments of the inter-arrival times of a Markovian Arrival Process (MAP).
 *
 *
 * The raw moment of a given order provides a measure of the shape of the distribution of
 * inter-arrival times. The MAP is represented by two matrices: D0 and D1. D0 is the hidden transition matrix,
 * representing transitions without an observed event, while D1 is the visible transition matrix, representing
 * transitions with an observed event.
 *
 * @param D0    the hidden transition matrix of the MAP
 * @param D1    the visible transition matrix of the MAP
 * @param order the moment order, i.e., E[X^order], where X is the inter-arrival time
 * @return the raw moment of the inter-arrival times of the specified order
 */

fun map_moment(D0: Matrix, D1: Matrix, order: Int): Double {
    // Handle case where D0 is zero matrix or singular, similar to MATLAB implementation
    // Check if all elements are zero or if determinant is effectively zero
    var isZeroMatrix = true
    for (i in 0 until D0.numRows) {
        for (j in 0 until D0.numCols) {
            if (Math.abs(D0.get(i, j)) > 1e-14) {
                isZeroMatrix = false
                break
            }
        }
        if (!isZeroMatrix) break
    }

    if (isZeroMatrix || Math.abs(D0.det()) < 1e-12) {
        return 0.0
    }

    var pie = map_pie(D0, D1)
    var iD0 = D0.copy()
    iD0.scaleEq(-1.0)
    iD0 = iD0.inv()
    var iD0k = Matrix(iD0)
    for (i in 2..order) {
        iD0k = iD0k.mult(iD0) // scale(i) introduces the factorial
        iD0k.scaleEq(i.toDouble())
    }
    val e = Matrix.ones(D0.numRows, 1)
    pie = pie.mult(iD0k)
    pie = pie.mult(e)
    return pie.toDouble()
}

/**
 * Computes the raw moments of the inter-arrival times of a MAP
 * stored in a MatrixCell that contains the MAP's transition matrices.
 *
 * @param MAP   a MatrixCell containing the transition matrices D0 and D1 of the MAP
 * @param order the moment order, i.e., E[X^order], where X is the inter-arrival time
 * @return the raw moment of the inter-arrival times of the specified order
 */

fun map_moment(MAP: MatrixCell, order: Int): Double {
    return map_moment(MAP[0], MAP[1], order)
}
/**
 * MAP moment algorithms
 */
@Suppress("unused")
class MapMoment {
    companion object {
        // Class documentation marker for Dokka
    }
}