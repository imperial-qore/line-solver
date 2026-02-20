/**
 * @file Markovian Arrival Process embedded DTMC steady-state analysis
 * 
 * Computes steady-state probability vector of embedded discrete-time Markov chain.
 * Fundamental for analyzing long-term behavior and equilibrium properties of MAP processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the steady-state probability vector of the embedded Discrete Time Markov Chain (DTMC)
 * associated with a Markovian Arrival Process (MAP).
 *
 *
 * The MAP is represented by two matrices: D0 and D1. D0 is the hidden transition matrix, representing
 * transitions without an observed event, while D1 is the visible transition matrix, representing transitions
 * with an observed event.
 *
 * @param D0 the hidden transition matrix of the MAP
 * @param D1 the visible transition matrix of the MAP
 * @return the embedded steady-state probability vector
 */

fun map_pie(D0: Matrix, D1: Matrix): Matrix {
    val e = Matrix.ones(D0.numRows, 1)
    val A = map_piq(D0, D1).mult(D1) // piq*D1
    A.scaleEq(1.0 / A.mult(e).toDouble()) // divide by A*ones(length(MAP{2}),1)
    return A
}

/**
 * Computes the steady-state probability vector of the embedded DTMC of a MAP
 * stored in a MatrixCell that contains the MAP's transition matrices.
 *
 * @param MAP a MatrixCell containing the transition matrices D0 and D1 of the MAP
 * @return the embedded steady-state probability vector
 */

fun map_pie(MAP: MatrixCell): Matrix {
    return map_pie(MAP[0], MAP[1])
}
/**
 * MAP pie algorithms
 */
@Suppress("unused")
class MapPieAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}