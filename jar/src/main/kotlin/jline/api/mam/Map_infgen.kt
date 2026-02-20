/**
 * @file Markovian Arrival Process infinitesimal generator computation
 * 
 * Computes infinitesimal generator matrix of underlying CTMC by combining MAP transition
 * matrices. Fundamental operation for all MAP analysis requiring continuous-time Markov chains.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the infinitesimal generator matrix (Q) of the Continuous-Time Markov Chain (CTMC) underlying a Markovian Arrival Process (MAP).
 *
 *
 * This function calculates the Q matrix by summing the hidden (D0) and visible (D1) transition matrices of the MAP.
 * The resulting matrix Q represents the infinitesimal generator of the CTMC, which describes the rate of transitions between states.
 *
 *
 * @param D0 The hidden transition matrix of the MAP, representing transitions without visible events.
 * @param D1 The visible transition matrix of the MAP, representing transitions with visible events.
 * @return The CTMC infinitesimal generator matrix Q.
 */
fun map_infgen(D0: Matrix, D1: Matrix): Matrix {
    return D0.add(1.0, D1)
}

/**
 * Computes the infinitesimal generator matrix (Q) of the Continuous-Time Markov Chain (CTMC) underlying a Markovian Arrival Process (MAP).
 *
 *
 * This is a convenience method that extracts the D0 and D1 matrices from a given MAP stored in a MatrixCell and computes the Q matrix.
 * The resulting matrix Q represents the infinitesimal generator of the CTMC, describing the rate of transitions between states.
 *
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices.
 * @return The CTMC infinitesimal generator matrix Q.
 */
fun map_infgen(MAP: MatrixCell): Matrix {
    return map_infgen(MAP[0], MAP[1])
}
/**
 * MAP infgen algorithms
 */
@Suppress("unused")
class MapInfgenAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}