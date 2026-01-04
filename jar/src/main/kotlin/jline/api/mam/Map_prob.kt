/**
 * @file Markovian Arrival Process equilibrium distribution computation
 * 
 * Computes equilibrium probability distribution of underlying CTMC for MAP analysis.
 * Essential for steady-state analysis and fundamental MAP performance calculations.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.api.mc.ctmc_solve
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Computes the equilibrium distribution of the underlying continuous-time Markov chain for a MAP.
 *
 * This function calculates the steady-state probabilities of the continuous-time Markov chain
 * underlying a Markovian Arrival Process (MAP). The equilibrium distribution is computed
 * by solving the system Q*π = 0, where Q = D0 + D1 is the infinitesimal generator matrix.
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices
 * @return The equilibrium distribution as a Matrix (row vector)
 */
fun map_prob(MAP: MatrixCell): Matrix {
    return ctmc_solve(map_infgen(MAP))
}

/**
 * Computes the equilibrium distribution of the underlying continuous-time Markov chain for a MAP.
 *
 * This function calculates the steady-state probabilities of the continuous-time Markov chain
 * underlying a Markovian Arrival Process (MAP) given the D0 and D1 matrices directly.
 *
 * @param D0 The hidden transition matrix of the MAP
 * @param D1 The visible transition matrix of the MAP
 * @return The equilibrium distribution as a Matrix (row vector)
 */
fun map_prob(D0: Matrix, D1: Matrix): Matrix {
    return ctmc_solve(map_infgen(D0, D1))
}
/**
 * MAP prob algorithms
 */
@Suppress("unused")
class MapProbAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}