/**
 * @file Discrete-time Markov chain steady-state solver
 * 
 * Computes the steady-state probability distribution for DTMCs by converting the
 * transition matrix problem (P-I)x = 0 into a CTMC-equivalent system and leveraging
 * the robust CTMC solver with automatic reducibility handling.
 *
 * @since LINE 3.0
 */
package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Returns the steady-state solution of a DTMC.
 *
 * @param P Transition matrix of the DTMC
 * @return Steady-state solution vector of the DTMC
 */

fun dtmc_solve(P: Matrix): Matrix {
    val Plocal = P.copy()
    //P-eye(size(P))
    for (i in 0..<Plocal.numRows) {
        Plocal[i, i] = Plocal[i, i] - 1.0
    }
    return ctmc_solve(Plocal)
}
/**
 * DTMC solve algorithms
 */
@Suppress("unused")
class DtmcSolveAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}