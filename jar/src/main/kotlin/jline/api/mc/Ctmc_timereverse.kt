/**
 * @file CTMC time-reversal transformation
 * 
 * Computes the infinitesimal generator of the time-reversed continuous-time 
 * Markov chain using detailed balance equations. Time-reversal is fundamental 
 * in queueing theory and statistical mechanics for analyzing equilibrium properties.
 * 
 * @since LINE 3.0
 */
package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Compute the infinitesimal generator of the time-reserved CTMC
 *
 * @param Q Infinitesimal generator of the CTMC
 * @return Infinitesimal generator of the time-reversed CTMC
 */

fun ctmc_timereverse(Q: Matrix): Matrix {
    val piq = ctmc_solve(Q)
    val Qrev = Matrix(Q.numCols, Q.numRows)
    for (i in 0..<Q.numRows) {
        for (j in 0..<Q.numCols) {
            Qrev[i, j] = Q[i, j] * piq[i] / piq[j]
        }
    }
    return Qrev.transpose()
}
/**
 * CTMC timereverse algorithms
 */
@Suppress("unused")
class CtmcTimereverseAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}