package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Compute the infinitesimal generator of the time-reversed DTMC.
 *
 * @param P Infinitesimal generator of the DTMC
 * @return Infinitesimal generator of the time-reversed DTMC
 */

fun dtmc_timereverse(P: Matrix): Matrix {
    val pie = dtmc_solve(P)
    val Prev = Matrix(P.numCols, P.numRows)

    for (i in 0..<P.numRows) {
        for (j in 0..<P.numCols) {
            Prev[i, j] = P[i, j] * pie[i] / pie[j]
        }
    }
    return Prev.transpose()
}
/**
 * DTMC timereverse algorithms
 */
@Suppress("unused")
class DtmcTimereverseAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}