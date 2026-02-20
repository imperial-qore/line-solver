/**
 * @file DTMC stochastic matrix normalization
 * 
 * Converts non-negative matrices into valid discrete-time Markov chain transition 
 * matrices by row normalization. Ensures row stochasticity property (row sums = 1) 
 * required for valid DTMC transition matrices.
 * 
 * @since LINE 3.0
 */
package jline.api.mc

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.math.max

/**
 * Normalize a given non-negative matrix into a DTMC
 *
 * @param P nonegative matrix
 * @return Transition matrix of the DTMC
 */

fun dtmc_makestochastic(P: Matrix): Matrix {
    var P = P
    P = P.copy()
    val n = P.length()
    for (i in 0..<n) {
        var rowSum = P.getRow(i).elementSum()
        if (rowSum > 0) {
            for (j in 0..<n) {
                P[i, j] = P[i, j] / rowSum
            }
            rowSum = P.getRow(i).elementSum()
            P[i, i] = FastMath.min(max(0.0, 1.0 - (rowSum - P[i, i])), 1.0)
        } else {
            for (j in 0..<n) {
                P[i, j] = 0.0
            }
            P[i, i] = 1.0
        }
    }
    return P
}
/**
 * DTMC makestochastic algorithms
 */
@Suppress("unused")
class DtmcMakestochasticAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}