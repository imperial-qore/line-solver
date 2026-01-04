/**
 * @file Cache Probability Analysis via Fixed Point Iteration
 * 
 * Computes cache state probabilities using fixed-point iteration algorithms.
 * Provides iterative solution methods for cache systems where direct
 * analytical solutions are not feasible due to system complexity.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.util.matrix.Matrix

/**
 * Estimate asymptotic values of the cache state probabilities at steady-state.
 *
 * @param gamma - Matrix representing the cache access factors.
 * @param m     - Matrix representing the cache capacity vector.
 * @return Matrix - A matrix containing the estimated cache state probabilities at equilibrium.
 */

fun cache_prob_fpi(gamma: Matrix, m: Matrix): Matrix {
    // FPI method
    val n = gamma.numRows
    val h = gamma.numCols

    val xi = cache_xi_fp(gamma, m, null).xi
    val prob = Matrix(n, h + 1)
    for (i in 0..<n) {
        val mul = Matrix.extractRows(gamma, i, i + 1, null).mult(xi.columnMajorOrder(), null)
        prob[i, 0] = 1 / (1 + mul[0])
        for (j in 1..<h + 1) {
            prob[i, j] = mul[0] / (1 + mul[0])
        }
    }
    return prob
}
/**
 * Cache prob fpi algorithms
 */
@Suppress("unused")
class CacheProbFpiAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}