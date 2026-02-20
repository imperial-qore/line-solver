/**
 * @file Load-dependent scaling vector position shifting utility
 * 
 * Provides utility function for shifting load-dependent scaling vectors by one position,
 * commonly used in recursive normalizing constant computations. Supports the manipulation
 * of state-dependent service rate vectors during iterative solution algorithms.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.util.matrix.Matrix

/**
 * Shifts a load-dependent scaling vector by one position
 *
 * @param mu - load-dependent scalings
 * @return normalizing constant and its logarithm
 */

fun pfqn_mushift(mu: Matrix, k: Int): Matrix {
    val M = mu.numRows
    val N = mu.numCols
    val mushift = Matrix(M, N - 1)
    Matrix.extract(mu, 0, M, 0, N - 1, mushift, 0, 0)
    for (j in 0..<N - 1) {
        mushift[k, j] = mu[k, j + 1]
    }
    return mushift
}
/**
 * PFQN mushift algorithms
 */
@Suppress("unused")
class PfqnMushiftAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}