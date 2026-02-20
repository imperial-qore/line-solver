/**
 * @file Normalizing Constant Approximation (NCA) for single-class networks
 * 
 * Implements the Normalizing Constant Approximation method for single-class closed
 * queueing networks with load-dependent service rates. Uses iterative computation
 * with service rate functions to approximate the normalizing constant.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.nc

import jline.util.matrix.Matrix

fun pfqn_nca(L: Matrix, N: Int, alpha: Matrix): Matrix {
    val M = L.numRows
    val R = L.numCols
    if (R > 1 && M > 1) {
        throw RuntimeException("pfqn_nca supports only single class queuing networks; L must be a column vector.")
    }
    val X = Matrix.ones(M, N)
    val pi = Matrix(M, 1 + N)
    pi.zero()
    for (j in 0..<N) {
        X[0, j] = alpha[0, j] / L[0]
    }
    for (i in 1..<M) {
        pi[i, 0] = 1.0
        for (k in 0..<N) {
            for (j in 0..<k) {
                pi[i, 1 + j] = pi[i - 1, j] * L[i] / alpha[i, j]
            }
            pi[i, 1] = pi[i, 1] / X[i - 1, k]
            X[i, k] = 1 / pi.sumSubMatrix(i, i, 1, k + 1)
            for (j in -1..<k) {
                pi[i, j + 1] = X[i, k] * pi[i, j + 1]
            }
        }
    }
    return X
}
/**
 * PFQN nca algorithms
 */
@Suppress("unused")
class PfqnNcaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}