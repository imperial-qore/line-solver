/**
 * @file TTL-based LRUM Cache Analysis
 * 
 * Implementation of Time-To-Live approximation for Least Recently Used with
 * Multiple servers (LRUM) cache systems. Extends TTL analysis to handle
 * multi-server cache architectures with LRU replacement policies.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Solve multi-list least-recently-used caches LRU(m) using the TTL approximation.
 *
 * @param gamma - Array of matrices representing access factors for items across different lists.
 * @param m     - Matrix representing cache capacity vector.
 * @return Matrix - A matrix containing the computed probabilities for the cache.
 */

fun cache_ttl_lrum(gamma: Array<Matrix?>, m: Matrix): Matrix {
    val gamma2 = Matrix(gamma[0]!!.numRows, gamma[0]!!.numCols - 1)
    for (i in 0..<gamma[0]!!.numRows) {  //10
        for (j in 0..<gamma[0]!!.numCols - 1) {  //2
            gamma2[i, j] = gamma[0]!![i, j]
        }
    }


    val n = gamma2.numRows
    val h = gamma2.numCols

    val t = cache_t_lrum(gamma2, m)

    val probh = Matrix(n, h) // steady state 1,...,h
    val prob0 = Matrix(n, 1) // steady state 0
    val trans = Matrix(n, h)
    val denom = Matrix(1, n)

    // Calculate transiction probabilities and denominators
    for (k in 0..<n) {
        for (j in 0..<h) {
            gamma2[k, j]
            t[0, j]
            trans[k, j] = FastMath.exp(gamma2[k, j] * t[0, j]) - 1 // birth and death rates
        }
        denom[0, k] = Matrix.sumCumprod(Matrix.extract(trans, k, k + 1, 0, h))
    }

    // Calculate steady-state probabilities
    for (k in 0..<n) {
        for (l in 0..<h) {
            probh[k, l] = (Matrix.extract(trans, k, k + 1, 0, l + 1)).elementMult() / (1 + denom[0, k])
        }
        val probhK = Matrix(1, h) //probh(k,:)
        for (i in 0..<h) {
            probhK[0, i] = probh[k, i]
        }
        prob0[k, 0] = 1 - probhK.elementSum()
    }

    return Matrix.concatColumns(prob0, probh, Matrix(n, h + 1))
}
/**
 * Cache ttl lrum algorithms
 */
@Suppress("unused")
class CacheTtlLrumAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}