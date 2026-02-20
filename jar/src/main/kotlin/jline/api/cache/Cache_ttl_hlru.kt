/**
 * @file TTL-based H-LRU Cache Analysis
 * 
 * Implements TTL approximation for Hierarchical LRU (H-LRU) cache systems.
 * Combines time-to-live eviction modeling with hierarchical LRU replacement
 * policies for multi-level cache performance analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Solve hierarchical least-recently-used caches h-LRU using the TTL approximation.
 *
 * @param lambda - Array of matrices representing request arrival rates from users to items of individual lists.
 * @param m      - Matrix representing cache capacity vector.
 * @return Matrix - A matrix containing the computed steady-state probabilities for the cache.
 */

fun cache_ttl_hlru(lambda: Array<Matrix>, m: Matrix): Matrix {
    val lambda2 = Matrix(lambda[0].numRows, lambda[0].numCols - 1)
    for (i in 0..<lambda[0].numRows) {
        for (j in 0..<lambda[0].numCols - 1) {
            lambda2[i, j] = lambda[0][i, j]
        }
    }


    val n = lambda2.numRows
    val h = m.numCols

    val t = cache_t_hlru(lambda2, m)

    val probh = Matrix(n, h) // steady state at list h
    val prob0 = Matrix(n, 1) // steady state at list 0
    val temp1 = Matrix(n, 1)
    temp1.ones()
    val temp2 = Matrix(n, 1)


    for (k in 0..<n) {
        for (s in 0..<h) {
            temp1[k, 0] = temp1[k, 0] * (1 - FastMath.exp(-lambda2[k, s] * t[s]))
        }

        var middtemp = 1.0
        var middtemp2 = 0.0
        for (l in 0..<h - 1) {
            for (s in 0..l) {
                middtemp = middtemp * (1 - FastMath.exp(-lambda2[k, s] * t[s]))
            }
            middtemp2 += middtemp
        }
        temp2[k, 0] = FastMath.exp(-lambda2[k, h - 1] * t[h - 1]) * (1 + middtemp2)

        probh[k, 0] = temp1[k, 0] / (temp1[k, 0] + temp2[k, 0])
        prob0[k, 0] = 1 / (temp1[k, 0] + temp2[k, 0]) / FastMath.exp(lambda2[k, h - 1] * t[0, h - 1])
    }


    return Matrix.concatColumns(prob0, probh, Matrix(n, h + 1))
}
/**
 * Cache ttl hlru algorithms
 */
@Suppress("unused")
class CacheTtlHlruAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}