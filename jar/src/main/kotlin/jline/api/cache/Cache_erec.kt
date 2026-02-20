/**
 * @file Exact Recursive Cache Analysis
 * 
 * Implements exact recursive (EREC) algorithms for cache system analysis.
 * Provides numerically exact solutions for cache performance metrics in
 * systems where computational complexity allows for precise evaluation.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.util.matrix.Matrix

/**
 * Computes the cache miss rate using an exact recursive method.
 * This method serves as a wrapper that calls the auxiliary function to perform the actual computation.
 *
 * @param gamma - Matrix representing the cache access factors.
 * @param m     - Matrix representing the cache capacity vector.
 * @return Matrix - A matrix containing the computed cache miss rates.
 */

fun cache_erec(gamma: Matrix, m: Matrix): Matrix {
    return cache_erec_aux(gamma, m, gamma.numRows)
}

/**
 * Auxiliary method for computing the cache miss rate using an exact recursive method.
 * This method performs the core computation recursively, adjusting the size of the input matrix.
 *
 * @param gamma - Matrix representing the cache access factors.
 * @param m     - Matrix representing the cache capacity vector.
 * @param k     - Integer representing the current number of rows in the recursive step.
 * @return Matrix - A matrix containing the computed cache miss rates.
 */
fun cache_erec_aux(gamma: Matrix, m: Matrix, k: Int): Matrix {
    val h = m.numElements

    if (m.elementSum() == 0.0) {
        return Matrix.singleton(1.0)
    }

    if (m.elementSum() > k || m.elementMin() < 0) {
        return Matrix.singleton(0.0)
    }

    if (k == 1 && m.elementSum() == 1.0) {
        val j = m.find().value().toInt() // Find the index of the non-zero element in m
        return Matrix.singleton(gamma[0, j])
    }

    var E = cache_erec_aux(gamma, m, k - 1)
    for (j in 0..<h) {
        if (m[j] > 0) {
            val onerM = Matrix.oner(m, j)
            val term = cache_erec_aux(gamma, onerM, k - 1).scale(gamma[k - 1, j] * m[j])
            E = E.add(term)
        }
    }
    return E
}
/**
 * Cache erec algorithms
 */
@Suppress("unused")
class CacheErecAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}