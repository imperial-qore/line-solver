/**
 * @file Exact Recursive Cache Probability Computation
 * 
 * Computes exact cache state probabilities using recursive methods based on 
 * exact recursive (EREC) algorithms. Provides precise analysis for small to 
 * medium-sized cache systems where computational complexity is manageable.
 * 
 * @since LINE 3.0
 */
package jline.api.cache

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.stream.IntStream

/**
 * Computes the cache state probabilities using an exact recursive method.
 * This method calculates the probabilities of the cache being in different states based on the cache
 * access factors and capacity.
 *
 * @param gamma Matrix representing the cache access factors.
 * @param m Matrix representing the cache capacity vector.
 * @return matrix containing the computed cache state probabilities for each item and level.
 */
fun cache_prob_erec(gamma: Matrix, m: Matrix): Matrix {
    val n = gamma.numRows
    val h = gamma.numCols
    val E = cache_erec(gamma, m)
    val prob = Matrix(n, h + 1)

    for (i in 0..<n) {
        for (j in 0..<h) {
            val subGamma = gamma.copy()
            subGamma.removeRows(setOf(i))
            val Ei = cache_erec(subGamma, Matrix.oner(m, j))
            val value = m[j] * gamma[i, j] * Ei.value() / E.value()
            prob[i, j + 1] = value
        }
        val finalI = i
        val rowSum = IntStream.range(1, h + 1).mapToDouble { j: Int -> prob[finalI, j] }.sum()
        prob[i, 0] = FastMath.abs(1 - rowSum)
    }
    return prob
}
/**
 * Cache prob erec algorithms
 */
@Suppress("unused")
class CacheProbErecAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}