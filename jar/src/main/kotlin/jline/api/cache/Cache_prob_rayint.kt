package jline.api.cache

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import java.util.stream.IntStream

/**
 * Computes the cache state probabilities using the ray method.
 * This method calculates the probabilities of the cache being in different states based on the access factors
 * and capacity, utilizing the logarithm of the partition function obtained through the ray method.
 *
 * @param gamma - Matrix representing the cache access factors.
 * @param m     - Matrix representing the cache capacity vector.
 * @return Matrix - A matrix containing the computed cache state probabilities for each item and level.
 */

fun cache_prob_rayint(gamma: Matrix, m: Matrix): Matrix {
    val n = gamma.numRows
    val h = gamma.numCols
    val lE = cache_rayint(gamma, m).lZ
    val prob = Matrix(n, h + 1)

    for (i in 0..<n) {
        for (j in 0..<h) {
            val subGamma = gamma.copy()
            subGamma.removeRows(setOf(i))
            val lEi = cache_rayint(subGamma, Matrix.oner(m, j)).lZ
            val value = m[j] * gamma[i, j] * FastMath.exp(lEi - lE)
            prob[i, j + 1] = value
        }
        val finalI = i
        val rowSum = IntStream.range(1, h + 1).mapToDouble { j: Int -> prob[finalI, j] }.sum()
        prob[i, 0] = FastMath.abs(1 - rowSum)
    }
    return prob
}