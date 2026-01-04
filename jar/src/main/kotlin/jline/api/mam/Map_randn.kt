/**
 * @file Markovian Arrival Process random generation with noise
 * 
 * Generates random MAP samples with added numerical noise for robustness testing.
 * Used for MAP validation, sensitivity analysis, and stochastic simulation studies.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.*
import kotlin.math.abs

/**
 * Generates a random Markovian Arrival Process (MAP) with 2 states using normal distribution
 * with mean 1 and standard deviation 2.
 *
 * @return a MatrixCell representing the random MAP transition matrices
 */
fun map_randn(): MatrixCell {
    return map_randn(2, 1.0, 2.0)
}

/**
 * Generates a random Markovian Arrival Process (MAP) with K states using normal distribution.
 *
 * @param K     the number of states
 * @param mu    the mean of the normal distribution
 * @param sigma the standard deviation of the normal distribution
 * @return a MatrixCell representing the random MAP transition matrices
 */

fun map_randn(K: Int, mu: Double, sigma: Double): MatrixCell {
    val random = Random()
    // Randomize D0
    val D0 = ArrayList<List<Double>>()
    for (i in 0..<K) {
        val row: MutableList<Double> = ArrayList()
        for (j in 0..<K) {
            row.add(abs(mu + sigma * random.nextGaussian()))
        }
        D0.add(row)
    }
    // Randomize D1
    val D1 = ArrayList<List<Double>>()
    for (i in 0..<K) {
        val row: MutableList<Double> = ArrayList()
        for (j in 0..<K) {
            row.add(abs(mu + sigma * random.nextGaussian()))
        }
        D1.add(row)
    }
    return map_normalize(Matrix(D0), Matrix(D1))
}
/**
 * MAP randn algorithms
 */
@Suppress("unused")
class MapRandnAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}