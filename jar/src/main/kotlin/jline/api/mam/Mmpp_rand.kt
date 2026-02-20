/**
 * @file Markov Modulated Poisson Process random generation
 * 
 * Generates random MMPP models with diagonal D1 matrices for testing and simulation.
 * Essential for validation and stochastic analysis of MMPP algorithms.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.*

/**
 * Generates a random Markov Modulated Poisson Process (MMPP) with K states.
 *
 * This function creates a random MMPP by generating random D0 and D1 matrices,
 * where D1 is constrained to be diagonal (characteristic of MMPPs).
 *
 * @param K the number of states
 * @return a MatrixCell representing the random MMPP transition matrices
 */
fun mmpp_rand(K: Int): MatrixCell {
    val random = Random()
    
    // Generate random D0 matrix
    val D0 = Matrix(K, K)
    for (i in 0..<K) {
        for (j in 0..<K) {
            D0.set(i, j, random.nextDouble())
        }
    }
    
    // Generate random D1 matrix (diagonal for MMPP)
    val D1 = Matrix(K, K)
    for (i in 0..<K) {
        for (j in 0..<K) {
            if (i == j) {
                D1.set(i, j, random.nextDouble())
            } else {
                D1.set(i, j, 0.0)
            }
        }
    }
    
    return map_normalize(D0, D1)
}

/**
 * Generates a random Markov Modulated Poisson Process (MMPP) with 2 states.
 *
 * @return a MatrixCell representing the random MMPP transition matrices
 */
fun mmpp_rand(): MatrixCell {
    return mmpp_rand(2)
}
/**
 * Mmpp Rand algorithms
 */
@Suppress("unused")
class MmppRandAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}