/**
 * @file Hellinger distance for probability distributions
 * 
 * Computes the Hellinger distance measuring dissimilarity between probability distributions
 * based on the square root of their densities. Bounded between 0 and √2, it's commonly
 * used in statistics for comparing continuous and discrete probability distributions.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.sqrt

/**
 * Hellinger distance between two probability distributions.
 * Part of the fidelity family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Hellinger distance
 */
fun ms_hellinger(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val diff = sqrt(P.get(i)) - sqrt(Q.get(i))
        sum += diff * diff
    }
    
    return sqrt(2.0 * sum)
}
/**
 * Hellinger metric algorithms
 */
@Suppress("unused")
class MsHellingerAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}