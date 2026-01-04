/**
 * @file Bhattacharyya distance for probability distributions
 * 
 * Computes the Bhattacharyya distance measuring the similarity between probability
 * distributions based on the Bhattacharyya coefficient. Widely used in pattern
 * recognition, classification, and feature selection for comparing statistical models.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.ln
import kotlin.math.sqrt

/**
 * Bhattacharyya distance between two probability distributions.
 * Part of the fidelity family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Bhattacharyya distance
 */
fun ms_bhattacharyya(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        sum += sqrt(P.get(i) * Q.get(i))
    }
    
    return -ln(sum)
}
/**
 * Bhattacharyya metric algorithms
 */
@Suppress("unused")
class MsBhattacharyyaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}