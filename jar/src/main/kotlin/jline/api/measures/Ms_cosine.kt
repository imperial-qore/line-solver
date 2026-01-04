/**
 * @file Cosine distance for probability distributions
 * 
 * Computes cosine distance (1 - cosine similarity) measuring the angle between two
 * vectors in high-dimensional space. Particularly useful for text analysis and
 * information retrieval where magnitude is less important than direction.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.sqrt

/**
 * Cosine distance between two probability distributions.
 * Part of the inner product family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Cosine distance
 */
fun ms_cosine(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var dotProduct = 0.0
    var sumP2 = 0.0
    var sumQ2 = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        dotProduct += pi * qi
        sumP2 += pi * pi
        sumQ2 += qi * qi
    }
    
    val s = dotProduct / sqrt(sumP2) / sqrt(sumQ2)
    return 1.0 - s
}
/**
 * Cosine metric algorithms
 */
@Suppress("unused")
class MsCosineAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}