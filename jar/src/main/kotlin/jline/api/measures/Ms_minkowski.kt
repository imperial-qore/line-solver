/**
 * Minkowski Distance for Probability Distributions
 * 
 * Implements the generalized Minkowski distance (Lp norm) between probability distributions.
 * Special cases include Manhattan distance (p=1) and Euclidean distance (p=2). Widely
 * used in machine learning for measuring dissimilarity between feature vectors.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.pow

/**
 * Minkowski distance between two probability distributions.
 * Generalization of Euclidean (p=2) and Manhattan (p=1) distances.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution  
 * @param p order parameter (p >= 1)
 * @return Minkowski distance of order p
 */
fun ms_minkowski(P: Matrix, Q: Matrix, p: Double): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    require(p >= 1.0) { "Order parameter p must be >= 1" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        sum += abs(P.get(i) - Q.get(i)).pow(p)
    }
    
    return sum.pow(1.0 / p)
}
/**
 * Minkowski metric algorithms
 */
@Suppress("unused")
class MsMinkowskiAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}