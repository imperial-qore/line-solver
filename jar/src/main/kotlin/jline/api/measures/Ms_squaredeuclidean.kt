package jline.api.measures

import jline.util.matrix.Matrix

/**
 * Squared Euclidean distance between two probability distributions.
 * Part of the squared L2 or chi-squared family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Squared Euclidean distance
 */
fun ms_squaredeuclidean(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val diff = P.get(i) - Q.get(i)
        sum += diff * diff
    }
    
    return sum
}
/**
 * Squaredeuclidean metric algorithms
 */
@Suppress("unused")
class MsSquaredeuclideanAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}