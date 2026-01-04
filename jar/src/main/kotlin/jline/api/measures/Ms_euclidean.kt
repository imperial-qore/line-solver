/**
 * @file Euclidean distance metric
 * 
 * Implements the standard Euclidean distance d(p,q) = √∑(pᵢ-qᵢ)² between probability 
 * distributions or vectors. The most commonly used distance metric providing geometric 
 * interpretation of dissimilarity in statistical analysis and machine learning.
 * 
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.sqrt

/**
 * Euclidean distance between two probability distributions.
 * Part of the Minkowski family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Euclidean distance
 */
fun ms_euclidean(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val diff = P.get(i) - Q.get(i)
        sum += diff * diff
    }
    
    return sqrt(sum)
}

/**
 * Euclidean distance metric algorithms
 */
@Suppress("unused")
class MsEuclideanAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}