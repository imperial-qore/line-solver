package jline.api.measures

import jline.util.matrix.Matrix

/**
 * Divergence distance between two probability distributions.
 * Part of the squared L2 or chi-squared family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Divergence distance
 */
fun ms_divergence(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val sumPQ = pi + qi
        if (sumPQ > 0) {
            val diff = pi - qi
            sum += (diff * diff) / (sumPQ * sumPQ)
        }
    }
    
    return 2.0 * sum
}
/**
 * Divergence metric algorithms
 */
@Suppress("unused")
class MsDivergenceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}