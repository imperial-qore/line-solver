package jline.api.measures

import jline.util.matrix.Matrix

/**
 * Probabilistic symmetry chi-squared distance between two probability distributions.
 * Part of the squared L2 or chi-squared family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Probabilistic symmetry chi-squared distance
 */
fun ms_probsymmchisquared(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val denom = pi + qi
        if (denom > 0) {
            val diff = pi - qi
            sum += (diff * diff) / denom
        }
    }
    
    return 2.0 * sum
}