package jline.api.measures

import jline.util.matrix.Matrix

/**
 * Pearson chi-squared distance between two probability distributions.
 * Part of the squared L2 or chi-squared family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Pearson chi-squared distance
 */
fun ms_pearsonchisquared(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        if (qi > 0) {
            val diff = pi - qi
            sum += (diff * diff) / qi
        }
    }
    
    return sum
}
/**
 * Pearsonchisquared metric algorithms
 */
@Suppress("unused")
class MsPearsonchisquaredAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}