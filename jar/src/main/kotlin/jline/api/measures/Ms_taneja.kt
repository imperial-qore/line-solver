package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.ln
import kotlin.math.sqrt

/**
 * Taneja distance between two probability distributions.
 * Part of the combination family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Taneja distance
 */
fun ms_taneja(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val m = (pi + qi) / 2.0
        val sqrtPQ = sqrt(pi * qi)
        if (m > 0 && sqrtPQ > 0) {
            sum += m * ln(m / sqrtPQ)
        }
    }
    
    return sum
}
/**
 * Taneja metric algorithms
 */
@Suppress("unused")
class MsTanejaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}