package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.sqrt

/**
 * Clark distance between two probability distributions.
 * Part of the squared L2 or chi-squared family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Clark distance
 */
fun ms_clark(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val denom = pi + qi
        if (denom > 0) {
            val term = abs(pi - qi) / denom
            sum += term * term
        }
    }
    
    return sqrt(sum)
}
/**
 * Clark metric algorithms
 */
@Suppress("unused")
class MsClarkAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}