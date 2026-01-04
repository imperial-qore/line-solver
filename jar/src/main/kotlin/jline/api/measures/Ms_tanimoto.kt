package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.max
import kotlin.math.min

/**
 * Tanimoto distance between two probability distributions.
 * Part of the intersection family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Tanimoto distance
 */
fun ms_tanimoto(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sumDiff = 0.0
    var sumMax = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val maxVal = max(pi, qi)
        val minVal = min(pi, qi)
        sumDiff += maxVal - minVal
        sumMax += maxVal
    }
    
    return if (sumMax > 0) sumDiff / sumMax else 0.0
}
/**
 * Tanimoto metric algorithms
 */
@Suppress("unused")
class MsTanimotoAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}