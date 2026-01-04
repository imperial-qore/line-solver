package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Czekanowski distance between two probability distributions.
 * Part of the intersection family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Czekanowski distance
 */
fun ms_czekanowski(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sumDiff = 0.0
    var sumTotal = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        sumDiff += abs(pi - qi)
        sumTotal += pi + qi
    }
    
    return if (sumTotal > 0) sumDiff / sumTotal else 0.0
}
/**
 * Czekanowski metric algorithms
 */
@Suppress("unused")
class MsCzekanowskiAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}