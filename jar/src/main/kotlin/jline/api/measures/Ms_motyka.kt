package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.max

/**
 * Motyka distance between two probability distributions.
 * Part of the intersection family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Motyka distance
 */
fun ms_motyka(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sumMax = 0.0
    var sumTotal = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        sumMax += max(pi, qi)
        sumTotal += pi + qi
    }
    
    return if (sumTotal > 0) sumMax / sumTotal else 0.0
}
/**
 * Motyka metric algorithms
 */
@Suppress("unused")
class MsMotykaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}