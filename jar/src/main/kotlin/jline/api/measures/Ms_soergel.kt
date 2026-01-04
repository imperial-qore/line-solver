package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max

/**
 * Soergel distance between two probability distributions.
 * Part of the L1 family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Soergel distance
 */
fun ms_soergel(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sumDiff = 0.0
    var sumMax = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        sumDiff += abs(pi - qi)
        sumMax += max(pi, qi)
    }
    
    return if (sumMax > 0) sumDiff / sumMax else 0.0
}
/**
 * Soergel metric algorithms
 */
@Suppress("unused")
class MsSoergelAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}