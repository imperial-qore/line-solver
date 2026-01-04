package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max

/**
 * Wave-Hedges distance between two probability distributions.
 * Part of the intersection family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Wave-Hedges distance
 */
fun ms_wavehegdes(P: Matrix, Q: Matrix): Double {
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
 * Wavehegdes metric algorithms
 */
@Suppress("unused")
class MsWavehegdesAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}