package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.max
import kotlin.math.min

/**
 * Ruzicka distance between two probability distributions.
 * Part of the intersection family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Ruzicka distance
 */
fun ms_ruzicka(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sumMin = 0.0
    var sumMax = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        sumMin += min(pi, qi)
        sumMax += max(pi, qi)
    }
    
    return if (sumMax > 0) 1.0 - sumMin / sumMax else 0.0
}
/**
 * Ruzicka metric algorithms
 */
@Suppress("unused")
class MsRuzickaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}