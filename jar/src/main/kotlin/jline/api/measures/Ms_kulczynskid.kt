package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.min

/**
 * Kulczynski d distance between two probability distributions.
 * Part of the L1 family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Kulczynski d distance
 */
fun ms_kulczynskid(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sumDiff = 0.0
    var sumMin = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        sumDiff += abs(pi - qi)
        sumMin += min(pi, qi)
    }
    
    return if (sumMin > 0) sumDiff / sumMin else 0.0
}
/**
 * Kulczynskid metric algorithms
 */
@Suppress("unused")
class MsKulczynskidAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}