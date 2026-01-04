package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.min

/**
 * Kulczynski s distance between two probability distributions.
 * Part of the intersection family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Kulczynski s distance
 */
fun ms_kulczynskis(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sumMin = 0.0
    var sumDiff = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        sumMin += min(pi, qi)
        sumDiff += abs(pi - qi)
    }
    
    return if (sumMin > 0) sumDiff / sumMin else 0.0
}
/**
 * Kulczynskis metric algorithms
 */
@Suppress("unused")
class MsKulczynskisAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}