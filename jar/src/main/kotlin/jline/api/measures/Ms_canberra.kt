package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Canberra distance between two probability distributions.
 * Part of the L1 family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Canberra distance
 */
fun ms_canberra(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val denom = pi + qi
        if (denom > 0) {
            sum += abs(pi - qi) / denom
        }
    }
    
    return sum
}
/**
 * Canberra metric algorithms
 */
@Suppress("unused")
class MsCanberraAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}