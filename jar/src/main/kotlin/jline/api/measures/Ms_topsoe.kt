package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.ln

/**
 * Topsoe distance between two probability distributions.
 * Part of Shannon's entropy family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Topsoe distance
 */
fun ms_topsoe(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val denom = pi + qi
        
        if (pi > 0 && denom > 0) {
            sum += pi * ln(2.0 * pi / denom)
        }
        
        if (qi > 0 && denom > 0) {
            sum += qi * ln(2.0 * qi / denom)
        }
    }
    
    return sum
}
/**
 * Topsoe metric algorithms
 */
@Suppress("unused")
class MsTopsoeAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}