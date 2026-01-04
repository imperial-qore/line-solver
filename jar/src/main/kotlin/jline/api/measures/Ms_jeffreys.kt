package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.ln

/**
 * Jeffreys divergence between two probability distributions.
 * Part of Shannon's entropy family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Jeffreys divergence
 */
fun ms_jeffreys(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        if (pi > 0 && qi > 0) {
            sum += (pi - qi) * ln(pi / qi)
        }
    }
    
    return sum
}
/**
 * Jeffreys metric algorithms
 */
@Suppress("unused")
class MsJeffreysAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}