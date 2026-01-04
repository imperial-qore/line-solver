package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.ln

/**
 * Jensen difference divergence between two probability distributions.
 * Part of Shannon's entropy family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Jensen difference
 */
fun ms_jensendifference(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val m = (pi + qi) / 2.0
        
        var term = 0.0
        if (pi > 0) term += pi * ln(pi)
        if (qi > 0) term += qi * ln(qi)
        term /= 2.0
        
        if (m > 0) term -= m * ln(m)
        
        sum += term
    }
    
    return sum
}
/**
 * Jensendifference metric algorithms
 */
@Suppress("unused")
class MsJensendifferenceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}