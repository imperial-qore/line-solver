package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.ln

/**
 * Lorentzian distance between two probability distributions.
 * Part of the L1 family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Lorentzian distance
 */
fun ms_lorentzian(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        sum += ln(1.0 + abs(P.get(i) - Q.get(i)))
    }
    
    return sum
}
/**
 * Lorentzian metric algorithms
 */
@Suppress("unused")
class MsLorentzianAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}