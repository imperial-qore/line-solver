package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.sqrt

/**
 * Fidelity distance between two probability distributions.
 * Part of the fidelity family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Fidelity distance
 */
fun ms_fidelity(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        sum += sqrt(P.get(i) * Q.get(i))
    }
    
    return 1.0 - sum
}
/**
 * Fidelity metric algorithms
 */
@Suppress("unused")
class MsFidelityAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}