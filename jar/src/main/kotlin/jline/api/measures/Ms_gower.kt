package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Gower distance between two probability distributions.
 * Part of the L1 family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Gower distance
 */
fun ms_gower(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        sum += abs(P.get(i) - Q.get(i))
    }
    
    return sum / P.length()
}
/**
 * Gower metric algorithms
 */
@Suppress("unused")
class MsGowerAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}