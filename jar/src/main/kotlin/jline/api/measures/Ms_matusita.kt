package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.sqrt

/**
 * Matusita distance between two probability distributions.
 * Part of the fidelity family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Matusita distance
 */
fun ms_matusita(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val diff = sqrt(P.get(i)) - sqrt(Q.get(i))
        sum += diff * diff
    }
    
    return sqrt(sum)
}
/**
 * Matusita metric algorithms
 */
@Suppress("unused")
class MsMatusitaAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}