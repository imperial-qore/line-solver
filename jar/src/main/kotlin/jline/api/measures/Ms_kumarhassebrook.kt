package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Kumar-Hassebrook distance between two probability distributions.
 * Part of the inner product family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Kumar-Hassebrook distance
 */
fun ms_kumarhassebrook(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var dotProduct = 0.0
    var sumDiff2 = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        dotProduct += pi * qi
        val diff = abs(pi - qi)
        sumDiff2 += diff * diff
    }
    
    val s = dotProduct / (sumDiff2 + dotProduct)
    return 1.0 - s
}
/**
 * Kumarhassebrook metric algorithms
 */
@Suppress("unused")
class MsKumarhassebrookAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}