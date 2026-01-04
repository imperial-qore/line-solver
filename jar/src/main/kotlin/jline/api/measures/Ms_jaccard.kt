package jline.api.measures

import jline.util.matrix.Matrix

/**
 * Jaccard distance between two probability distributions.
 * Part of the inner product family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Jaccard distance
 */
fun ms_jaccard(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var dotProduct = 0.0
    var sumP2 = 0.0
    var sumQ2 = 0.0
    
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        dotProduct += pi * qi
        sumP2 += pi * pi
        sumQ2 += qi * qi
    }
    
    val s = dotProduct / (sumP2 + sumQ2 - dotProduct)
    return 1.0 - s
}
/**
 * Jaccard metric algorithms
 */
@Suppress("unused")
class MsJaccardAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}