package jline.api.measures

import jline.util.matrix.Matrix

/**
 * Product distance between two probability distributions.
 * Part of the inner product family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Product distance
 */
fun ms_product(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        sum += P.get(i) * Q.get(i)
    }
    
    return 1.0 - sum
}
/**
 * Product metric algorithms
 */
@Suppress("unused")
class MsProductAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}