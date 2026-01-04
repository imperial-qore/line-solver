package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Kumar-Johnson distance between two probability distributions.
 * Part of the combination family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Kumar-Johnson distance
 */
fun ms_kumarjohnson(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val denom = 2.0 * (pi * qi).pow(3.0/2.0)
        if (denom > 0) {
            val diff = pi * pi - qi * qi
            sum += (diff * diff) / denom
        }
    }
    
    return sum
}
/**
 * Kumarjohnson metric algorithms
 */
@Suppress("unused")
class MsKumarjohnsonAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}