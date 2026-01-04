package jline.api.measures

import jline.util.matrix.Matrix

/**
 * Harmonic mean distance between two probability distributions.
 * Part of the inner product family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Harmonic mean distance
 */
fun ms_harmonicmean(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        val pi = P.get(i)
        val qi = Q.get(i)
        val denom = pi + qi
        if (denom > 0) {
            sum += pi * qi / denom
        }
    }
    
    return 1.0 - 2.0 * sum
}
/**
 * Harmonicmean metric algorithms
 */
@Suppress("unused")
class MsHarmonicmeanAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}