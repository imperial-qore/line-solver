package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max

/**
 * Average L1 L-infinity distance between two probability distributions.
 * Part of the combination family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Average L1 L-infinity distance
 */
fun ms_avgl1linfty(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sumL1 = 0.0
    var maxLinfty = 0.0
    
    for (i in 0 until P.length()) {
        val diff = abs(P.get(i) - Q.get(i))
        sumL1 += diff
        maxLinfty = max(maxLinfty, diff)
    }
    
    return (sumL1 + maxLinfty) / 2.0
}
/**
 * Avgl1Linfty metric algorithms
 */
@Suppress("unused")
class MsAvgl1linftyAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}