/**
 * @file City block (Manhattan) distance metric
 * 
 * Implements the City block distance d(p,q) = ∑|pᵢ-qᵢ| between probability 
 * distributions. Also known as Manhattan distance or L¹ norm, commonly used 
 * in clustering and nearest neighbor applications for robust distance measurement.
 * 
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * City block (Manhattan) distance between two probability distributions.
 * Part of the Minkowski family. Also known as L1 distance.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return City block distance
 */
fun ms_cityblock(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        sum += abs(P.get(i) - Q.get(i))
    }
    
    return sum
}
/**
 * Cityblock metric algorithms
 */
@Suppress("unused")
class MsCityblockAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}