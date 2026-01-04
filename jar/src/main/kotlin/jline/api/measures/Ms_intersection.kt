package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Intersection distance between two probability distributions.
 * Part of the intersection family.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Intersection distance
 */
fun ms_intersection(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var sum = 0.0
    for (i in 0 until P.length()) {
        sum += abs(P.get(i) - Q.get(i))
    }
    
    return 0.5 * sum
}
/**
 * Intersection metric algorithms
 */
@Suppress("unused")
class MsIntersectionAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}