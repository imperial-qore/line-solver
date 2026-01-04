/**
 * Chebyshev Distance for Probability Distributions
 * 
 * Implements the Chebyshev distance (L∞ norm) computing the maximum absolute difference
 * between corresponding elements. Also known as chess-board distance, it represents
 * the limiting case of Minkowski distance as p approaches infinity.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max

/**
 * Chebyshev distance between two probability distributions.
 * Part of the Minkowski family. Also known as L-infinity distance.
 * 
 * @param P first probability distribution
 * @param Q second probability distribution
 * @return Chebyshev distance
 */
fun ms_chebyshev(P: Matrix, Q: Matrix): Double {
    require(P.length() == Q.length()) { "Distributions must have the same size" }
    
    var maxDist = 0.0
    for (i in 0 until P.length()) {
        maxDist = max(maxDist, abs(P.get(i) - Q.get(i)))
    }
    
    return maxDist
}
/**
 * Chebyshev metric algorithms
 */
@Suppress("unused")
class MsChebyshevAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}