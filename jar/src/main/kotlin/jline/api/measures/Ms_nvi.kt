/**
 * @file Normalized variation of information for discrete variables
 * 
 * Computes normalized variation information (1-I(X,Y)/H(X,Y)) measuring the normalized
 * amount of information lost in going from the joint distribution to individual variables.
 * Used as a distance metric in clustering and classification validation.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.max

/**
 * Compute normalized variation information z=(1-I(x,y)/H(x,y)) of two discrete variables x and y.
 * 
 * @param x first matrix
 * @param y second matrix of the same length
 * @return normalized variation information z=(1-I(x,y)/H(x,y))
 */
fun ms_nvi(x: Matrix, y: Matrix): Double {
    require(x.length() == y.length()) { "Input matrices must have the same length" }
    
    if (x.isEmpty) return 0.0
    
    val hx = ms_entropy(x)
    val hy = ms_entropy(y)
    val hxy = ms_jointentropy(x, y)
    
    // Handle edge case where joint entropy is zero
    if (hxy == 0.0) return 0.0
    
    // ms_nvi = 2 - (H(x) + H(y)) / H(x,y)
    val ms_nvi = 2.0 - (hx + hy) / hxy
    return max(0.0, ms_nvi)
}
/**
 * Ms_nvi algorithms
 */
@Suppress("unused")
class Ms_nvi {
    companion object {
        // Class documentation marker for Dokka
    }
}