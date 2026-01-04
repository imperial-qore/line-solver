/**
 * @file Normalized mutual information for discrete variables
 * 
 * Computes normalized mutual information I(X,Y)/√(H(X)×H(Y)) providing a scale-invariant
 * measure of dependence between variables. Commonly used in clustering evaluation and
 * feature selection where normalized measures are preferred over raw mutual information.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.max
import kotlin.math.sqrt

/**
 * Compute normalized mutual information I(x,y)/sqrt(H(x)*H(y)) of two discrete variables x and y.
 * 
 * @param x first matrix
 * @param y second matrix of the same length
 * @return normalized mutual information z=I(x,y)/sqrt(H(x)*H(y))
 */
fun ms_nmi(x: Matrix, y: Matrix): Double {
    require(x.length() == y.length()) { "Input matrices must have the same length" }
    
    if (x.isEmpty) return 0.0
    
    val hx = ms_entropy(x)
    val hy = ms_entropy(y)
    
    // Handle edge cases where entropy is zero
    if (hx == 0.0 || hy == 0.0) return 0.0
    
    val mi = ms_mutinfo(x, y)
    
    // Normalized mutual information
    val ms_nmi = sqrt((mi / hx) * (mi / hy))
    return max(0.0, ms_nmi)
}
/**
 * Ms_nmi algorithms
 */
@Suppress("unused")
class Ms_nmi {
    companion object {
        // Class documentation marker for Dokka
    }
}