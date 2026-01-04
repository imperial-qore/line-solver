/**
 * @file Mutual information for discrete variables
 * 
 * Computes mutual information I(X,Y) = H(X) + H(Y) - H(X,Y) measuring the amount of 
 * information obtained about one variable through observing the other. Fundamental
 * measure in feature selection, clustering validation, and dependency analysis.
 *
 * @since LINE 3.0
 */
package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.max

/**
 * Compute mutual information I(x,y) of two discrete variables x and y.
 * 
 * @param x first matrix
 * @param y second matrix of the same length
 * @return mutual information z=I(x,y)
 */
fun ms_mutinfo(x: Matrix, y: Matrix): Double {
    require(x.length() == y.length()) { "Input matrices must have the same length" }
    
    if (x.isEmpty) return 0.0
    
    // I(x,y) = H(x) + H(y) - H(x,y)
    val hx = ms_entropy(x)
    val hy = ms_entropy(y)
    val hxy = ms_jointentropy(x, y)
    
    val mutualInfo = hx + hy - hxy
    return max(0.0, mutualInfo)
}