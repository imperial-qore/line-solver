package jline.api.measures

import jline.util.matrix.Matrix
import kotlin.math.max

/**
 * Compute conditional entropy z=H(x|y) of two discrete variables x and y.
 * 
 * @param x first matrix
 * @param y second matrix of the same length
 * @return conditional entropy z=H(x|y)
 */
fun ms_condentropy(x: Matrix, y: Matrix): Double {
    require(x.length() == y.length()) { "Input matrices must have the same length" }
    
    if (x.isEmpty) return 0.0
    
    // H(x|y) = H(x,y) - H(y)
    val hxy = ms_jointentropy(x, y)
    val hy = ms_entropy(y)
    
    val condEntropy = hxy - hy
    return max(0.0, condEntropy)
}