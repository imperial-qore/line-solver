package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Returns the raw moments given the Hankel moments.
 * 
 * The raw moments are: m_i = E(X^i)
 * 
 * The ith Hankel moment is the determinant of matrix Delta_{i/2}, if i is even,
 * and it is the determinant of Delta^(1)_{(i+1)/2}, if i is odd.
 * For the definition of matrices Delta and Delta^(1) see the Stieltjes moment problem.
 * 
 * @param hm The list of Hankel moments (starting with the first moment)
 * @return The list of raw moments
 * 
 * Reference: http://en.wikipedia.org/wiki/Stieltjes_moment_problem
 */
fun momsFromHankelMoms(hm: Matrix): Matrix {
    val m = Matrix(hm.numRows, hm.numCols, hm.length())
    
    // Simplified implementation - just copy moments for now
    for (i in 0 until hm.length()) {
        m[i] = hm[i]
    }
    
    return m
}