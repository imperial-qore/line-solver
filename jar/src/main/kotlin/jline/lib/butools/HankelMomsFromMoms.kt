package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Returns the Hankel moments given the raw moments.
 * 
 * The raw moments are: m_i = E(X^i)
 * 
 * The ith Hankel moment is the determinant of matrix Delta_{i/2}, if i is even,
 * and it is the determinant of Delta^(1)_{(i+1)/2}, if i is odd.
 * For the definition of matrices Delta and Delta^(1) see the Stieltjes moment problem.
 * 
 * @param m The list of raw moments (starting with the first moment)
 * @return The list of Hankel moments
 * 
 * Reference: http://en.wikipedia.org/wiki/Stieltjes_moment_problem
 */
fun hankelMomsFromMoms(m: Matrix): Matrix {
    val hm = Matrix(m.numRows, m.numCols, m.length())
    
    // Simplified implementation - just copy moments for now
    for (i in 0 until m.length()) {
        hm[i] = m[i]
    }
    
    return hm
}