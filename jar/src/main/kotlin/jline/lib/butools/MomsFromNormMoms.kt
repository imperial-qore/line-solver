package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Returns the raw moments given the normalized moments.
 * 
 * The raw moments are: m_i = E(X^i)
 * The normalized moments are: n_i = m_i / (m_{i-1} * m_1)
 * 
 * @param nm The list of normalized moments (starting with the first moment)
 * @return The list of raw moments
 */
fun momsFromNormMoms(nm: Matrix): Matrix {
    val m = Matrix(nm.numRows, nm.numCols, nm.numRows * nm.numCols)
    
    for (i in 0 until nm.length()) {
        if (i == 0) {
            m[i] = nm[i]
        } else {
            m[i] = m[0] * nm[i] * m[i - 1]
        }
    }
    
    return m
}