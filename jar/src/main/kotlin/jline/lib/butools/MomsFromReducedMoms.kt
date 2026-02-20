package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Returns the raw moments given the reduced moments.
 * 
 * The raw moments are: m_i = E(X^i)
 * The reduced moments are: r_i = m_i / i!
 * 
 * @param rm The list of reduced moments (starting with the first moment)
 * @return The list of raw moments
 */
fun momsFromReducedMoms(rm: Matrix): Matrix {
    val m = Matrix(rm.numRows, rm.numCols, rm.numRows * rm.numCols)
    var factorial = 1.0
    
    for (i in 0 until rm.length()) {
        factorial *= (i + 1) // Calculate (i+1)!
        m[i] = rm[i] * factorial
    }
    
    return m
}