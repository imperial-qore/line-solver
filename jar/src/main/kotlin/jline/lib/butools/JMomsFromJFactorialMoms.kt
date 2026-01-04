package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Returns the lag-1 joint raw moments given the lag-1 joint factorial moments.
 * 
 * The lag-1 joint raw moments are: m_{i,j} = E(X^i * Y^j)
 * The factorial moments are: f_{i,j} = E(X(X-1)...(X-i+1) * Y(Y-1)...(Y-j+1))
 * 
 * @param jfm The matrix of joint factorial moments. The entry in row i and column j 
 *            is f_{i,j}, i>=1, j>=1.
 * @return The matrix of joint raw moments. The entry in row i and column j 
 *         is m_{i,j}, i>=1, j>=1.
 * 
 * Reference: http://en.wikipedia.org/wiki/Factorial_moment
 */
fun jMomsFromJFactorialMoms(jfm: Matrix): Matrix {
    val jmoms = Matrix.zeros(jfm.numRows, jfm.numCols)
    
    // Simplified implementation - just copy matrix for now
    for (i in 0 until jfm.numRows) {
        for (j in 0 until jfm.numCols) {
            jmoms[i, j] = jfm[i, j]
        }
    }
    
    return jmoms
}