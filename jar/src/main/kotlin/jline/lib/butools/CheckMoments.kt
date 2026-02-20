package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Checks if the given moment sequence is valid in the sense
 * that it belongs to a distribution with support (0,inf).
 * 
 * This procedure checks the determinant of Delta_n
 * and Delta_n^(1) according to the Stieltjes moment problem.
 * 
 * @param m The (raw) moments to check (starts with the first moment).
 *          Its length must be odd.
 * @param prec Entries with absolute value less than prec are 
 *             considered to be zeros. Default value is 1e-14.
 * @return The result of the check.
 * 
 * Reference: http://en.wikipedia.org/wiki/Stieltjes_moment_problem
 */
fun checkMoments(m: Matrix, prec: Double = 1e-14): Boolean {
    if (m.length() % 2 == 0) {
        throw IllegalArgumentException("CheckMoments: the number of moments must be odd!")
    }
    
    // Prepend 1.0 to moments (zeroth moment)
    val moments = Matrix(1, m.length() + 1, m.length() + 1)
    moments[0] = 1.0
    for (i in 0 until m.length()) {
        moments[i + 1] = m[i]
    }
    val N = (moments.length() / 2) - 1
    
    for (n in 0..N) {
        // Create Hankel matrices
        val H = hankelMatrix(moments, n + 1, 0, n, 2 * n)
        val H0 = hankelMatrix(moments, n + 1, 1, n + 1, 2 * n + 1)
        
        if (H.det() < -prec || H0.det() < -prec) {
            return false
        }
    }
    
    return true
}

/**
 * Creates a Hankel matrix from the given moment vector.
 * A Hankel matrix has constant anti-diagonals.
 */
private fun hankelMatrix(moments: Matrix, size: Int, firstRowStart: Int, lastColStart: Int, lastColEnd: Int): Matrix {
    val result = Matrix.zeros(size, size)
    
    for (i in 0 until size) {
        for (j in 0 until size) {
            val idx = i + j
            if (idx < size) {
                result[i, j] = moments[firstRowStart + idx]
            } else {
                result[i, j] = moments[lastColStart + idx - size + 1]
            }
        }
    }
    
    return result
}