package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Returns the raw moments given the factorial moments.
 *
 * The raw moments are: m_i = E(X^i)
 * The factorial moments are: f_i = E(X(X-1)...(X-i+1))
 *
 * @param fm The list of factorial moments (starting with the first moment)
 * @return The list of raw moments
 *
 * Reference: http://en.wikipedia.org/wiki/Factorial_moment
 */
fun MomsFromFactorialMoms(fm: Matrix): Matrix {
    val n = fm.length()
    val m = Matrix(1, n, n)
    m[0] = fm[0]

    for (i in 1 until n) {
        // Compute polynomial with roots [0, 1, 2, ..., i]
        val rootsMatrix = Matrix(1, i + 1, i + 1)
        for (k in 0..i) {
            rootsMatrix[k] = k.toDouble()
        }

        // Get polynomial coefficients and negate
        // poly([0,1,...,i]) gives i+2 coefficients
        val ehMatrix = poly(rootsMatrix)
        val eh = DoubleArray(ehMatrix.numCols) { k -> -ehMatrix.get(0, k) }

        // MATLAB: eh(end-1:-1:2) for array of length i+2
        // 1-based: from index i+1 down to 2 = i elements
        // 0-based: from index i down to 1 = i elements
        val eh_coeff = DoubleArray(i) { k -> eh[i - k] }

        // Compute dot product: eh_coeff[0..i-1] dot m[0..i-1]
        var sum = 0.0
        for (k in 0 until i) {
            sum += eh_coeff[k] * m[k]
        }

        m[i] = fm[i] + sum
    }

    if (fm.numRows > fm.numCols) {
        return m.transpose()
    }
    return m
}
