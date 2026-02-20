package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Returns the factorial moments given the raw moments.
 *
 * The raw moments are: m_i = E(X^i)
 * The factorial moments are: f_i = E(X(X-1)...(X-i+1))
 *
 * @param m The list of raw moments (starting with the first moment)
 * @return The list of factorial moments
 *
 * Reference: http://en.wikipedia.org/wiki/Factorial_moment
 */
fun factorialMomsFromMoms(m: Matrix): Matrix {
    val n = m.length()
    val fm = Matrix(m.numRows, m.numCols, m.length())

    for (i in 0 until n) {
        // Build roots matrix [0, 1, ..., i]
        val rootsMatrix = Matrix(1, i + 1, i + 1)
        for (k in 0..i) {
            rootsMatrix[k] = k.toDouble()
        }

        // poly([0,1,...,i]) gives coefficients of (x-0)(x-1)...(x-i)
        // which has i+2 coefficients [c_{i+1}, c_i, ..., c_1, c_0]
        val ehMatrix = poly(rootsMatrix)

        // MATLAB: eh = eh(end-1:-1:1)
        // Take from 2nd-to-last down to 1st element (1-based),
        // i.e., from index i down to 0 (0-based), giving i+1 elements
        val eh = DoubleArray(i + 1) { k -> ehMatrix.get(0, i - k) }

        // dot product: eh * m(1:i+1)
        var sum = 0.0
        for (k in 0..i) {
            sum += eh[k] * m[k]
        }
        fm[i] = sum
    }

    return fm
}
