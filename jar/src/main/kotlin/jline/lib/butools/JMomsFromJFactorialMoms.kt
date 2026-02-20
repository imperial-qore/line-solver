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
    val s1 = jfm.numRows
    val s2 = jfm.numCols
    val jmoms = Matrix.zeros(s1, s2)

    for (i in 0 until s1) {
        for (j in 0 until s2) {
            // xCoeff = poly(0:i) coefficients reversed (skip last)
            val xRoots = Matrix(1, i + 1, i + 1)
            for (k in 0..i) xRoots[k] = k.toDouble()
            val xPoly = poly(xRoots) // i+2 coefficients
            val xCoeff = DoubleArray(i + 1) { k -> xPoly.get(0, i - k) }

            val yRoots = Matrix(1, j + 1, j + 1)
            for (k in 0..j) yRoots[k] = k.toDouble()
            val yPoly = poly(yRoots) // j+2 coefficients
            val yCoeff = DoubleArray(j + 1) { k -> yPoly.get(0, j - k) }

            // eh = -(xCoeff' * yCoeff) (negated outer product)
            val eh = Matrix.zeros(i + 1, j + 1)
            for (r in 0..i) {
                for (c in 0..j) {
                    eh[r, c] = -xCoeff[r] * yCoeff[c]
                }
            }

            // jmoms(i,j) = jfmoms(i,j) + trace(jmoms(1:i, 1:j) * eh')
            // trace(A * B^T) = sum(A .* B) element-wise
            var traceVal = 0.0
            for (r in 0..i) {
                for (k in 0..j) {
                    traceVal += jmoms[r, k] * eh[r, k]
                }
            }
            jmoms[i, j] = jfm[i, j] + traceVal
        }
    }

    return jmoms
}
