package jline.lib.butools

import jline.util.matrix.Matrix

/**
 * Returns the lag-1 joint factorial moments given the lag-1 joint raw moments.
 *
 * The lag-1 joint raw moments are: m_{i,j} = E(X^i * Y^j)
 * The factorial moments are: f_{i,j} = E(X(X-1)...(X-i+1) * Y(Y-1)...(Y-j+1))
 *
 * @param jm The matrix of joint raw moments. The entry in row i and column j
 *           is m_{i,j}, i>=1, j>=1.
 * @return The matrix of joint factorial moments. The entry in row i and column j
 *         is f_{i,j}, i>=1, j>=1.
 *
 * Reference: http://en.wikipedia.org/wiki/Factorial_moment
 */
fun jFactorialMomsFromJMoms(jm: Matrix): Matrix {
    val s1 = jm.numRows
    val s2 = jm.numCols
    val jfmoms = Matrix.zeros(s1, s2)

    for (i in 0 until s1) {
        for (j in 0 until s2) {
            // xCoeff = poly(0:i) coefficients reversed (skip last)
            // MATLAB: xCoeff = poly(0:i-1); xCoeff = xCoeff(end-1:-1:1);
            // 1-based i, roots 0..i-1 → i+1 coefficients
            // Here 0-based i corresponds to MATLAB i+1, so roots 0..i → i+1 values
            val xRoots = Matrix(1, i + 1, i + 1)
            for (k in 0..i) xRoots[k] = k.toDouble()
            val xPoly = poly(xRoots) // i+2 coefficients
            // MATLAB: xCoeff = xCoeff(end-1:-1:1) → 0-based: indices i down to 0 → i+1 elements
            val xCoeff = DoubleArray(i + 1) { k -> xPoly.get(0, i - k) }

            val yRoots = Matrix(1, j + 1, j + 1)
            for (k in 0..j) yRoots[k] = k.toDouble()
            val yPoly = poly(yRoots) // j+2 coefficients
            val yCoeff = DoubleArray(j + 1) { k -> yPoly.get(0, j - k) }

            // eh = xCoeff' * yCoeff (outer product, (i+1) x (j+1) matrix)
            val eh = Matrix.zeros(i + 1, j + 1)
            for (r in 0..i) {
                for (c in 0..j) {
                    eh[r, c] = xCoeff[r] * yCoeff[c]
                }
            }

            // jfmoms(i,j) = trace(jm(1:i, 1:j) * eh')
            // 0-based: trace(jm[0..i, 0..j] * eh^T)
            // eh^T is (j+1) x (i+1), jm sub is (i+1) x (j+1)
            // product is (i+1) x (i+1), trace = sum of diagonal
            var traceVal = 0.0
            val sz = i + 1 // square trace size
            for (r in 0..i) {
                for (k in 0..j) {
                    traceVal += jm[r, k] * eh[r, k] // trace(A * B^T) = sum(A .* B)
                }
            }
            jfmoms[i, j] = traceVal
        }
    }

    return jfmoms
}
