/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.util.matrix.Matrix

/**
 * Returns the probability mass function of a matrix-geometric distribution.
 *
 * @param alpha The initial vector of the matrix-geometric distribution.
 *        The sum of the entries of alpha is less or equal to 1.
 * @param A The matrix parameter of the matrix-geometric distribution.
 * @param x Vector of non-negative integers at which to compute the PMF.
 * @return The probabilities that the matrix-geometrically distributed random variable
 *         takes the corresponding "x" values.
 */
fun pmfFromMG(alpha: Matrix, A: Matrix, x: IntArray): DoubleArray {
    val N = A.numRows

    // a = 1 - sum(A, 2) = closing vector
    val a = Matrix(N, 1)
    for (i in 0 until N) {
        var rowSum = 0.0
        for (j in 0 until N) {
            rowSum += A[i, j]
        }
        a[i, 0] = 1.0 - rowSum
    }

    val pmf = DoubleArray(x.size)

    for (i in x.indices) {
        if (x[i] == 0) {
            // PMF(0) = 1 - sum(alpha)
            pmf[i] = 1.0 - alpha.elementSum()
        } else {
            // PMF(x) = alpha * A^(x-1) * a
            val APower = Matrix.pow(A, x[i] - 1)
            val term = alpha.mult(APower).mult(a)
            pmf[i] = term.elementSum()
        }
    }

    return pmf
}

/**
 * Overload for DoubleArray alpha.
 */
fun pmfFromMG(alpha: DoubleArray, A: Matrix, x: IntArray): DoubleArray {
    return pmfFromMG(Matrix(alpha), A, x)
}

/**
 * Overload for single integer x.
 */
fun pmfFromMG(alpha: Matrix, A: Matrix, x: Int): Double {
    return pmfFromMG(alpha, A, intArrayOf(x))[0]
}

/**
 * Overload for DoubleArray alpha and single integer x.
 */
fun pmfFromMG(alpha: DoubleArray, A: Matrix, x: Int): Double {
    return pmfFromMG(Matrix(alpha), A, intArrayOf(x))[0]
}
