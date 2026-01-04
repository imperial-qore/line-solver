/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.util.matrix.Matrix

/**
 * Returns the cumulative distribution function of a matrix-geometric distribution.
 *
 * @param alpha The initial vector of the matrix-geometric distribution.
 * @param A The matrix parameter of the matrix-geometric distribution.
 * @param x Vector of non-negative integers at which to compute the CDF.
 * @return The probabilities that the matrix-geometrically distributed random variable
 *         is less or equal to the corresponding "x" values.
 */
fun cdfFromMG(alpha: Matrix, A: Matrix, x: IntArray): DoubleArray {
    val cdf = DoubleArray(x.size)

    for (i in x.indices) {
        // CDF(x) = 1 - alpha * A^x * ones
        val APower = Matrix.pow(A, x[i])
        val term = alpha.mult(APower)
        cdf[i] = 1.0 - term.elementSum()
    }

    return cdf
}

/**
 * Overload for DoubleArray alpha.
 */
fun cdfFromMG(alpha: DoubleArray, A: Matrix, x: IntArray): DoubleArray {
    return cdfFromMG(Matrix(alpha), A, x)
}

/**
 * Overload for single integer x.
 */
fun cdfFromMG(alpha: Matrix, A: Matrix, x: Int): Double {
    return cdfFromMG(alpha, A, intArrayOf(x))[0]
}

/**
 * Overload for DoubleArray alpha and single integer x.
 */
fun cdfFromMG(alpha: DoubleArray, A: Matrix, x: Int): Double {
    return cdfFromMG(Matrix(alpha), A, intArrayOf(x))[0]
}
