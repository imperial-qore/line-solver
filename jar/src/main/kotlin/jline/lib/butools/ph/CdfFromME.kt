/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph

import jline.util.matrix.Matrix

/**
 * Returns the cumulative distribution function of a
 * matrix-exponential distribution.
 *
 * @param alpha The initial vector of the matrix-exponential distribution.
 * @param A The matrix parameter of the matrix-exponential distribution.
 * @param x The points at which the CDF will be computed.
 * @return The values of the CDF at the corresponding "x" values.
 */
fun cdfFromME(alpha: Matrix, A: Matrix, x: DoubleArray): DoubleArray {
    val cdf = DoubleArray(x.size)
    for (i in x.indices) {
        if (x[i] < 0) {
            cdf[i] = 0.0
        } else {
            val expAt = (A.scale(x[i])).expm()
            cdf[i] = 1.0 - alpha.mult(expAt).elementSum()
        }
    }
    return cdf
}

/**
 * Returns the cumulative distribution function of a
 * phase-type distribution.
 */
fun cdfFromPH(alpha: Matrix, A: Matrix, x: DoubleArray): DoubleArray {
    return cdfFromME(alpha, A, x)
}

/**
 * Overload for DoubleArray alpha.
 */
fun cdfFromME(alpha: DoubleArray, A: Matrix, x: DoubleArray): DoubleArray {
    return cdfFromME(Matrix(alpha), A, x)
}

/**
 * Overload for DoubleArray alpha.
 */
fun cdfFromPH(alpha: DoubleArray, A: Matrix, x: DoubleArray): DoubleArray {
    return cdfFromME(Matrix(alpha), A, x)
}
