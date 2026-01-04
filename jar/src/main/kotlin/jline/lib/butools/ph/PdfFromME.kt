/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph

import jline.util.matrix.Matrix

/**
 * Returns the probability density function of a matrix-
 * exponential distribution.
 *
 * @param alpha The initial vector of the matrix-exponential distribution.
 * @param A The matrix parameter of the matrix-exponential distribution.
 * @param x The points at which the density function will be computed.
 * @return The values of the density function at the corresponding "x" values.
 */
fun pdfFromME(alpha: Matrix, A: Matrix, x: DoubleArray): DoubleArray {
    val pdf = DoubleArray(x.size)
    val negA = A.neg()
    for (i in x.indices) {
        if (x[i] < 0) {
            pdf[i] = 0.0
        } else {
            val expAt = (A.scale(x[i])).expm()
            pdf[i] = alpha.mult(expAt).mult(negA).elementSum()
        }
    }
    return pdf
}

/**
 * Returns the probability density function of a phase-type distribution.
 */
fun pdfFromPH(alpha: Matrix, A: Matrix, x: DoubleArray): DoubleArray {
    return pdfFromME(alpha, A, x)
}

/**
 * Overload for DoubleArray alpha.
 */
fun pdfFromME(alpha: DoubleArray, A: Matrix, x: DoubleArray): DoubleArray {
    return pdfFromME(Matrix(alpha), A, x)
}

/**
 * Overload for DoubleArray alpha.
 */
fun pdfFromPH(alpha: DoubleArray, A: Matrix, x: DoubleArray): DoubleArray {
    return pdfFromME(Matrix(alpha), A, x)
}
