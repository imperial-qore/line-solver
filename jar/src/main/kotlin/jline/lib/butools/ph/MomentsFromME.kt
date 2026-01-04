/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph

import jline.util.matrix.Matrix

/**
 * Returns the first K moments of a matrix-exponential
 * distribution.
 *
 * @param alpha The initial vector of the matrix-exponential distribution.
 * @param A The matrix parameter of the matrix-exponential distribution.
 * @param K Number of moments to compute. If K=0, 2*M-1 moments
 *        are computed. The default value is 0.
 * @return The vector of moments.
 */
fun momentsFromME(alpha: Matrix, A: Matrix, K: Int = 0): DoubleArray {
    val numMoments = if (K == 0) 2 * alpha.length() - 1 else K

    val moms = DoubleArray(numMoments)
    val iA = A.neg().inv()

    var iAPower = iA.copy()
    for (i in 0 until numMoments) {
        val factorial = (1..i + 1).fold(1L) { acc, n -> acc * n }
        moms[i] = factorial.toDouble() * alpha.mult(iAPower).elementSum()
        iAPower = iAPower.mult(iA)
    }

    return moms
}

/**
 * Returns the first K moments of a phase-type distribution.
 * This is equivalent to momentsFromME for PH distributions.
 *
 * @param alpha The initial vector of the phase-type distribution.
 * @param A The transient generator of the phase-type distribution.
 * @param K Number of moments to compute. If K=0, 2*M-1 moments
 *        are computed. The default value is 0.
 * @return The vector of moments.
 */
fun momentsFromPH(alpha: Matrix, A: Matrix, K: Int = 0): DoubleArray {
    return momentsFromME(alpha, A, K)
}

/**
 * Overload for DoubleArray alpha.
 */
fun momentsFromME(alpha: DoubleArray, A: Matrix, K: Int = 0): DoubleArray {
    return momentsFromME(Matrix(alpha), A, K)
}

/**
 * Overload for DoubleArray alpha.
 */
fun momentsFromPH(alpha: DoubleArray, A: Matrix, K: Int = 0): DoubleArray {
    return momentsFromME(Matrix(alpha), A, K)
}
