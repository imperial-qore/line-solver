/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.lib.butools.MomsFromFactorialMoms
import jline.util.matrix.Matrix

/**
 * Returns the first K moments of a matrix-geometric distribution.
 *
 * @param alpha The initial vector of the matrix-geometric distribution.
 *        The sum of the entries of alpha is less or equal to 1.
 * @param A The matrix parameter of the matrix-geometric distribution.
 * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
 *        The default value is 0.
 * @return The vector of moments.
 */
fun momentsFromMG(alpha: Matrix, A: Matrix, K: Int = 0): DoubleArray {
    val N = A.numRows
    
    var rowAlpha = alpha
    if (alpha.numRows == N && alpha.numCols == 1) {
        rowAlpha = alpha.transpose()
    }
    
    val numMoments = if (K == 0) 2 * N - 1 else K

    // iA = inv(I - A)
    val I = Matrix.eye(N)
    val iA = I.sub(A).inv()

    // Compute factorial moments: fmoms(i) = i! * alpha * iA^i * A^(i-1)
    val fmoms = DoubleArray(numMoments)
    var factorial = 1.0
    var iAPower = iA.copy()
    var APower = Matrix.eye(N)

    for (i in 1..numMoments) {
        factorial *= i
        // fmoms[i-1] = factorial * sum(alpha * iA^i * A^(i-1))
        val term = rowAlpha.mult(iAPower).mult(APower)
        fmoms[i - 1] = factorial * term.elementSum()

        // Update powers for next iteration
        iAPower = iAPower.mult(iA)
        APower = APower.mult(A)
    }

    // Convert factorial moments to raw moments
    val fmomsMatrix = Matrix(fmoms)
    val moms = MomsFromFactorialMoms(fmomsMatrix)

    val result = DoubleArray(numMoments)
    for (i in 0 until numMoments) {
        result[i] = moms[i]
    }
    return result
}

/**
 * Overload for DoubleArray alpha.
 */
fun momentsFromMG(alpha: DoubleArray, A: Matrix, K: Int = 0): DoubleArray {
    return momentsFromMG(Matrix(alpha), A, K)
}
