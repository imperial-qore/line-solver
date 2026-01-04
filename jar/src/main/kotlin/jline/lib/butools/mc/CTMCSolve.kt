/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc

import jline.util.matrix.Matrix

/**
 * Computes the stationary solution of a continuous time
 * Markov chain.
 *
 * @param Q The generator matrix of the Markov chain.
 * @param prec Numerical precision. The default value is 1e-14.
 * @return The stationary probability vector (row vector).
 */
@JvmOverloads
fun ctmcSolve(Q: Matrix, prec: Double = 1e-14): Matrix {
    val n = Q.numRows

    // Create a copy of Q and modify the last column
    val A = Q.copy()
    for (i in 0 until n) {
        A[i, n - 1] = 1.0
    }

    // Create right-hand side vector
    val b = Matrix.zeros(1, n)
    b[0, n - 1] = 1.0

    // Solve the linear system b = pi * A, i.e., pi = b * inv(A)
    val pi = b.mult(A.inv())

    // Normalize to ensure sum is 1
    val sum = pi.elementSum()
    if (sum > 0) {
        for (i in 0 until n) {
            pi[0, i] = pi[0, i] / sum
        }
    }

    return pi
}
