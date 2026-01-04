/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc

import jline.util.matrix.Matrix
import kotlin.math.abs

/**
 * Computes the stationary solution of a continuous time
 * rational process (CRP).
 *
 * @param Q The generator matrix of the rational process
 * @param prec Numerical precision. The default value is 1e-14.
 * @return The vector that satisfies pi*Q = 0, sum(pi) = 1
 *
 * Note: Continuous time rational processes are like continuous
 * time Markov chains, but the generator does not have to pass
 * the checkGenerator test (but the rowsums still have to be zeros).
 */
fun crpSolve(Q: Matrix, prec: Double = 1e-14): Matrix {
    val n = Q.numRows

    // Check rowsums
    for (i in 0 until n) {
        var rowSum = 0.0
        for (j in 0 until n) {
            rowSum += Q[i, j]
        }
        if (abs(rowSum) > prec) {
            throw IllegalArgumentException("CRPSolve: The matrix has a rowsum which isn't zero!")
        }
    }

    // M = Q with first column replaced by ones
    val M = Q.copy()
    for (i in 0 until n) {
        M[i, 0] = 1.0
    }

    // m = [1, 0, 0, ...]
    val m = Matrix.zeros(1, n)
    m[0, 0] = 1.0

    // pi = m / M = m * inv(M)
    val pi = m.mult(M.inv())

    return pi
}

/**
 * Computes the stationary solution of a discrete time
 * rational process (DRP).
 *
 * @param P The matrix parameter of the rational process
 * @return The vector that satisfies pi*P = pi, sum(pi) = 1
 *
 * Note: Discrete time rational processes are like discrete
 * time Markov chains, but the P matrix does not have to pass
 * the checkProbMatrix test (but the rowsums still have to be ones).
 */
fun drpSolve(P: Matrix): Matrix {
    val n = P.numRows
    // Q = P - I
    val Q = P.sub(Matrix.eye(n))
    return crpSolve(Q)
}
