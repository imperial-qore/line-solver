/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc

import jline.util.matrix.Matrix

/**
 * Computes the stationary solution of a discrete time
 * Markov chain.
 *
 * @param P The transition probability matrix of the Markov chain.
 * @param prec Numerical precision. The default value is 1e-14.
 * @return The stationary probability vector (row vector).
 */
@JvmOverloads
fun dtmcSolve(P: Matrix, prec: Double = 1e-14): Matrix {
    val n = P.numRows

    // Compute Q = P - I
    val Q = P.sub(Matrix.eye(n))

    // Use CTMC solve on the embedded generator
    return ctmcSolve(Q, prec)
}
