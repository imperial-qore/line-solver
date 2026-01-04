/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph

import jline.util.matrix.Matrix

/**
 * Returns the first K moments of a discrete phase-type distribution.
 *
 * @param alpha The initial probability vector of the discrete phase-type distribution.
 *        The sum of the entries of alpha is less or equal to 1.
 * @param A The transient generator matrix of the discrete phase-type distribution.
 * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
 *        The default value is 0.
 * @return The vector of moments.
 */
@JvmOverloads
fun momentsFromDPH(alpha: Matrix, A: Matrix, K: Int = 0): DoubleArray {
    return momentsFromMG(alpha, A, K)
}

/**
 * Overload for DoubleArray alpha.
 */
@JvmOverloads
fun momentsFromDPH(alpha: DoubleArray, A: Matrix, K: Int = 0): DoubleArray {
    return momentsFromDPH(Matrix(alpha), A, K)
}
