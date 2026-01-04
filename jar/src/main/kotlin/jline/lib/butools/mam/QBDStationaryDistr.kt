/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam

import jline.util.matrix.Matrix

/**
 * Returns the stationary distribution of a QBD up to a given level K.
 *
 * @param pi0 The stationary probability vector of level zero
 * @param R The matrix parameter of the matrix geometrical distribution of the QBD
 * @param K The stationary distribution is returned up to this level
 * @return The stationary probability vector up to level K (length (K+1)*N)
 */
fun qbdStationaryDistr(pi0: Matrix, R: Matrix, K: Int): Matrix {
    val m = R.numRows
    val pi = Matrix.zeros(1, (K + 1) * m)

    // Copy pi0 to first m elements
    for (i in 0 until m) {
        pi[0, i] = pi0[0, i]
    }

    var pix = pi0.copy()
    for (k in 1..K) {
        pix = pix.mult(R)
        for (i in 0 until m) {
            pi[0, k * m + i] = pix[0, i]
        }
    }

    return pi
}

/**
 * Overload for DoubleArray pi0.
 */
fun qbdStationaryDistr(pi0: DoubleArray, R: Matrix, K: Int): Matrix {
    return qbdStationaryDistr(Matrix(pi0), R, K)
}
