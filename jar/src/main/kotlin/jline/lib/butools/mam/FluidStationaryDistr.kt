/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mam

import jline.util.matrix.Matrix

/**
 * Returns the stationary distribution of a Markovian fluid model
 * at specific points.
 *
 * @param mass0 The stationary probability vector of zero level
 * @param ini The initial vector of the stationary density
 * @param K The matrix parameter of the stationary density
 * @param clo The closing matrix of the stationary density
 * @param x The points at which the distribution is evaluated
 * @return The stationary distribution at the specified points
 */
fun fluidStationaryDistr(
    mass0: Matrix,
    ini: Matrix,
    K: Matrix,
    clo: Matrix,
    x: DoubleArray
): Matrix {
    val N = mass0.numCols
    val Np = K.numRows
    val result = Matrix.zeros(x.size, N)

    // closing = -K^{-1} * clo
    val closing = K.neg().inv().mult(clo)
    val I = Matrix.eye(Np)

    for (i in x.indices) {
        // y(x) = mass0 + ini * (I - e^(K*x)) * closing
        val expKx = K.scale(x[i]).expm()
        val IminusExp = I.sub(expKx)
        val contrib = ini.mult(IminusExp).mult(closing)
        for (j in 0 until N) {
            result[i, j] = mass0[0, j] + contrib[0, j]
        }
    }

    return result
}
