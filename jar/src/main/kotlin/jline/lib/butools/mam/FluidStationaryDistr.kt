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
    val result = Matrix.zeros(x.size, N)

    for (i in x.indices) {
        if (x[i] == 0.0) {
            // At zero, return mass0
            for (j in 0 until N) {
                result[i, j] = mass0[0, j]
            }
        } else {
            // For x > 0, compute ini * e^(K*x) * clo
            val expKx = K.scale(x[i]).expm()
            val density = ini.mult(expKx).mult(clo)
            for (j in 0 until N) {
                result[i, j] = density[0, j]
            }
        }
    }

    return result
}
