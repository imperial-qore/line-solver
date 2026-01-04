/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.dph.momentsFromMG
import jline.lib.butools.mc.drpSolve
import jline.util.matrix.Matrix

/**
 * Returns the lag autocorrelations of a discrete rational arrival process.
 *
 * @param H0 The H0 matrix of the discrete rational arrival process
 * @param H1 The H1 matrix of the discrete rational arrival process
 * @param L The number of lags to compute. The default value is 1.
 * @param prec Numerical precision to check if the input is valid
 * @return The lag autocorrelation function up to lag L
 */
fun lagCorrelationsFromDRAP(H0: Matrix, H1: Matrix, L: Int = 1, prec: Double = 1e-14): DoubleArray {
    if (!checkDRAPRepresentation(H0, H1, prec)) {
        throw IllegalArgumentException("LagCorrelationsFromDRAP: Input isn't a valid DRAP representation!")
    }

    val n = H0.numRows
    val I = Matrix.eye(n)

    // H0i = inv(I - H0)
    val H0i = I.sub(H0).inv()

    // P = H0i * H1
    val P = H0i.mult(H1)

    // pi = DRPSolve(P)
    var pi = drpSolve(P)

    // moms = MomentsFromMG(pi, H0, 2)
    val moms = momentsFromMG(pi, H0, 2)

    // pi = pi * H0i * P
    pi = pi.mult(H0i).mult(P)

    // Compute lag correlations
    val acf = DoubleArray(L)
    for (i in 0 until L) {
        // acf(i) = (sum(pi*H0i) - moms(1)^2) / (moms(2) - moms(1)^2)
        val piH0i = pi.mult(H0i)
        val sumPiH0i = piH0i.elementSum()
        acf[i] = (sumPiH0i - moms[0] * moms[0]) / (moms[1] - moms[0] * moms[0])

        // pi = pi * P
        pi = pi.mult(P)
    }

    return acf
}
