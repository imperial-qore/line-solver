/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map

import jline.lib.butools.mc.drpSolve
import jline.lib.butools.ph.momentsFromME
import jline.util.matrix.Matrix

/**
 * Returns the lag autocorrelations of a Markovian arrival process.
 *
 * @param D0 The D0 matrix of the Markovian arrival process
 * @param D1 The D1 matrix of the Markovian arrival process
 * @param L The number of lags to compute. The default value is 1
 * @return The lag autocorrelation function up to lag L
 */
@JvmOverloads
fun lagCorrelationsFromMAP(D0: Matrix, D1: Matrix, L: Int = 1): DoubleArray {
    return lagCorrelationsFromRAP(D0, D1, L)
}

/**
 * Returns the lag autocorrelations of a rational arrival process.
 *
 * @param H0 The H0 matrix of the rational arrival process
 * @param H1 The H1 matrix of the rational arrival process
 * @param L The number of lags to compute. The default value is 1
 * @return The lag autocorrelation function up to lag L
 */
@JvmOverloads
fun lagCorrelationsFromRAP(H0: Matrix, H1: Matrix, L: Int = 1): DoubleArray {
    val H0i = H0.neg().inv()
    val P = H0i.mult(H1)
    var pi = drpSolve(P)
    val moms = momentsFromME(pi, H0, 2)

    pi = pi.mult(H0i).mult(P)

    val acf = DoubleArray(L)
    for (i in 0 until L) {
        val sum = pi.mult(H0i).elementSum()
        acf[i] = (sum - moms[0] * moms[0]) / (moms[1] - moms[0] * moms[0])
        pi = pi.mult(P)
    }

    return acf
}
