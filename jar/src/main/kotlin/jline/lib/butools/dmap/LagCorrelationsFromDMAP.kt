/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix

/**
 * Returns the lag autocorrelations of a discrete Markovian arrival process.
 *
 * @param D0 The D0 matrix of the discrete Markovian arrival process
 * @param D1 The D1 matrix of the discrete Markovian arrival process
 * @param L The number of lags to compute. The default value is 1.
 * @param prec Numerical precision to check if the input is valid
 * @return The lag autocorrelation function up to lag L
 */
fun lagCorrelationsFromDMAP(D0: Matrix, D1: Matrix, L: Int = 1, prec: Double = 1e-14): DoubleArray {
    if (!checkDMAPRepresentation(D0, D1, prec)) {
        throw IllegalArgumentException("LagCorrelationsFromDMAP: Input isn't a valid DMAP representation!")
    }

    return lagCorrelationsFromDRAP(D0, D1, L, prec)
}
