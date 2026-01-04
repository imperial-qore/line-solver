/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.dph.momentsFromDPH
import jline.util.matrix.Matrix

/**
 * Returns the moments of the marginal distribution of a discrete
 * Markovian arrival process.
 *
 * @param D0 The D0 matrix of the discrete Markovian arrival process
 * @param D1 The D1 matrix of the discrete Markovian arrival process
 * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
 * @param prec Numerical precision for checking if the input is valid
 * @return The vector of moments
 */
@JvmOverloads
fun marginalMomentsFromDMAP(D0: Matrix, D1: Matrix, K: Int = 0, prec: Double = 1e-14): DoubleArray {
    val numMoments = if (K == 0) 2 * D0.numRows - 1 else K

    if (!checkDMAPRepresentation(D0, D1, prec)) {
        throw IllegalArgumentException("MarginalMomentsFromDMAP: Input isn't a valid DMAP representation!")
    }

    val margDist = marginalDistributionFromDMAP(D0, D1, prec)
    return momentsFromDPH(margDist.alpha, margDist.A, numMoments)
}
