/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.dph.momentsFromMG
import jline.util.matrix.Matrix

/**
 * Returns the moments of the marginal distribution of a discrete
 * rational arrival process.
 *
 * @param H0 The H0 matrix of the discrete rational arrival process
 * @param H1 The H1 matrix of the discrete rational arrival process
 * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
 * @param prec Numerical precision for checking if the input is valid
 * @return The vector of moments
 */
fun marginalMomentsFromDRAP(H0: Matrix, H1: Matrix, K: Int = 0, prec: Double = 1e-14): DoubleArray {
    val numMoments = if (K == 0) 2 * H0.numRows - 1 else K

    if (!checkDRAPRepresentation(H0, H1, prec)) {
        throw IllegalArgumentException("MarginalMomentsFromDRAP: Input isn't a valid DRAP representation!")
    }

    val margDist = marginalDistributionFromDRAP(H0, H1, prec)
    return momentsFromMG(margDist.alpha, margDist.A, numMoments)
}
