/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.dph.momentsFromMG
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Returns the moments of the marginal distribution of a discrete
 * marked rational arrival process.
 *
 * @param H The H0...HN matrices of the DMRAP (as MatrixCell)
 * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
 * @param prec Numerical precision for checking if the input is valid
 * @return The vector of moments
 */
fun marginalMomentsFromDMRAP(H: MatrixCell, K: Int = 0, prec: Double = 1e-14): DoubleArray {
    val numMoments = if (K == 0) 2 * H[0].numRows - 1 else K

    if (!checkDMRAPRepresentation(H, prec)) {
        throw IllegalArgumentException("MarginalMomentsFromDMRAP: Input isn't a valid DMRAP representation!")
    }

    val margDist = marginalDistributionFromDMRAP(H, prec)
    return momentsFromMG(margDist.alpha, margDist.A, numMoments)
}

/**
 * Overload for Array<Matrix>.
 */
fun marginalMomentsFromDMRAP(H: Array<Matrix>, K: Int = 0, prec: Double = 1e-14): DoubleArray {
    val numMoments = if (K == 0) 2 * H[0].numRows - 1 else K

    if (!checkDMRAPRepresentation(H, prec)) {
        throw IllegalArgumentException("MarginalMomentsFromDMRAP: Input isn't a valid DMRAP representation!")
    }

    val margDist = marginalDistributionFromDMRAP(H, prec)
    return momentsFromMG(margDist.alpha, margDist.A, numMoments)
}
