/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.dph.momentsFromDPH
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Returns the moments of the marginal distribution of a discrete
 * marked Markovian arrival process.
 *
 * @param D The D0...DN matrices of the DMMAP (as MatrixCell)
 * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
 * @param prec Numerical precision for checking if the input is valid
 * @return The vector of moments
 */
fun marginalMomentsFromDMMAP(D: MatrixCell, K: Int = 0, prec: Double = 1e-14): DoubleArray {
    val numMoments = if (K == 0) 2 * D[0].numRows - 1 else K

    if (!checkDMMAPRepresentation(D, prec)) {
        throw IllegalArgumentException("MarginalMomentsFromDMMAP: Input isn't a valid DMMAP representation!")
    }

    val margDist = marginalDistributionFromDMMAP(D, prec)
    return momentsFromDPH(margDist.alpha, margDist.A, numMoments)
}

/**
 * Overload for Array<Matrix>.
 */
fun marginalMomentsFromDMMAP(D: Array<Matrix>, K: Int = 0, prec: Double = 1e-14): DoubleArray {
    val numMoments = if (K == 0) 2 * D[0].numRows - 1 else K

    if (!checkDMMAPRepresentation(D, prec)) {
        throw IllegalArgumentException("MarginalMomentsFromDMMAP: Input isn't a valid DMMAP representation!")
    }

    val margDist = marginalDistributionFromDMMAP(D, prec)
    return momentsFromDPH(margDist.alpha, margDist.A, numMoments)
}
