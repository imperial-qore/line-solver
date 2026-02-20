/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map

import jline.lib.butools.ph.momentsFromPH
import jline.util.matrix.Matrix

/**
 * Returns the moments of the marginal distribution of a
 * Markovian arrival process.
 *
 * @param D0 The D0 matrix of the Markovian arrival process
 * @param D1 The D1 matrix of the Markovian arrival process
 * @param K Number of moments to compute. If K=0, 2*M-1 moments
 *        are computed. The default value is 0.
 * @return The vector of moments
 */
@JvmOverloads
fun marginalMomentsFromMAP(D0: Matrix, D1: Matrix, K: Int = 0): DoubleArray {
    val numMoments = if (K == 0) 2 * D0.numRows - 1 else K
    val (alpha, A) = marginalDistributionFromMAP(D0, D1)
    return momentsFromPH(alpha, A, numMoments)
}

/**
 * Returns the moments of the marginal distribution of a
 * rational arrival process.
 *
 * @param H0 The H0 matrix of the rational arrival process
 * @param H1 The H1 matrix of the rational arrival process
 * @param K Number of moments to compute. If K=0, 2*M-1 moments
 *        are computed. The default value is 0.
 * @return The vector of moments
 */
@JvmOverloads
fun marginalMomentsFromRAP(H0: Matrix, H1: Matrix, K: Int = 0): DoubleArray {
    val numMoments = if (K == 0) 2 * H0.numRows - 1 else K
    val (alpha, A) = marginalDistributionFromRAP(H0, H1)
    return momentsFromPH(alpha, A, numMoments)
}
