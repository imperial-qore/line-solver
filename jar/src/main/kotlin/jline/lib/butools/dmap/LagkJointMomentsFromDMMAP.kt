/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Returns the lag-L joint moments of a discrete marked Markovian arrival process.
 *
 * @param D The D0...DN matrices of the DMMAP (as MatrixCell)
 * @param K The dimension of the matrix of joint moments to compute.
 *          If K=0, the MxM joint moments will be computed.
 * @param L The lag at which the joint moments are computed. Default is 1.
 * @param prec Numerical precision to check if the input is valid.
 * @return List of matrices containing the lag-L joint moments
 */
fun lagkJointMomentsFromDMMAP(
    D: MatrixCell,
    K: Int = 0,
    L: Int = 1,
    prec: Double = 1e-14
): MatrixCell {
    if (!checkDMMAPRepresentation(D, prec)) {
        throw IllegalArgumentException("LagkJointMomentsFromDMMAP: Input isn't a valid DMMAP representation!")
    }

    val actualK = if (K == 0) D[0].numRows - 1 else K

    return lagkJointMomentsFromDMRAP(D, actualK, L, prec)
}

/**
 * Overload for Array<Matrix>.
 */
fun lagkJointMomentsFromDMMAP(
    D: Array<Matrix>,
    K: Int = 0,
    L: Int = 1,
    prec: Double = 1e-14
): MatrixCell {
    val cell = MatrixCell(D.size)
    for (i in D.indices) {
        cell[i] = D[i]
    }
    return lagkJointMomentsFromDMMAP(cell, K, L, prec)
}
