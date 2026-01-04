/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Creates a discrete rational arrival process that has the
 * same marginal and lag-1 joint moments as given.
 *
 * @param moms The list of marginal moments. To obtain a rational
 *             process of order M, 2*M-1 marginal moments are required.
 * @param Nm The matrix of lag-1 joint moments
 * @return Pair of (H0, H1) matrices of the discrete rational process
 */
fun drapFromMoments(moms: DoubleArray, Nm: Matrix): Pair<Matrix, Matrix> {
    val NmCell = MatrixCell(1)
    NmCell[0] = Nm

    val H = dmrapFromMoments(moms, NmCell)

    return Pair(H[0], H[1])
}
