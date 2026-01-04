/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.util.matrix.Matrix

/**
 * Returns a discrete MAP(2) which has the same 3 marginal
 * moments and lag-1 autocorrelation as given.
 *
 * @param moms First three marginal moments of the inter-arrival times
 * @param corr1 The lag-1 autocorrelation of the inter-arrival times
 * @return Pair of (D0, D1) matrices of the discrete MAP(2)
 *
 * Note: Raises an exception if the moments are not feasible with a DMAP(2).
 */
fun dmap2FromMoments(moms: DoubleArray, corr1: Double): Pair<Matrix, Matrix> {
    val m1 = moms[0]
    val m2 = moms[1]

    // Nm = [1, moms(1); moms(1), corr1*(moms(2)-moms(1)^2)+moms(1)^2]
    val Nm = Matrix(2, 2)
    Nm[0, 0] = 1.0
    Nm[0, 1] = m1
    Nm[1, 0] = m1
    Nm[1, 1] = corr1 * (m2 - m1 * m1) + m1 * m1

    val (H0, H1) = drapFromMoments(moms, Nm)

    return canonicalFromDMAP2(H0, H1)
}
