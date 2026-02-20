/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map

import jline.util.matrix.Matrix

/**
 * Returns the canonical form of an order-2 Markovian arrival process.
 *
 * @param D0 The D0 matrix of the MAP(2)
 * @param D1 The D1 matrix of the MAP(2)
 * @param prec Numerical precision to check the input
 * @return Pair of (G0, G1) matrices in canonical form
 */
@JvmOverloads
fun canonicalFromMAP2(D0: Matrix, D1: Matrix, prec: Double = 1e-14): Pair<Matrix, Matrix> {
    if (D0.numRows != 2) {
        throw IllegalArgumentException("CanonicalFromMAP2: size is not 2!")
    }
    if (!checkMAPRepresentation(D0, D1, prec)) {
        throw IllegalArgumentException("CanonicalFromMAP2: Input isn't a valid MAP representation!")
    }

    val moms = marginalMomentsFromMAP(D0, D1, 3)
    val corr1 = lagCorrelationsFromMAP(D0, D1, 1)

    return map2FromMoments(moms, corr1[0])
}
