/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.dph.MGRepresentation
import jline.lib.butools.mc.drpSolve
import jline.util.matrix.Matrix

/**
 * Returns the matrix geometrically distributed marginal distribution
 * of a discrete rational arrival process.
 *
 * @param H0 The H0 matrix of the discrete rational arrival process
 * @param H1 The H1 matrix of the discrete rational arrival process
 * @param prec Numerical precision for checking if the input is valid
 * @return The MGRepresentation containing alpha (initial vector) and A (matrix parameter)
 */
fun marginalDistributionFromDRAP(H0: Matrix, H1: Matrix, prec: Double = 1e-14): MGRepresentation {
    if (!checkDRAPRepresentation(H0, H1, prec)) {
        throw IllegalArgumentException("MarginalDistributionFromDRAP: Input isn't a valid DRAP representation!")
    }

    val n = H0.numRows
    val I = Matrix.eye(n)

    // alpha = DRPSolve(inv(I - H0) * H1)
    val IminusH0 = I.sub(H0)
    val IminusH0inv = IminusH0.inv()
    val P = IminusH0inv.mult(H1)
    val alpha = drpSolve(P)

    return MGRepresentation(alpha, H0)
}
