/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.dph.MGRepresentation
import jline.lib.butools.mc.drpSolve
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Returns the matrix geometrically distributed marginal distribution
 * of a discrete marked rational arrival process.
 *
 * @param H The H0...HN matrices of the DMRAP (as MatrixCell)
 * @param prec Numerical precision for checking if the input is valid
 * @return The MGRepresentation containing alpha (initial vector) and A (matrix parameter)
 */
fun marginalDistributionFromDMRAP(H: MatrixCell, prec: Double = 1e-14): MGRepresentation {
    if (!checkDMRAPRepresentation(H, prec)) {
        throw IllegalArgumentException("MarginalDistributionFromDMRAP: Input isn't a valid DMRAP representation!")
    }

    val n = H[0].numRows
    val I = Matrix.eye(n)

    // Sum H1...HN
    var sumH = H[1].copy()
    for (i in 2 until H.size()) {
        sumH = sumH.add(H[i])
    }

    // alpha = DRPSolve(inv(I - H0) * sumH)
    val IminusH0 = I.sub(H[0])
    val IminusH0inv = IminusH0.inv()
    val P = IminusH0inv.mult(sumH)
    val alpha = drpSolve(P)

    return MGRepresentation(alpha, H[0])
}

/**
 * Overload for Array<Matrix>.
 */
fun marginalDistributionFromDMRAP(H: Array<Matrix>, prec: Double = 1e-14): MGRepresentation {
    if (!checkDMRAPRepresentation(H, prec)) {
        throw IllegalArgumentException("MarginalDistributionFromDMRAP: Input isn't a valid DMRAP representation!")
    }

    val n = H[0].numRows
    val I = Matrix.eye(n)

    // Sum H1...HN
    var sumH = H[1].copy()
    for (i in 2 until H.size) {
        sumH = sumH.add(H[i])
    }

    // alpha = DRPSolve(inv(I - H0) * sumH)
    val IminusH0 = I.sub(H[0])
    val IminusH0inv = IminusH0.inv()
    val P = IminusH0inv.mult(sumH)
    val alpha = drpSolve(P)

    return MGRepresentation(alpha, H[0])
}
