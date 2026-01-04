/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.dph.MGRepresentation
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Returns the discrete phase type distributed marginal distribution
 * of a discrete marked Markovian arrival process.
 *
 * @param D The D0...DN matrices of the DMMAP (as MatrixCell)
 * @param prec Numerical precision for checking if the input is valid
 * @return The MGRepresentation containing alpha (initial vector) and A (transient generator)
 */
fun marginalDistributionFromDMMAP(D: MatrixCell, prec: Double = 1e-14): MGRepresentation {
    if (!checkDMMAPRepresentation(D, prec)) {
        throw IllegalArgumentException("MarginalDistributionFromDMMAP: Input isn't a valid DMMAP representation!")
    }

    return marginalDistributionFromDMRAP(D, prec)
}

/**
 * Overload for Array<Matrix>.
 */
fun marginalDistributionFromDMMAP(D: Array<Matrix>, prec: Double = 1e-14): MGRepresentation {
    if (!checkDMMAPRepresentation(D, prec)) {
        throw IllegalArgumentException("MarginalDistributionFromDMMAP: Input isn't a valid DMMAP representation!")
    }

    return marginalDistributionFromDMRAP(D, prec)
}
