/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap

import jline.lib.butools.dph.MGRepresentation
import jline.util.matrix.Matrix

/**
 * Returns the discrete phase type distributed marginal distribution
 * of a discrete Markovian arrival process.
 *
 * @param D0 The D0 matrix of the discrete Markovian arrival process
 * @param D1 The D1 matrix of the discrete Markovian arrival process
 * @param prec Numerical precision for checking if the input is valid
 * @return The MGRepresentation containing alpha (initial vector) and A (transient generator)
 */
fun marginalDistributionFromDMAP(D0: Matrix, D1: Matrix, prec: Double = 1e-14): MGRepresentation {
    if (!checkDMAPRepresentation(D0, D1, prec)) {
        throw IllegalArgumentException("MarginalDistributionFromDMAP: Input isn't a valid DMAP representation!")
    }

    return marginalDistributionFromDRAP(D0, D1, prec)
}
