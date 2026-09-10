/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.util.matrix.Matrix;

public final class MarginalDistributionFromDMAP {
    private MarginalDistributionFromDMAP() {}

    /**
     * Returns the discrete phase type distributed marginal distribution
     * of a discrete Markovian arrival process.
     *
     * @param D0 The D0 matrix of the discrete Markovian arrival process
     * @param D1 The D1 matrix of the discrete Markovian arrival process
     * @param prec Numerical precision for checking if the input is valid
     * @return The MGRepresentation containing alpha (initial vector) and A (transient generator)
     */
    public static MGRepresentation marginalDistributionFromDMAP(Matrix D0, Matrix D1, double prec) {
        if (!CheckDMAPRepresentation.checkDMAPRepresentation(D0, D1, prec)) {
            throw new IllegalArgumentException("MarginalDistributionFromDMAP: Input isn't a valid DMAP representation!");
        }

        return MarginalDistributionFromDRAP.marginalDistributionFromDRAP(D0, D1, prec);
    }

    public static MGRepresentation marginalDistributionFromDMAP(Matrix D0, Matrix D1) {
        return marginalDistributionFromDMAP(D0, D1, 1e-14);
    }
}
