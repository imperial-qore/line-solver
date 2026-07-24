/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MarginalDistributionFromDMMAP {
    private MarginalDistributionFromDMMAP() {}

    /**
     * Returns the discrete phase type distributed marginal distribution
     * of a discrete marked Markovian arrival process.
     *
     * @param D The D0...DN matrices of the DMMAP (as MatrixCell)
     * @param prec Numerical precision for checking if the input is valid
     * @return The MGRepresentation containing alpha (initial vector) and A (transient generator)
     */
    public static MGRepresentation marginalDistributionFromDMMAP(MatrixCell D, double prec) {
        if (!CheckDMMAPRepresentation.checkDMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("MarginalDistributionFromDMMAP: Input isn't a valid DMMAP representation!");
        }

        return MarginalDistributionFromDMRAP.marginalDistributionFromDMRAP(D, prec);
    }

    public static MGRepresentation marginalDistributionFromDMMAP(MatrixCell D) {
        return marginalDistributionFromDMMAP(D, 1e-14);
    }

    public static MGRepresentation marginalDistributionFromDMMAP(Matrix[] D, double prec) {
        if (!CheckDMMAPRepresentation.checkDMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("MarginalDistributionFromDMMAP: Input isn't a valid DMMAP representation!");
        }

        return MarginalDistributionFromDMRAP.marginalDistributionFromDMRAP(D, prec);
    }

    public static MGRepresentation marginalDistributionFromDMMAP(Matrix[] D) {
        return marginalDistributionFromDMMAP(D, 1e-14);
    }
}
