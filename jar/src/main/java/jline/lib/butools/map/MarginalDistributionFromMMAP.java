/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.lib.butools.ph.PHRepresentation;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MarginalDistributionFromMMAP {
    private MarginalDistributionFromMMAP() {}

    /**
     * Returns the phase type distributed marginal distribution
     * of a continuous marked Markovian arrival process.
     *
     * @param D The D0...DN matrices of the MMAP (as MatrixCell)
     * @param prec Numerical precision for checking if the input is valid
     * @return The PHRepresentation containing alpha (initial vector) and A (transient generator)
     */
    public static PHRepresentation marginalDistributionFromMMAP(MatrixCell D, double prec) {
        if (!CheckMMAPRepresentation.checkMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("MarginalDistributionFromMMAP: Input isn't a valid MMAP representation!");
        }

        return MarginalDistributionFromMRAP.marginalDistributionFromMRAP(D, prec);
    }

    public static PHRepresentation marginalDistributionFromMMAP(MatrixCell D) {
        return marginalDistributionFromMMAP(D, 1e-14);
    }

    public static PHRepresentation marginalDistributionFromMMAP(Matrix[] D, double prec) {
        if (!CheckMMAPRepresentation.checkMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("MarginalDistributionFromMMAP: Input isn't a valid MMAP representation!");
        }

        return MarginalDistributionFromMRAP.marginalDistributionFromMRAP(D, prec);
    }

    public static PHRepresentation marginalDistributionFromMMAP(Matrix[] D) {
        return marginalDistributionFromMMAP(D, 1e-14);
    }
}
