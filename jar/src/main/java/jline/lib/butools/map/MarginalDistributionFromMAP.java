/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.lib.butools.mc.CRPSolve;
import jline.lib.butools.ph.PHRepresentation;
import jline.util.matrix.Matrix;

public final class MarginalDistributionFromMAP {
    private MarginalDistributionFromMAP() {}

    /**
     * Returns the phase type distributed marginal distribution
     * of a Markovian arrival process.
     *
     * @param D0 The D0 matrix of the Markovian arrival process
     * @param D1 The D1 matrix of the Markovian arrival process
     * @return PHRepresentation containing alpha (initial vector) and A (transient generator)
     */
    public static PHRepresentation marginalDistributionFromMAP(Matrix D0, Matrix D1) {
        return marginalDistributionFromRAP(D0, D1);
    }

    /**
     * Returns the matrix exponential distributed marginal distribution
     * of a rational arrival process.
     *
     * @param H0 The H0 matrix of the rational arrival process
     * @param H1 The H1 matrix of the rational arrival process
     * @return PHRepresentation containing alpha (initial vector) and A (matrix parameter)
     */
    public static PHRepresentation marginalDistributionFromRAP(Matrix H0, Matrix H1) {
        // P = inv(-H0) * H1
        Matrix P = H0.neg().inv().mult(H1);

        // alpha = DRPSolve(P)
        Matrix alpha = CRPSolve.drpSolve(P);

        return new PHRepresentation(alpha, H0);
    }
}
