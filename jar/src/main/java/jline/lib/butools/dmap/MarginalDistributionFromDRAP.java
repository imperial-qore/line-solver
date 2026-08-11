/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.lib.butools.mc.CRPSolve;
import jline.util.matrix.Matrix;

public final class MarginalDistributionFromDRAP {
    private MarginalDistributionFromDRAP() {}

    /**
     * Returns the matrix geometrically distributed marginal distribution
     * of a discrete rational arrival process.
     *
     * @param H0 The H0 matrix of the discrete rational arrival process
     * @param H1 The H1 matrix of the discrete rational arrival process
     * @param prec Numerical precision for checking if the input is valid
     * @return The MGRepresentation containing alpha (initial vector) and A (matrix parameter)
     */
    public static MGRepresentation marginalDistributionFromDRAP(Matrix H0, Matrix H1, double prec) {
        if (!CheckDRAPRepresentation.checkDRAPRepresentation(H0, H1, prec)) {
            throw new IllegalArgumentException("MarginalDistributionFromDRAP: Input isn't a valid DRAP representation!");
        }

        int n = H0.getNumRows();
        Matrix I = Matrix.eye(n);

        // alpha = DRPSolve(inv(I - H0) * H1)
        Matrix IminusH0 = I.sub(H0);
        Matrix IminusH0inv = IminusH0.inv();
        Matrix P = IminusH0inv.mult(H1);
        Matrix alpha = CRPSolve.drpSolve(P);

        return new MGRepresentation(alpha, H0);
    }

    public static MGRepresentation marginalDistributionFromDRAP(Matrix H0, Matrix H1) {
        return marginalDistributionFromDRAP(H0, H1, 1e-14);
    }
}
