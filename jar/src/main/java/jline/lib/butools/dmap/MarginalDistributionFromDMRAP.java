/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.lib.butools.mc.CRPSolve;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MarginalDistributionFromDMRAP {
    private MarginalDistributionFromDMRAP() {}

    /**
     * Returns the matrix geometrically distributed marginal distribution
     * of a discrete marked rational arrival process.
     *
     * @param H The H0...HN matrices of the DMRAP (as MatrixCell)
     * @param prec Numerical precision for checking if the input is valid
     * @return The MGRepresentation containing alpha (initial vector) and A (matrix parameter)
     */
    public static MGRepresentation marginalDistributionFromDMRAP(MatrixCell H, double prec) {
        if (!CheckDMRAPRepresentation.checkDMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException(
                    "MarginalDistributionFromDMRAP: Input isn't a valid DMRAP representation!");
        }

        int n = H.get(0).getNumRows();
        Matrix I = Matrix.eye(n);

        // Sum H1...HN
        Matrix sumH = H.get(1).copy();
        for (int i = 2; i < H.size(); i++) {
            sumH = sumH.add(H.get(i));
        }

        // alpha = DRPSolve(inv(I - H0) * sumH)
        Matrix IminusH0 = I.sub(H.get(0));
        Matrix IminusH0inv = IminusH0.inv();
        Matrix P = IminusH0inv.mult(sumH);
        Matrix alpha = CRPSolve.drpSolve(P);

        return new MGRepresentation(alpha, H.get(0));
    }

    public static MGRepresentation marginalDistributionFromDMRAP(MatrixCell H) {
        return marginalDistributionFromDMRAP(H, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static MGRepresentation marginalDistributionFromDMRAP(Matrix[] H, double prec) {
        if (!CheckDMRAPRepresentation.checkDMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException(
                    "MarginalDistributionFromDMRAP: Input isn't a valid DMRAP representation!");
        }

        int n = H[0].getNumRows();
        Matrix I = Matrix.eye(n);

        Matrix sumH = H[1].copy();
        for (int i = 2; i < H.length; i++) {
            sumH = sumH.add(H[i]);
        }

        Matrix IminusH0 = I.sub(H[0]);
        Matrix IminusH0inv = IminusH0.inv();
        Matrix P = IminusH0inv.mult(sumH);
        Matrix alpha = CRPSolve.drpSolve(P);

        return new MGRepresentation(alpha, H[0]);
    }

    public static MGRepresentation marginalDistributionFromDMRAP(Matrix[] H) {
        return marginalDistributionFromDMRAP(H, 1e-14);
    }
}
