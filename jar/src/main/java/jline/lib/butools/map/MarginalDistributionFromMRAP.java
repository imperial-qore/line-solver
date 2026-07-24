/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.lib.butools.mc.CRPSolve;
import jline.lib.butools.ph.PHRepresentation;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MarginalDistributionFromMRAP {
    private MarginalDistributionFromMRAP() {}

    /**
     * Returns the phase type distributed marginal distribution
     * of a continuous marked rational arrival process.
     *
     * @param H The H0...HN matrices of the MRAP (as MatrixCell)
     * @param prec Numerical precision for checking if the input is valid
     * @return The PHRepresentation containing alpha (initial vector) and A (matrix parameter)
     */
    public static PHRepresentation marginalDistributionFromMRAP(MatrixCell H, double prec) {
        if (!CheckMRAPRepresentation.checkMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException(
                    "MarginalDistributionFromMRAP: Input isn't a valid MRAP representation!");
        }

        // Sum H1...HN
        Matrix sumH = H.get(1).copy();
        for (int i = 2; i < H.size(); i++) {
            sumH = sumH.add(H.get(i));
        }

        // KEY CHANGE from DMRAP: iH0 = (-H0)^{-1} instead of (I-H0)^{-1}
        Matrix iH0 = H.get(0).neg().inv();
        Matrix P = iH0.mult(sumH);

        Matrix alpha = CRPSolve.drpSolve(P);

        return new PHRepresentation(alpha, H.get(0));
    }

    public static PHRepresentation marginalDistributionFromMRAP(MatrixCell H) {
        return marginalDistributionFromMRAP(H, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static PHRepresentation marginalDistributionFromMRAP(Matrix[] H, double prec) {
        if (!CheckMRAPRepresentation.checkMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException(
                    "MarginalDistributionFromMRAP: Input isn't a valid MRAP representation!");
        }

        Matrix sumH = H[1].copy();
        for (int i = 2; i < H.length; i++) {
            sumH = sumH.add(H[i]);
        }

        Matrix iH0 = H[0].neg().inv();
        Matrix P = iH0.mult(sumH);

        Matrix alpha = CRPSolve.drpSolve(P);

        return new PHRepresentation(alpha, H[0]);
    }

    public static PHRepresentation marginalDistributionFromMRAP(Matrix[] H) {
        return marginalDistributionFromMRAP(H, 1e-14);
    }
}
