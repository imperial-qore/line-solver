/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.lib.butools.ph.MomentsFromME;
import jline.lib.butools.ph.PHRepresentation;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MarginalMomentsFromMRAP {
    private MarginalMomentsFromMRAP() {}

    /**
     * Returns the moments of the marginal distribution of a continuous
     * marked rational arrival process.
     *
     * @param H The H0...HN matrices of the MRAP (as MatrixCell)
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     * @param prec Numerical precision for checking if the input is valid
     * @return The vector of moments
     */
    public static double[] marginalMomentsFromMRAP(MatrixCell H, int K, double prec) {
        int numMoments = (K == 0) ? 2 * H.get(0).getNumRows() - 1 : K;

        if (!CheckMRAPRepresentation.checkMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException("MarginalMomentsFromMRAP: Input isn't a valid MRAP representation!");
        }

        PHRepresentation margDist = MarginalDistributionFromMRAP.marginalDistributionFromMRAP(H, prec);
        return MomentsFromME.momentsFromME(margDist.getAlpha(), margDist.getA(), numMoments);
    }

    public static double[] marginalMomentsFromMRAP(MatrixCell H, int K) {
        return marginalMomentsFromMRAP(H, K, 1e-14);
    }

    public static double[] marginalMomentsFromMRAP(MatrixCell H) {
        return marginalMomentsFromMRAP(H, 0, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static double[] marginalMomentsFromMRAP(Matrix[] H, int K, double prec) {
        int numMoments = (K == 0) ? 2 * H[0].getNumRows() - 1 : K;

        if (!CheckMRAPRepresentation.checkMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException("MarginalMomentsFromMRAP: Input isn't a valid MRAP representation!");
        }

        PHRepresentation margDist = MarginalDistributionFromMRAP.marginalDistributionFromMRAP(H, prec);
        return MomentsFromME.momentsFromME(margDist.getAlpha(), margDist.getA(), numMoments);
    }

    public static double[] marginalMomentsFromMRAP(Matrix[] H, int K) {
        return marginalMomentsFromMRAP(H, K, 1e-14);
    }

    public static double[] marginalMomentsFromMRAP(Matrix[] H) {
        return marginalMomentsFromMRAP(H, 0, 1e-14);
    }
}
