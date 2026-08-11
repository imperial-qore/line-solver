/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.lib.butools.dph.MomentsFromMG;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MarginalMomentsFromDMRAP {
    private MarginalMomentsFromDMRAP() {}

    /**
     * Returns the moments of the marginal distribution of a discrete
     * marked rational arrival process.
     *
     * @param H The H0...HN matrices of the DMRAP (as MatrixCell)
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     * @param prec Numerical precision for checking if the input is valid
     * @return The vector of moments
     */
    public static double[] marginalMomentsFromDMRAP(MatrixCell H, int K, double prec) {
        int numMoments = (K == 0) ? 2 * H.get(0).getNumRows() - 1 : K;

        if (!CheckDMRAPRepresentation.checkDMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException("MarginalMomentsFromDMRAP: Input isn't a valid DMRAP representation!");
        }

        MGRepresentation margDist = MarginalDistributionFromDMRAP.marginalDistributionFromDMRAP(H, prec);
        return MomentsFromMG.momentsFromMG(margDist.getAlpha(), margDist.getA(), numMoments);
    }

    public static double[] marginalMomentsFromDMRAP(MatrixCell H, int K) {
        return marginalMomentsFromDMRAP(H, K, 1e-14);
    }

    public static double[] marginalMomentsFromDMRAP(MatrixCell H) {
        return marginalMomentsFromDMRAP(H, 0, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static double[] marginalMomentsFromDMRAP(Matrix[] H, int K, double prec) {
        int numMoments = (K == 0) ? 2 * H[0].getNumRows() - 1 : K;

        if (!CheckDMRAPRepresentation.checkDMRAPRepresentation(H, prec)) {
            throw new IllegalArgumentException("MarginalMomentsFromDMRAP: Input isn't a valid DMRAP representation!");
        }

        MGRepresentation margDist = MarginalDistributionFromDMRAP.marginalDistributionFromDMRAP(H, prec);
        return MomentsFromMG.momentsFromMG(margDist.getAlpha(), margDist.getA(), numMoments);
    }

    public static double[] marginalMomentsFromDMRAP(Matrix[] H, int K) {
        return marginalMomentsFromDMRAP(H, K, 1e-14);
    }

    public static double[] marginalMomentsFromDMRAP(Matrix[] H) {
        return marginalMomentsFromDMRAP(H, 0, 1e-14);
    }
}
