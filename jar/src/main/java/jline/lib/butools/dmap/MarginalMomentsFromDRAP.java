/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.lib.butools.dph.MomentsFromMG;
import jline.util.matrix.Matrix;

public final class MarginalMomentsFromDRAP {
    private MarginalMomentsFromDRAP() {}

    /**
     * Returns the moments of the marginal distribution of a discrete
     * rational arrival process.
     *
     * @param H0 The H0 matrix of the discrete rational arrival process
     * @param H1 The H1 matrix of the discrete rational arrival process
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     * @param prec Numerical precision for checking if the input is valid
     * @return The vector of moments
     */
    public static double[] marginalMomentsFromDRAP(Matrix H0, Matrix H1, int K, double prec) {
        int numMoments = (K == 0) ? 2 * H0.getNumRows() - 1 : K;

        if (!CheckDRAPRepresentation.checkDRAPRepresentation(H0, H1, prec)) {
            throw new IllegalArgumentException("MarginalMomentsFromDRAP: Input isn't a valid DRAP representation!");
        }

        MGRepresentation margDist = MarginalDistributionFromDRAP.marginalDistributionFromDRAP(H0, H1, prec);
        return MomentsFromMG.momentsFromMG(margDist.getAlpha(), margDist.getA(), numMoments);
    }

    public static double[] marginalMomentsFromDRAP(Matrix H0, Matrix H1, int K) {
        return marginalMomentsFromDRAP(H0, H1, K, 1e-14);
    }

    public static double[] marginalMomentsFromDRAP(Matrix H0, Matrix H1) {
        return marginalMomentsFromDRAP(H0, H1, 0, 1e-14);
    }
}
