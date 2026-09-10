/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.lib.butools.ph.MomentsFromME;
import jline.lib.butools.ph.PHRepresentation;
import jline.util.matrix.Matrix;

public final class MarginalMomentsFromMAP {
    private MarginalMomentsFromMAP() {}

    /**
     * Returns the moments of the marginal distribution of a Markovian arrival process.
     *
     * @param D0 The D0 matrix of the Markovian arrival process
     * @param D1 The D1 matrix of the Markovian arrival process
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     * @return The vector of moments
     */
    public static double[] marginalMomentsFromMAP(Matrix D0, Matrix D1, int K) {
        int numMoments = (K == 0) ? 2 * D0.getNumRows() - 1 : K;
        PHRepresentation ph = MarginalDistributionFromMAP.marginalDistributionFromMAP(D0, D1);
        return MomentsFromME.momentsFromME(ph.getAlpha(), ph.getA(), numMoments);
    }

    public static double[] marginalMomentsFromMAP(Matrix D0, Matrix D1) {
        return marginalMomentsFromMAP(D0, D1, 0);
    }

    /**
     * Returns the moments of the marginal distribution of a rational arrival process.
     *
     * @param H0 The H0 matrix of the rational arrival process
     * @param H1 The H1 matrix of the rational arrival process
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     * @return The vector of moments
     */
    public static double[] marginalMomentsFromRAP(Matrix H0, Matrix H1, int K) {
        int numMoments = (K == 0) ? 2 * H0.getNumRows() - 1 : K;
        PHRepresentation ph = MarginalDistributionFromMAP.marginalDistributionFromRAP(H0, H1);
        return MomentsFromME.momentsFromME(ph.getAlpha(), ph.getA(), numMoments);
    }

    public static double[] marginalMomentsFromRAP(Matrix H0, Matrix H1) {
        return marginalMomentsFromRAP(H0, H1, 0);
    }
}
