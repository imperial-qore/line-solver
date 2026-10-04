/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.lib.butools.dph.MomentsFromDPH;
import jline.util.matrix.Matrix;

public final class MarginalMomentsFromDMAP {
    private MarginalMomentsFromDMAP() {}

    /**
     * Returns the moments of the marginal distribution of a discrete
     * Markovian arrival process.
     *
     * @param D0 The D0 matrix of the discrete Markovian arrival process
     * @param D1 The D1 matrix of the discrete Markovian arrival process
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     * @param prec Numerical precision for checking if the input is valid
     * @return The vector of moments
     */
    public static double[] marginalMomentsFromDMAP(Matrix D0, Matrix D1, int K, double prec) {
        int numMoments = (K == 0) ? 2 * D0.getNumRows() - 1 : K;

        if (!CheckDMAPRepresentation.checkDMAPRepresentation(D0, D1, prec)) {
            throw new IllegalArgumentException("MarginalMomentsFromDMAP: Input isn't a valid DMAP representation!");
        }

        MGRepresentation margDist = MarginalDistributionFromDMAP.marginalDistributionFromDMAP(D0, D1, prec);
        return MomentsFromDPH.momentsFromDPH(margDist.getAlpha(), margDist.getA(), numMoments);
    }

    public static double[] marginalMomentsFromDMAP(Matrix D0, Matrix D1, int K) {
        return marginalMomentsFromDMAP(D0, D1, K, 1e-14);
    }

    public static double[] marginalMomentsFromDMAP(Matrix D0, Matrix D1) {
        return marginalMomentsFromDMAP(D0, D1, 0, 1e-14);
    }
}
