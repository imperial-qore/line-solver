/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.dph.MGFromMoments.MGRepresentation;
import jline.lib.butools.dph.MomentsFromDPH;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MarginalMomentsFromDMMAP {
    private MarginalMomentsFromDMMAP() {}

    /**
     * Returns the moments of the marginal distribution of a discrete
     * marked Markovian arrival process.
     *
     * @param D The D0...DN matrices of the DMMAP (as MatrixCell)
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     * @param prec Numerical precision for checking if the input is valid
     * @return The vector of moments
     */
    public static double[] marginalMomentsFromDMMAP(MatrixCell D, int K, double prec) {
        int numMoments = (K == 0) ? 2 * D.get(0).getNumRows() - 1 : K;

        if (!CheckDMMAPRepresentation.checkDMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("MarginalMomentsFromDMMAP: Input isn't a valid DMMAP representation!");
        }

        MGRepresentation margDist = MarginalDistributionFromDMMAP.marginalDistributionFromDMMAP(D, prec);
        return MomentsFromDPH.momentsFromDPH(margDist.getAlpha(), margDist.getA(), numMoments);
    }

    public static double[] marginalMomentsFromDMMAP(MatrixCell D, int K) {
        return marginalMomentsFromDMMAP(D, K, 1e-14);
    }

    public static double[] marginalMomentsFromDMMAP(MatrixCell D) {
        return marginalMomentsFromDMMAP(D, 0, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static double[] marginalMomentsFromDMMAP(Matrix[] D, int K, double prec) {
        int numMoments = (K == 0) ? 2 * D[0].getNumRows() - 1 : K;

        if (!CheckDMMAPRepresentation.checkDMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("MarginalMomentsFromDMMAP: Input isn't a valid DMMAP representation!");
        }

        MGRepresentation margDist = MarginalDistributionFromDMMAP.marginalDistributionFromDMMAP(D, prec);
        return MomentsFromDPH.momentsFromDPH(margDist.getAlpha(), margDist.getA(), numMoments);
    }

    public static double[] marginalMomentsFromDMMAP(Matrix[] D, int K) {
        return marginalMomentsFromDMMAP(D, K, 1e-14);
    }

    public static double[] marginalMomentsFromDMMAP(Matrix[] D) {
        return marginalMomentsFromDMMAP(D, 0, 1e-14);
    }
}
