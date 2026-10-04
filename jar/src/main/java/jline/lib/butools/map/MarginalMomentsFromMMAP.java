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

public final class MarginalMomentsFromMMAP {
    private MarginalMomentsFromMMAP() {}

    /**
     * Returns the moments of the marginal distribution of a continuous
     * marked Markovian arrival process.
     *
     * @param D The D0...DN matrices of the MMAP (as MatrixCell)
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     * @param prec Numerical precision for checking if the input is valid
     * @return The vector of moments
     */
    public static double[] marginalMomentsFromMMAP(MatrixCell D, int K, double prec) {
        int numMoments = (K == 0) ? 2 * D.get(0).getNumRows() - 1 : K;

        if (!CheckMMAPRepresentation.checkMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("MarginalMomentsFromMMAP: Input isn't a valid MMAP representation!");
        }

        PHRepresentation margDist = MarginalDistributionFromMMAP.marginalDistributionFromMMAP(D, prec);
        return MomentsFromME.momentsFromME(margDist.getAlpha(), margDist.getA(), numMoments);
    }

    public static double[] marginalMomentsFromMMAP(MatrixCell D, int K) {
        return marginalMomentsFromMMAP(D, K, 1e-14);
    }

    public static double[] marginalMomentsFromMMAP(MatrixCell D) {
        return marginalMomentsFromMMAP(D, 0, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static double[] marginalMomentsFromMMAP(Matrix[] D, int K, double prec) {
        int numMoments = (K == 0) ? 2 * D[0].getNumRows() - 1 : K;

        if (!CheckMMAPRepresentation.checkMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("MarginalMomentsFromMMAP: Input isn't a valid MMAP representation!");
        }

        PHRepresentation margDist = MarginalDistributionFromMMAP.marginalDistributionFromMMAP(D, prec);
        return MomentsFromME.momentsFromME(margDist.getAlpha(), margDist.getA(), numMoments);
    }

    public static double[] marginalMomentsFromMMAP(Matrix[] D, int K) {
        return marginalMomentsFromMMAP(D, K, 1e-14);
    }

    public static double[] marginalMomentsFromMMAP(Matrix[] D) {
        return marginalMomentsFromMMAP(D, 0, 1e-14);
    }
}
