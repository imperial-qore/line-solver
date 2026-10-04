/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class DMAP2FromMoments {
    private DMAP2FromMoments() {}

    /**
     * Returns a discrete MAP(2) which has the same 3 marginal
     * moments and lag-1 autocorrelation as given.
     *
     * @param moms First three marginal moments of the inter-arrival times
     * @param corr1 The lag-1 autocorrelation of the inter-arrival times
     * @return Pair of (D0, D1) matrices of the discrete MAP(2)
     *
     * <p>Note: Raises an exception if the moments are not feasible with a DMAP(2).
     */
    public static Pair<Matrix, Matrix> dmap2FromMoments(double[] moms, double corr1) {
        double m1 = moms[0];
        double m2 = moms[1];

        // Nm = [1, moms(1); moms(1), corr1*(moms(2)-moms(1)^2)+moms(1)^2]
        Matrix Nm = new Matrix(2, 2);
        Nm.set(0, 0, 1.0);
        Nm.set(0, 1, m1);
        Nm.set(1, 0, m1);
        Nm.set(1, 1, corr1 * (m2 - m1 * m1) + m1 * m1);

        Pair<Matrix, Matrix> hh = DRAPFromMoments.drapFromMoments(moms, Nm);
        Matrix H0 = hh.getLeft();
        Matrix H1 = hh.getRight();

        return CanonicalFromDMAP2.canonicalFromDMAP2(H0, H1);
    }
}
