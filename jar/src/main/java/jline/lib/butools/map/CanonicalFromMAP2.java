/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class CanonicalFromMAP2 {
    private CanonicalFromMAP2() {}

    /**
     * Returns the canonical form of an order-2 Markovian arrival process.
     *
     * @param D0 The D0 matrix of the MAP(2)
     * @param D1 The D1 matrix of the MAP(2)
     * @param prec Numerical precision to check the input
     * @return Pair of (G0, G1) matrices in canonical form
     */
    public static Pair<Matrix, Matrix> canonicalFromMAP2(Matrix D0, Matrix D1, double prec) {
        if (D0.getNumRows() != 2) {
            throw new IllegalArgumentException("CanonicalFromMAP2: size is not 2!");
        }
        if (!CheckMAPRepresentation.checkMAPRepresentation(D0, D1, prec)) {
            throw new IllegalArgumentException("CanonicalFromMAP2: Input isn't a valid MAP representation!");
        }

        double[] moms = MarginalMomentsFromMAP.marginalMomentsFromMAP(D0, D1, 3);
        double[] corr1 = LagCorrelationsFromMAP.lagCorrelationsFromMAP(D0, D1, 1);

        return MAP2FromMoments.map2FromMoments(moms, corr1[0]);
    }

    public static Pair<Matrix, Matrix> canonicalFromMAP2(Matrix D0, Matrix D1) {
        return canonicalFromMAP2(D0, D1, 1e-14);
    }
}
