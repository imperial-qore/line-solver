/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.util.matrix.Matrix;
import jline.lib.butools.ph.PH2From3Moments.PH2Representation;

public final class CanonicalFromPH2 {
    private CanonicalFromPH2() {}

    /**
     * Returns the canonical form of an order-2 phase-type distribution.
     *
     * @param alpha Initial vector of the phase-type distribution.
     * @param A Transient generator of the phase-type distribution.
     * @param prec Numerical precision, default value is 1e-14.
     * @return The PH2Representation containing beta (canonical initial vector) and B (canonical generator).
     *
     * Note: This procedure calculates 3 moments of the input and calls 'ph2From3Moments'.
     */
    public static PH2Representation canonicalFromPH2(Matrix alpha, Matrix A, double prec) {
        if (A.getNumRows() != 2) {
            throw new IllegalArgumentException("CanonicalFromPH2: Dimension is not 2!");
        }

        if (!CheckMERepresentation.checkMERepresentation(alpha, A)) {
            throw new IllegalArgumentException("CanonicalFromPH2: Input isn't a valid ME distribution!");
        }

        double[] moms = MomentsFromME.momentsFromME(alpha, A, 3);
        return PH2From3Moments.ph2From3Moments(moms, prec);
    }

    public static PH2Representation canonicalFromPH2(Matrix alpha, Matrix A) {
        return canonicalFromPH2(alpha, A, 1e-14);
    }

    public static PH2Representation canonicalFromPH2(double[] alpha, Matrix A, double prec) {
        return canonicalFromPH2(new Matrix(alpha), A, prec);
    }

    public static PH2Representation canonicalFromPH2(double[] alpha, Matrix A) {
        return canonicalFromPH2(new Matrix(alpha), A, 1e-14);
    }
}
