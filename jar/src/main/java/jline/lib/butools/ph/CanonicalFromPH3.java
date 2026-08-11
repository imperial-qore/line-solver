/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.util.matrix.Matrix;

public final class CanonicalFromPH3 {
    private CanonicalFromPH3() {}

    /**
     * Returns the canonical form of an order-3 phase-type distribution.
     *
     * @param alpha Initial vector of the phase-type distribution.
     * @param A Transient generator of the phase-type distribution.
     * @param prec Numerical precision, default value is 1e-10.
     * @return The PH3Representation containing beta (canonical initial vector) and B (canonical generator).
     *
     * Note: This procedure calculates 5 moments of the input and calls 'ph3From5Moments'.
     */
    public static PH3Representation canonicalFromPH3(Matrix alpha, Matrix A, double prec) {
        if (A.getNumRows() != 3) {
            throw new IllegalArgumentException("CanonicalFromPH3: Dimension is not 3!");
        }

        if (!CheckMERepresentation.checkMERepresentation(alpha, A)) {
            throw new IllegalArgumentException("CanonicalFromPH3: Input isn't a valid ME distribution!");
        }

        double[] moms = MomentsFromME.momentsFromME(alpha, A, 5);
        return PH3From5Moments.ph3From5Moments(moms, prec);
    }

    public static PH3Representation canonicalFromPH3(Matrix alpha, Matrix A) {
        return canonicalFromPH3(alpha, A, 1e-10);
    }

    public static PH3Representation canonicalFromPH3(double[] alpha, Matrix A, double prec) {
        return canonicalFromPH3(new Matrix(alpha), A, prec);
    }

    public static PH3Representation canonicalFromPH3(double[] alpha, Matrix A) {
        return canonicalFromPH3(new Matrix(alpha), A, 1e-10);
    }
}
