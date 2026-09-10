/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.Matrix;

public final class CheckRAPRepresentation {
    private CheckRAPRepresentation() {}

    /**
     * Checks if the given matrices define a valid Rational Arrival Process (RAP) representation.
     *
     * The dominance test follows BuTools: the eigenvalues are split into the
     * real ones and the complex ones, and the representation is rejected only
     * when the largest real part among the COMPLEX eigenvalues strictly exceeds
     * the largest real eigenvalue. A tie is accepted, which is what the MATLAB
     * and Python ports do (they emit a verbose note and return true). Picking
     * the single eigenvalue of smallest modulus of the real part and demanding
     * that it be real, as an earlier version of this port did, rejects valid
     * RAPs whose dominant real eigenvalue is tied in real part with a complex
     * pair, and does so nondeterministically because the tie is broken by the
     * order in which eig() happens to return the eigenvalues.
     *
     * @param H0   The H0 matrix of the RAP
     * @param H1   The H1 matrix of the RAP
     * @param prec Numerical precision. Default 1e-12, as in BuTools.
     * @return True if the representation is valid, false otherwise
     */
    public static boolean checkRAPRepresentation(Matrix H0, Matrix H1, double prec) {
        if (H0.getNumRows() != H0.getNumCols()) {
            return false;
        }
        if (H1.getNumRows() != H1.getNumCols()) {
            return false;
        }
        if (H0.getNumRows() != H1.getNumRows()) {
            return false;
        }

        int n = H0.getNumRows();

        Matrix Q = H0.add(H1);
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += Q.get(i, j);
            }
            if (Math.abs(rowSum) > prec) {
                return false;
            }
        }

        List<Complex> eigenvalues = H0.eig();
        for (Complex ev : eigenvalues) {
            if (ev.getReal() >= -prec) {
                return false;
            }
        }

        double maxReal = Double.NEGATIVE_INFINITY;
        double maxComplex = Double.NEGATIVE_INFINITY;
        for (Complex ev : eigenvalues) {
            if (Math.abs(ev.getImaginary()) <= prec) {
                maxReal = Math.max(maxReal, ev.getReal());
            } else {
                maxComplex = Math.max(maxComplex, ev.getReal());
            }
        }

        // Reject only when the complex eigenvalues dominate by more than the
        // working precision. BuTools compares the two maxima directly, which is
        // exact in MATLAB for this example but not here: eig() returns the real
        // eigenvalue as -1.0 and the complex pair with real part
        // -0.9999999999999999, so a bare comparison rejects the pair on a 1e-16
        // tie. The tolerance is the same one the row-sum and eigenvalue-sign
        // tests above use, so a genuinely dominant complex pair is still caught.
        if (maxReal < maxComplex - prec) {
            return false;
        }

        return true;
    }

    public static boolean checkRAPRepresentation(Matrix H0, Matrix H1) {
        return checkRAPRepresentation(H0, H1, 1e-12);
    }
}
