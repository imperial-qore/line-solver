/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.Matrix;

public final class CheckMERepresentation {
    private CheckMERepresentation() {}

    /**
     * Checks if the given vector and matrix define a valid matrix-
     * exponential representation.
     *
     * @param alpha Initial vector of the matrix-exponential distribution to check
     * @param A Matrix parameter of the matrix-exponential distribution to check
     * @param prec Numerical precision. The default value is 1e-14.
     * @return True, if the matrix is a square matrix, the vector and
     *         the matrix have the same size, the dominant eigenvalue
     *         is negative and real
     *
     * Note: This procedure does not check the positivity of the density!
     * Call 'checkMEPositiveDensity' if it is needed, but keep in mind
     * that it can be time-consuming, while this procedure is fast.
     */
    public static boolean checkMERepresentation(Matrix alpha, Matrix A, double prec) {
        // Check if matrix is square
        if (A.getNumRows() != A.getNumCols()) {
            return false;
        }

        // Check if vector and matrix have same size
        if (alpha.length() != A.getNumRows()) {
            return false;
        }

        // Check if sum of alpha is between 0 and 1
        double sumAlpha = alpha.elementSum();
        if (sumAlpha < -prec * alpha.length() || sumAlpha > 1 + prec * alpha.length()) {
            return false;
        }

        // Check if all eigenvalues have negative real parts
        List<Complex> eigenvalues = A.eig();
        for (Complex ev : eigenvalues) {
            if (ev.getReal() >= prec) {
                return false;
            }
        }

        // Find dominant eigenvalue (closest to zero)
        Complex maxEv = eigenvalues.get(0);
        for (Complex ev : eigenvalues) {
            if (Math.abs(ev.getReal()) < Math.abs(maxEv.getReal())) {
                maxEv = ev;
            }
        }

        // The dominant real part need not be attained by a single eigenvalue, and the
        // scan above then keeps whichever tied eigenvalue came first. A concentrated
        // matrix exponential of order 2n+1 puts its whole spectrum on the line
        // Re = -mu1, so the arbitrary pick returns a complex eigenvalue and the test
        // below rejects a valid ME even though the real eigenvalue -mu1 is equally
        // dominant. Among the eigenvalues attaining the dominant real part, prefer a
        // real one: "the dominant eigenvalue is real" is a statement about the
        // spectrum, not about which tied eigenvalue the scan happened to keep. This
        // only widens the accepted set, and only in the tied case. Matches MATLAB
        // BUTools/ph/CheckMERepresentation.m and native Python butools.ph.check.
        double domRe = Math.abs(maxEv.getReal());
        double domTol = Math.max(prec, 1e-8 * domRe);
        for (Complex ev : eigenvalues) {
            if (Math.abs(Math.abs(ev.getReal()) - domRe) <= domTol
                    && Math.abs(ev.getImaginary()) <= domTol) {
                maxEv = ev;
                break;
            }
        }

        // Check if dominant eigenvalue is real
        if (Math.abs(maxEv.getImaginary()) > prec) {
            return false;
        }

        return true;
    }

    public static boolean checkMERepresentation(Matrix alpha, Matrix A) {
        return checkMERepresentation(alpha, A, 1e-14);
    }

    /**
     * Overload for double[] alpha.
     */
    public static boolean checkMERepresentation(double[] alpha, Matrix A, double prec) {
        return checkMERepresentation(new Matrix(alpha), A, prec);
    }

    public static boolean checkMERepresentation(double[] alpha, Matrix A) {
        return checkMERepresentation(new Matrix(alpha), A, 1e-14);
    }
}
