/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import jline.lib.butools.mc.CheckProbMatrix;
import jline.lib.butools.mc.CheckProbVector;
import jline.util.matrix.Matrix;

public final class CheckDPHRepresentation {
    private CheckDPHRepresentation() {}

    /**
     * Checks if the given vector and matrix define a valid discrete phase-type representation.
     *
     * @param alpha Initial vector of the phase-type distribution to check
     * @param A Transient generator of the phase-type distribution to check
     * @param prec Numerical precision.
     * @return True if vector alpha is a probability vector and matrix A is substochastic,
     *         and they have the same size.
     */
    public static boolean checkDPHRepresentation(Matrix alpha, Matrix A, double prec) {
        // Check if alpha and A have the same size
        if (alpha.length() != A.getNumRows()) {
            return false;
        }

        // Check if A is a valid transient (substochastic) probability matrix
        if (!CheckProbMatrix.checkProbMatrix(A, true, prec)) {
            return false;
        }

        // Check if alpha is a valid substochastic probability vector
        if (!CheckProbVector.checkProbVector(alpha, true, prec)) {
            return false;
        }

        return true;
    }

    public static boolean checkDPHRepresentation(Matrix alpha, Matrix A) {
        return checkDPHRepresentation(alpha, A, 1e-14);
    }

    /**
     * Overload for double[] alpha.
     */
    public static boolean checkDPHRepresentation(double[] alpha, Matrix A, double prec) {
        return checkDPHRepresentation(new Matrix(alpha), A, prec);
    }

    public static boolean checkDPHRepresentation(double[] alpha, Matrix A) {
        return checkDPHRepresentation(new Matrix(alpha), A, 1e-14);
    }
}
