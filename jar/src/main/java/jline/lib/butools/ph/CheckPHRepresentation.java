/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import static jline.lib.butools.mc.CheckGenerator.checkGenerator;
import static jline.lib.butools.mc.CheckProbVector.checkProbVector;

import jline.util.matrix.Matrix;

public final class CheckPHRepresentation {
    private CheckPHRepresentation() {}

    /**
     * Checks if the given vector and matrix define a valid phase-
     * type representation.
     *
     * @param alpha Initial vector of the phase-type distribution to check
     * @param A Transient generator of the phase-type distribution to check
     * @param prec Numerical precision. The default value is 1e-14.
     * @return True, if vector alpha is a probability vector and matrix
     *         A is a transient generator, and they have the same size.
     */
    public static boolean checkPHRepresentation(Matrix alpha, Matrix A, double prec) {
        // Check if vector and matrix have same size
        if (alpha.length() != A.getNumRows()) {
            return false;
        }

        // Check if A is a valid transient generator
        if (!checkGenerator(A, true, prec)) {
            return false;
        }

        // Check if alpha is a valid substochastic probability vector
        if (!checkProbVector(alpha, true, prec)) {
            return false;
        }

        return true;
    }

    public static boolean checkPHRepresentation(Matrix alpha, Matrix A) {
        return checkPHRepresentation(alpha, A, 1e-14);
    }

    /**
     * Overload for double[] alpha.
     */
    public static boolean checkPHRepresentation(double[] alpha, Matrix A, double prec) {
        return checkPHRepresentation(new Matrix(alpha), A, prec);
    }

    public static boolean checkPHRepresentation(double[] alpha, Matrix A) {
        return checkPHRepresentation(new Matrix(alpha), A, 1e-14);
    }
}
