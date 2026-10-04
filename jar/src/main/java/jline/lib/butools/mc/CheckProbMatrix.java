/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc;

import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.Matrix;

public final class CheckProbMatrix {
    private CheckProbMatrix() {}

    /**
     * Checks if the matrix is a valid probability matrix.
     *
     * @param P         The matrix to check.
     * @param transient_ If true, checks if P is a transient probability matrix.
     * @param prec      Numerical precision for zero comparisons.
     * @return The result of the check.
     */
    public static boolean checkProbMatrix(Matrix P, boolean transient_, double prec) {
        if (P.getNumRows() != P.getNumCols()) {
            return false;
        }

        int N = P.getNumRows();

        // Check for negative elements
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (P.get(i, j) < -prec) {
                    return false;
                }
            }
        }

        if (transient_) {
            for (int i = 0; i < N; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < N; j++) {
                    rowSum += P.get(i, j);
                }
                if (rowSum - 1.0 > N * prec) {
                    return false;
                }
            }

            List<Complex> eigenvalues = P.eig();
            for (Complex ev : eigenvalues) {
                if (ev.getReal() >= 1.0 - prec) {
                    return false;
                }
            }
        } else {
            for (int i = 0; i < N; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < N; j++) {
                    rowSum += P.get(i, j);
                }
                if (Math.abs(rowSum - 1.0) > N * prec) {
                    return false;
                }
            }
        }

        return true;
    }

    public static boolean checkProbMatrix(Matrix P, boolean transient_) {
        return checkProbMatrix(P, transient_, 1e-14);
    }

    public static boolean checkProbMatrix(Matrix P) {
        return checkProbMatrix(P, false, 1e-14);
    }
}
