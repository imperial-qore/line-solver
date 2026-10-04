/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.Matrix;

public final class CheckMGRepresentation {
    private CheckMGRepresentation() {}

    /**
     * Checks if the given vector and matrix define a valid matrix-geometric representation.
     *
     * @param alpha Initial vector of the matrix-geometric distribution to check
     * @param A Matrix parameter of the matrix-geometric distribution to check
     * @param prec Numerical precision.
     * @return True if the matrix is a square matrix, the vector and the matrix have the same size,
     *         the dominant eigenvalue is positive, less than 1 and real.
     *
     * Note: This procedure does not check the positivity of the density!
     */
    public static boolean checkMGRepresentation(Matrix alpha, Matrix A, double prec) {
        // Check if matrix is square
        if (A.getNumRows() != A.getNumCols()) {
            return false;
        }

        // Check if vector and matrix have matching sizes
        if (alpha.length() != A.getNumRows()) {
            return false;
        }

        // Check if sum of alpha is between 0 and 1
        double alphaSum = alpha.elementSum();
        if (alphaSum < -prec || alphaSum > 1 + prec) {
            return false;
        }

        // Get eigenvalues and sort by absolute value (descending)
        List<Complex> eigenvalues = A.eig();
        Collections.sort(eigenvalues, new Comparator<Complex>() {
            @Override
            public int compare(Complex a, Complex b) {
                double absA = Math.abs(a.getReal() * a.getReal() + a.getImaginary() * a.getImaginary());
                double absB = Math.abs(b.getReal() * b.getReal() + b.getImaginary() * b.getImaginary());
                return Double.compare(absB, absA);
            }
        });

        if (eigenvalues.isEmpty()) {
            return false;
        }

        Complex maxEv = eigenvalues.get(0);

        // Check if largest eigenvalue is real
        if (Math.abs(maxEv.getImaginary()) > prec) {
            return false;
        }

        // Check if largest eigenvalue is not greater than 1
        if (maxEv.getReal() > 1 + prec) {
            return false;
        }

        return true;
    }

    public static boolean checkMGRepresentation(Matrix alpha, Matrix A) {
        return checkMGRepresentation(alpha, A, 1e-14);
    }

    public static boolean checkMGRepresentation(double[] alpha, Matrix A, double prec) {
        return checkMGRepresentation(new Matrix(alpha), A, prec);
    }

    public static boolean checkMGRepresentation(double[] alpha, Matrix A) {
        return checkMGRepresentation(new Matrix(alpha), A, 1e-14);
    }
}
