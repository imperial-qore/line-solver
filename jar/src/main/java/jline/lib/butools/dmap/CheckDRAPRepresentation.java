/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.complex.Complex;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public final class CheckDRAPRepresentation {
    private CheckDRAPRepresentation() {}

    /**
     * Checks if the input matrices define a discrete time RAP.
     *
     * @param D0 The D0 matrix of the DRAP to check
     * @param D1 The D1 matrix of the DRAP to check
     * @param prec Numerical precision
     * @return True if the matrices define a valid DRAP
     */
    public static boolean checkDRAPRepresentation(Matrix D0, Matrix D1, double prec) {
        // Check if D0 and D1 have the same size and D0 is square
        if (D0.getNumRows() != D1.getNumRows() || D0.getNumCols() != D1.getNumCols()
                || D0.getNumRows() != D0.getNumCols()) {
            return false;
        }

        // Check if rowsums of D0+D1 equal 1
        Matrix D0D1 = D0.add(D1);
        for (int i = 0; i < D0D1.getNumRows(); i++) {
            double rowSum = 0.0;
            for (int j = 0; j < D0D1.getNumCols(); j++) {
                rowSum += D0D1.get(i, j);
            }
            if (Math.abs(rowSum - 1.0) > prec) {
                return false;
            }
        }

        // Get eigenvalues and sort by absolute value (descending)
        List<Complex> eigenvalues = D0.eig();
        Collections.sort(eigenvalues, new Comparator<Complex>() {
            @Override
            public int compare(Complex a, Complex b) {
                double absA = Math.abs(a.getReal() * a.getReal() + a.getImaginary() * a.getImaginary());
                double absB = Math.abs(b.getReal() * b.getReal() + b.getImaginary() * b.getImaginary());
                return Double.compare(absB, absA); // descending
            }
        });

        if (eigenvalues.isEmpty()) {
            return false;
        }

        Complex maxEv = eigenvalues.get(0);

        // Check if dominant eigenvalue is real
        if (Math.abs(maxEv.getImaginary()) > prec) {
            return false;
        }

        // Check if dominant eigenvalue is not greater than 1
        if (maxEv.getReal() > 1 + prec) {
            return false;
        }

        return true;
    }

    public static boolean checkDRAPRepresentation(Matrix D0, Matrix D1) {
        return checkDRAPRepresentation(D0, D1, 1e-14);
    }
}
