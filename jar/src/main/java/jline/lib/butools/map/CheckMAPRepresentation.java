/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.lib.butools.mc.CheckGenerator;
import jline.util.matrix.Matrix;

public final class CheckMAPRepresentation {
    private CheckMAPRepresentation() {}

    /**
     * Checks if the input matrices define a continuous time MAP.
     */
    public static boolean checkMAPRepresentation(Matrix D0, Matrix D1, double prec) {
        // Check if D0 is a transient generator
        if (!CheckGenerator.checkGenerator(D0, true, prec)) {
            return false;
        }

        // Check if D0 and D1 have same size
        if (D0.getNumRows() != D1.getNumRows() || D0.getNumCols() != D1.getNumCols()) {
            return false;
        }

        int n = D0.getNumRows();

        // Check if D1 has non-negative elements
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (D1.get(i, j) < -prec) {
                    return false;
                }
            }
        }

        // Check if rowsum of D0+D1 is 0
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += D0.get(i, j) + D1.get(i, j);
            }
            if (Math.abs(rowSum) > prec) {
                return false;
            }
        }

        return true;
    }

    public static boolean checkMAPRepresentation(Matrix D0, Matrix D1) {
        return checkMAPRepresentation(D0, D1, 1e-12);
    }

    /**
     * Checks if the input matrices define a valid RAP representation.
     */
    public static boolean checkRAPRepresentation(Matrix H0, Matrix H1, double prec) {
        if (H0.getNumRows() != H0.getNumCols() || H1.getNumRows() != H1.getNumCols()) {
            return false;
        }
        if (H0.getNumRows() != H1.getNumRows()) {
            return false;
        }

        int n = H0.getNumRows();

        // Check if rowsum of H0+H1 is 0
        for (int i = 0; i < n; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < n; j++) {
                rowSum += H0.get(i, j) + H1.get(i, j);
            }
            if (Math.abs(rowSum) > prec) {
                return false;
            }
        }

        // Check if dominant eigenvalue of H0 is negative and real
        List<Complex> eigenvalues = H0.eig();
        double maxRealEv = Double.NEGATIVE_INFINITY;
        for (Complex ev : eigenvalues) {
            if (ev.getReal() > maxRealEv) {
                maxRealEv = ev.getReal();
            }
        }

        return !(maxRealEv >= prec);
    }

    public static boolean checkRAPRepresentation(Matrix H0, Matrix H1) {
        return checkRAPRepresentation(H0, H1, 1e-12);
    }
}
