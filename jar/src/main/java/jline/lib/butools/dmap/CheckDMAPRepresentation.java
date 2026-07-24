/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.lib.butools.mc.CheckProbMatrix;
import jline.util.matrix.Matrix;

public final class CheckDMAPRepresentation {
    private CheckDMAPRepresentation() {}

    /**
     * Checks if the input matrices define a discrete time MAP.
     *
     * Matrices D0 and D1 must have the same size, D0 must be a transient
     * probability matrix, D1 has only non-negative elements, and the rowsum
     * of D0+D1 is 1 (up to the numerical precision).
     *
     * @param D0 The D0 matrix of the DMAP to check
     * @param D1 The D1 matrix of the DMAP to check
     * @param prec Numerical precision, default value is 1e-14
     * @return True if the matrices define a valid DMAP
     */
    public static boolean checkDMAPRepresentation(Matrix D0, Matrix D1, double prec) {
        // Check if D0 is a transient probability matrix
        if (!CheckProbMatrix.checkProbMatrix(D0, true, prec)) {
            return false;
        }

        // Check if D0 and D1 have the same size
        if (D0.getNumRows() != D1.getNumRows() || D0.getNumCols() != D1.getNumCols()) {
            return false;
        }

        // Check for negative elements
        if (D0.elementMin() < -prec || D1.elementMin() < -prec) {
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

        return true;
    }

    public static boolean checkDMAPRepresentation(Matrix D0, Matrix D1) {
        return checkDMAPRepresentation(D0, D1, 1e-14);
    }
}
