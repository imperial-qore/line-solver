/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc;

import java.util.List;

import org.apache.commons.math3.complex.Complex;

import jline.util.matrix.Matrix;

public final class CheckGenerator {
    private CheckGenerator() {}

    /**
     * Checks if the matrix is a valid generator matrix.
     *
     * @param Q         The generator to check.
     * @param transient_ If true, checks if Q is a transient generator.
     * @param prec      Precision threshold (default 1e-14).
     * @return The result of the check.
     */
    public static boolean checkGenerator(Matrix Q, boolean transient_, double prec) {
        if (Q.getNumRows() != Q.getNumCols()) {
            return false;
        }

        int N = Q.getNumRows();

        for (int i = 0; i < N; i++) {
            if (Q.get(i, i) >= prec) {
                return false;
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != j && Q.get(i, j) < -prec) {
                    return false;
                }
            }
        }

        if (transient_) {
            for (int i = 0; i < N; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < N; j++) {
                    rowSum += Q.get(i, j);
                }
                if (rowSum > prec) {
                    return false;
                }
            }

            List<Complex> eigenvalues = Q.eig();
            for (Complex ev : eigenvalues) {
                if (ev.getReal() >= prec) {
                    return false;
                }
            }
        } else {
            for (int i = 0; i < N; i++) {
                double rowSum = 0.0;
                for (int j = 0; j < N; j++) {
                    rowSum += Q.get(i, j);
                }
                if (Math.abs(rowSum) > prec) {
                    return false;
                }
            }
        }

        return true;
    }

    public static boolean checkGenerator(Matrix Q, boolean transient_) {
        return checkGenerator(Q, transient_, 1e-14);
    }

    public static boolean checkGenerator(Matrix Q) {
        return checkGenerator(Q, false, 1e-14);
    }
}
