/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class CheckMMAPRepresentation {
    private CheckMMAPRepresentation() {}

    /**
     * Checks if the input matrices define a continuous time MMAP.
     *
     * All matrices D0...DK must have the same size, D0 must be a
     * transient generator matrix, D1...DK have only non-negative
     * elements, and the rowsum of D0+D1+...+DK is 0 (up to the
     * numerical precision).
     *
     * @param D The D0...DK matrices of the MMAP to check (as MatrixCell)
     * @param prec Numerical precision
     * @return True if the matrices define a valid MMAP
     */
    public static boolean checkMMAPRepresentation(MatrixCell D, double prec) {
        if (D.size() < 2) {
            return false;
        }

        for (int i = 1; i < D.size(); i++) {
            if (D.get(i).elementMin() < -prec) {
                return false;
            }
        }

        Matrix sumD = D.get(1).copy();
        for (int i = 2; i < D.size(); i++) {
            sumD = sumD.add(D.get(i));
        }

        return CheckMAPRepresentation.checkMAPRepresentation(D.get(0), sumD, prec);
    }

    public static boolean checkMMAPRepresentation(MatrixCell D) {
        return checkMMAPRepresentation(D, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static boolean checkMMAPRepresentation(Matrix[] D, double prec) {
        if (D.length < 2) {
            return false;
        }

        for (int i = 1; i < D.length; i++) {
            if (D[i].elementMin() < -prec) {
                return false;
            }
        }

        Matrix sumD = D[1].copy();
        for (int i = 2; i < D.length; i++) {
            sumD = sumD.add(D[i]);
        }

        return CheckMAPRepresentation.checkMAPRepresentation(D[0], sumD, prec);
    }

    public static boolean checkMMAPRepresentation(Matrix[] D) {
        return checkMMAPRepresentation(D, 1e-14);
    }
}
