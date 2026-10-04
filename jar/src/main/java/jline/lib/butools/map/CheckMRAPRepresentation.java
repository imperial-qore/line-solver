/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.lib.butools.ph.CheckRAPRepresentation;

public final class CheckMRAPRepresentation {
    private CheckMRAPRepresentation() {}

    /**
     * Checks if the input matrices define a continuous time MRAP.
     *
     * All matrices H0...HK must have the same size, the dominant eigenvalue
     * of H0 has negative real part, and the rowsum of H0+H1+...+HK is 0
     * (up to the numerical precision).
     *
     * @param H The H0...HK matrices of the MRAP to check (as MatrixCell)
     * @param prec Numerical precision, default value is 1e-14
     * @return True if the matrices define a valid MRAP
     */
    public static boolean checkMRAPRepresentation(MatrixCell H, double prec) {
        if (H.size() < 2) {
            return false;
        }

        // Sum H1...HK
        Matrix sumH = H.get(1).copy();
        for (int i = 2; i < H.size(); i++) {
            sumH = sumH.add(H.get(i));
        }

        return CheckRAPRepresentation.checkRAPRepresentation(H.get(0), sumH, prec);
    }

    public static boolean checkMRAPRepresentation(MatrixCell H) {
        return checkMRAPRepresentation(H, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static boolean checkMRAPRepresentation(Matrix[] H, double prec) {
        if (H.length < 2) {
            return false;
        }

        // Sum H1...HK
        Matrix sumH = H[1].copy();
        for (int i = 2; i < H.length; i++) {
            sumH = sumH.add(H[i]);
        }

        return CheckRAPRepresentation.checkRAPRepresentation(H[0], sumH, prec);
    }

    public static boolean checkMRAPRepresentation(Matrix[] H) {
        return checkMRAPRepresentation(H, 1e-14);
    }
}
