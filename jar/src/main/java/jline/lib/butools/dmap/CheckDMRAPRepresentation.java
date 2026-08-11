/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class CheckDMRAPRepresentation {
    private CheckDMRAPRepresentation() {}

    /**
     * Checks if the input matrices define a discrete time MRAP.
     *
     * All matrices H0...HK must have the same size, the dominant eigenvalue
     * of H0 is real and less than 1, and the rowsum of H0+H1+...+HK is 1
     * (up to the numerical precision).
     *
     * @param H The H0...HK matrices of the DMRAP to check (as MatrixCell)
     * @param prec Numerical precision, default value is 1e-14
     * @return True if the matrices define a valid DMRAP
     */
    public static boolean checkDMRAPRepresentation(MatrixCell H, double prec) {
        if (H.size() < 2) {
            return false;
        }

        // Sum H1...HK
        Matrix sumH = H.get(1).copy();
        for (int i = 2; i < H.size(); i++) {
            sumH = sumH.add(H.get(i));
        }

        return CheckDRAPRepresentation.checkDRAPRepresentation(H.get(0), sumH, prec);
    }

    public static boolean checkDMRAPRepresentation(MatrixCell H) {
        return checkDMRAPRepresentation(H, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static boolean checkDMRAPRepresentation(Matrix[] H, double prec) {
        if (H.length < 2) {
            return false;
        }

        // Sum H1...HK
        Matrix sumH = H[1].copy();
        for (int i = 2; i < H.length; i++) {
            sumH = sumH.add(H[i]);
        }

        return CheckDRAPRepresentation.checkDRAPRepresentation(H[0], sumH, prec);
    }

    public static boolean checkDMRAPRepresentation(Matrix[] H) {
        return checkDMRAPRepresentation(H, 1e-14);
    }
}
