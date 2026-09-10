/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class DMAPFromDRAP {
    private DMAPFromDRAP() {}

    /**
     * Obtains a Markovian representation of a discrete rational
     * arrival process of the same size, if possible.
     *
     * @param H0 The H0 matrix of the discrete rational arrival process
     * @param H1 The H1 matrix of the discrete rational arrival process
     * @param prec A representation is considered to be Markovian if it is closer than this precision
     * @return Pair of (D0, D1) matrices of the discrete Markovian arrival process
     */
    public static Pair<Matrix, Matrix> dmapFromDRAP(Matrix H0, Matrix H1, double prec) {
        if (!CheckDRAPRepresentation.checkDRAPRepresentation(H0, H1, prec)) {
            throw new IllegalArgumentException("DMAPFromDRAP: Input isn't a valid DRAP representation!");
        }

        MatrixCell H = new MatrixCell(2);
        H.set(0, H0);
        H.set(1, H1);

        MatrixCell Y = DMMAPFromDMRAP.dmmapFromDMRAP(H, prec);

        return new Pair<Matrix, Matrix>(Y.get(0), Y.get(1));
    }

    public static Pair<Matrix, Matrix> dmapFromDRAP(Matrix H0, Matrix H1) {
        return dmapFromDRAP(H0, H1, 1e-14);
    }
}
