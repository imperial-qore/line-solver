/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 *
 * Reference:
 * P. Buchholz, M. Telek, "On minimal representation of rational arrival
 * processes." Madrid Conference on Queueing theory (MCQT), June 2010.
 */
package jline.lib.butools.map;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class MinimalRepFromRAP {
    private MinimalRepFromRAP() {}

    /**
     * Returns the minimal representation of a rational arrival process.
     *
     * @param H0 The H0 matrix of the rational arrival process
     * @param H1 The H1 matrix of the rational arrival process
     * @param how Determines how the representation is minimized:
     *            "cont" = controllability, "obs" = observability,
     *            "obscont" = both (default)
     * @param precision Precision used by the Staircase algorithm (default 1e-12)
     * @return Pair of (D0, D1) matrices of the minimal representation
     */
    public static Pair<Matrix, Matrix> minimalRepFromRAP(Matrix H0, Matrix H1, String how, double precision) {
        MatrixCell H = new MatrixCell(2);
        H.set(0, H0);
        H.set(1, H1);
        MatrixCell result = MinimalRepFromMRAP.minimalRepFromMRAP(H, how, precision);
        return new Pair<Matrix, Matrix>(result.get(0), result.get(1));
    }

    public static Pair<Matrix, Matrix> minimalRepFromRAP(Matrix H0, Matrix H1, String how) {
        return minimalRepFromRAP(H0, H1, how, 1e-12);
    }

    public static Pair<Matrix, Matrix> minimalRepFromRAP(Matrix H0, Matrix H1) {
        return minimalRepFromRAP(H0, H1, "obscont", 1e-12);
    }
}
