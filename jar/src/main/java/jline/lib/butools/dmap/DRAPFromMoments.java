/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class DRAPFromMoments {
    private DRAPFromMoments() {}

    /**
     * Creates a discrete rational arrival process that has the
     * same marginal and lag-1 joint moments as given.
     *
     * @param moms The list of marginal moments. To obtain a rational
     *             process of order M, 2*M-1 marginal moments are required.
     * @param Nm The matrix of lag-1 joint moments
     * @return Pair of (H0, H1) matrices of the discrete rational process
     */
    public static Pair<Matrix, Matrix> drapFromMoments(double[] moms, Matrix Nm) {
        MatrixCell NmCell = new MatrixCell(1);
        NmCell.set(0, Nm);

        MatrixCell H = DMRAPFromMoments.dmrapFromMoments(moms, NmCell);

        return new Pair<>(H.get(0), H.get(1));
    }
}
