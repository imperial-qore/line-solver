/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class LagkJointMomentsFromDMMAP {
    private LagkJointMomentsFromDMMAP() {}

    /**
     * Returns the lag-L joint moments of a discrete marked Markovian arrival process.
     *
     * @param D The D0...DN matrices of the DMMAP (as MatrixCell)
     * @param K The dimension of the matrix of joint moments to compute.
     *          If K=0, the MxM joint moments will be computed.
     * @param L The lag at which the joint moments are computed. Default is 1.
     * @param prec Numerical precision to check if the input is valid.
     * @return List of matrices containing the lag-L joint moments
     */
    public static MatrixCell lagkJointMomentsFromDMMAP(MatrixCell D, int K, int L, double prec) {
        if (!CheckDMMAPRepresentation.checkDMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("LagkJointMomentsFromDMMAP: Input isn't a valid DMMAP representation!");
        }

        int actualK = (K == 0) ? D.get(0).getNumRows() - 1 : K;

        return LagkJointMomentsFromDMRAP.lagkJointMomentsFromDMRAP(D, actualK, L, prec);
    }

    public static MatrixCell lagkJointMomentsFromDMMAP(MatrixCell D) {
        return lagkJointMomentsFromDMMAP(D, 0, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromDMMAP(MatrixCell D, int K) {
        return lagkJointMomentsFromDMMAP(D, K, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromDMMAP(MatrixCell D, int K, int L) {
        return lagkJointMomentsFromDMMAP(D, K, L, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static MatrixCell lagkJointMomentsFromDMMAP(Matrix[] D, int K, int L, double prec) {
        MatrixCell cell = new MatrixCell(D.length);
        for (int i = 0; i < D.length; i++) {
            cell.set(i, D[i]);
        }
        return lagkJointMomentsFromDMMAP(cell, K, L, prec);
    }

    public static MatrixCell lagkJointMomentsFromDMMAP(Matrix[] D) {
        return lagkJointMomentsFromDMMAP(D, 0, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromDMMAP(Matrix[] D, int K) {
        return lagkJointMomentsFromDMMAP(D, K, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromDMMAP(Matrix[] D, int K, int L) {
        return lagkJointMomentsFromDMMAP(D, K, L, 1e-14);
    }
}
