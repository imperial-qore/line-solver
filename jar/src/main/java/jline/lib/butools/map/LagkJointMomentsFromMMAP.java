/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class LagkJointMomentsFromMMAP {
    private LagkJointMomentsFromMMAP() {}

    /**
     * Returns the lag-L joint moments of a continuous marked Markovian arrival process.
     *
     * @param D The D0...DN matrices of the MMAP (as MatrixCell)
     * @param K The dimension of the matrix of joint moments to compute.
     *          If K=0, the MxM joint moments will be computed.
     * @param L The lag at which the joint moments are computed. Default is 1.
     * @param prec Numerical precision to check if the input is valid.
     * @return List of matrices containing the lag-L joint moments
     */
    public static MatrixCell lagkJointMomentsFromMMAP(MatrixCell D, int K, int L, double prec) {
        if (!CheckMMAPRepresentation.checkMMAPRepresentation(D, prec)) {
            throw new IllegalArgumentException("LagkJointMomentsFromMMAP: Input isn't a valid MMAP representation!");
        }

        int actualK = (K == 0) ? D.get(0).getNumRows() - 1 : K;

        return LagkJointMomentsFromMRAP.lagkJointMomentsFromMRAP(D, actualK, L, prec);
    }

    public static MatrixCell lagkJointMomentsFromMMAP(MatrixCell D) {
        return lagkJointMomentsFromMMAP(D, 0, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromMMAP(MatrixCell D, int K) {
        return lagkJointMomentsFromMMAP(D, K, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromMMAP(MatrixCell D, int K, int L) {
        return lagkJointMomentsFromMMAP(D, K, L, 1e-14);
    }

    /**
     * Overload for Matrix[].
     */
    public static MatrixCell lagkJointMomentsFromMMAP(Matrix[] D, int K, int L, double prec) {
        MatrixCell cell = new MatrixCell(D.length);
        for (int i = 0; i < D.length; i++) {
            cell.set(i, D[i]);
        }
        return lagkJointMomentsFromMMAP(cell, K, L, prec);
    }

    public static MatrixCell lagkJointMomentsFromMMAP(Matrix[] D) {
        return lagkJointMomentsFromMMAP(D, 0, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromMMAP(Matrix[] D, int K) {
        return lagkJointMomentsFromMMAP(D, K, 1, 1e-14);
    }

    public static MatrixCell lagkJointMomentsFromMMAP(Matrix[] D, int K, int L) {
        return lagkJointMomentsFromMMAP(D, K, L, 1e-14);
    }
}
