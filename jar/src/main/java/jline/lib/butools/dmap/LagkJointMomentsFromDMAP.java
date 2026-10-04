/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class LagkJointMomentsFromDMAP {
    private LagkJointMomentsFromDMAP() {}

    /**
     * Returns the lag-L joint moments of a discrete Markovian arrival process.
     *
     * @param D0 The D0 matrix of the discrete Markovian arrival process
     * @param D1 The D1 matrix of the discrete Markovian arrival process
     * @param K The dimension of the matrix of joint moments to compute.
     *          If K=0, the MxM joint moments will be computed.
     * @param L The lag at which the joint moments are computed. Default is 1.
     * @param prec Numerical precision to check if the input is valid.
     * @return Matrix containing the lag-L joint moments
     */
    public static Matrix lagkJointMomentsFromDMAP(Matrix D0, Matrix D1, int K, int L, double prec) {
        if (!CheckDMAPRepresentation.checkDMAPRepresentation(D0, D1, prec)) {
            throw new IllegalArgumentException("LagkJointMomentsFromDMAP: Input isn't a valid DMAP representation!");
        }

        MatrixCell H = new MatrixCell(2);
        H.set(0, D0);
        H.set(1, D1);

        int actualK = (K == 0) ? D0.getNumRows() - 1 : K;

        MatrixCell mom = LagkJointMomentsFromDMRAP.lagkJointMomentsFromDMRAP(H, actualK, L, prec);
        return mom.get(0);
    }

    public static Matrix lagkJointMomentsFromDMAP(Matrix D0, Matrix D1, int K, int L) {
        return lagkJointMomentsFromDMAP(D0, D1, K, L, 1e-14);
    }

    public static Matrix lagkJointMomentsFromDMAP(Matrix D0, Matrix D1, int K) {
        return lagkJointMomentsFromDMAP(D0, D1, K, 1, 1e-14);
    }

    public static Matrix lagkJointMomentsFromDMAP(Matrix D0, Matrix D1) {
        return lagkJointMomentsFromDMAP(D0, D1, 0, 1, 1e-14);
    }
}
