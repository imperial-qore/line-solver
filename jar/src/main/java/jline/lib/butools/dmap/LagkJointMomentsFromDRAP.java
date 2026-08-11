/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dmap;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class LagkJointMomentsFromDRAP {
    private LagkJointMomentsFromDRAP() {}

    /**
     * Returns the lag-L joint moments of a discrete rational arrival process.
     *
     * @param H0 The H0 matrix of the discrete rational arrival process
     * @param H1 The H1 matrix of the discrete rational arrival process
     * @param K The dimension of the matrix of joint moments to compute.
     *          If K=0, the MxM joint moments will be computed.
     * @param L The lag at which the joint moments are computed. Default is 1.
     * @param prec Numerical precision to check if the input is valid.
     * @return Matrix containing the lag-L joint moments
     */
    public static Matrix lagkJointMomentsFromDRAP(Matrix H0, Matrix H1, int K, int L, double prec) {
        if (!CheckDRAPRepresentation.checkDRAPRepresentation(H0, H1, prec)) {
            throw new IllegalArgumentException("LagkJointMomentsFromDRAP: Input isn't a valid DRAP representation!");
        }

        MatrixCell H = new MatrixCell(2);
        H.set(0, H0);
        H.set(1, H1);

        int actualK = (K == 0) ? H0.getNumRows() - 1 : K;

        MatrixCell mom = LagkJointMomentsFromDMRAP.lagkJointMomentsFromDMRAP(H, actualK, L, prec);
        return mom.get(0);
    }

    public static Matrix lagkJointMomentsFromDRAP(Matrix H0, Matrix H1, int K, int L) {
        return lagkJointMomentsFromDRAP(H0, H1, K, L, 1e-14);
    }

    public static Matrix lagkJointMomentsFromDRAP(Matrix H0, Matrix H1, int K) {
        return lagkJointMomentsFromDRAP(H0, H1, K, 1, 1e-14);
    }

    public static Matrix lagkJointMomentsFromDRAP(Matrix H0, Matrix H1) {
        return lagkJointMomentsFromDRAP(H0, H1, 0, 1, 1e-14);
    }
}
