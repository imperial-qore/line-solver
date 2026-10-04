/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.map;

import jline.util.matrix.Matrix;

public final class LagkJointMomentsFromMAP {
    private LagkJointMomentsFromMAP() {}

    /**
     * Returns the lag-L joint moments of a continuous Markovian arrival process.
     *
     * @param D0 The D0 matrix of the Markovian arrival process
     * @param D1 The D1 matrix of the Markovian arrival process
     * @param K The dimension of the matrix of joint moments to compute.
     *          If K=0, the MxM joint moments will be computed.
     * @param L The lag at which the joint moments are computed. Default is 1.
     * @param prec Numerical precision to check if the input is valid.
     * @return Matrix containing the lag-L joint moments
     */
    public static Matrix lagkJointMomentsFromMAP(Matrix D0, Matrix D1, int K, int L, double prec) {
        if (!CheckMAPRepresentation.checkMAPRepresentation(D0, D1, prec)) {
            throw new IllegalArgumentException("LagkJointMomentsFromMAP: Input isn't a valid MAP representation!");
        }

        return LagkJointMomentsFromRAP.lagkJointMomentsFromRAP(D0, D1, K, L, prec);
    }

    public static Matrix lagkJointMomentsFromMAP(Matrix D0, Matrix D1, int K, int L) {
        return lagkJointMomentsFromMAP(D0, D1, K, L, 1e-14);
    }

    public static Matrix lagkJointMomentsFromMAP(Matrix D0, Matrix D1, int K) {
        return lagkJointMomentsFromMAP(D0, D1, K, 1, 1e-14);
    }

    public static Matrix lagkJointMomentsFromMAP(Matrix D0, Matrix D1) {
        return lagkJointMomentsFromMAP(D0, D1, 0, 1, 1e-14);
    }
}
