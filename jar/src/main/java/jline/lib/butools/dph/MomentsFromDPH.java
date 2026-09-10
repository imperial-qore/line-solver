/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import jline.util.matrix.Matrix;

public final class MomentsFromDPH {
    private MomentsFromDPH() {}

    /**
     * Returns the first K moments of a discrete phase-type distribution.
     *
     * @param alpha The initial probability vector of the discrete phase-type distribution.
     *        The sum of the entries of alpha is less or equal to 1.
     * @param A The transient generator matrix of the discrete phase-type distribution.
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     *        The default value is 0.
     * @return The vector of moments.
     */
    public static double[] momentsFromDPH(Matrix alpha, Matrix A, int K) {
        return MomentsFromMG.momentsFromMG(alpha, A, K);
    }

    public static double[] momentsFromDPH(Matrix alpha, Matrix A) {
        return momentsFromDPH(alpha, A, 0);
    }

    /**
     * Overload for double[] alpha.
     */
    public static double[] momentsFromDPH(double[] alpha, Matrix A, int K) {
        return momentsFromDPH(new Matrix(alpha), A, K);
    }

    public static double[] momentsFromDPH(double[] alpha, Matrix A) {
        return momentsFromDPH(alpha, A, 0);
    }
}
