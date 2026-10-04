/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.util.matrix.Matrix;

public final class MomentsFromME {
    private MomentsFromME() {}

    /**
     * Returns the first K moments of a matrix-exponential distribution.
     *
     * @param alpha The initial vector of the matrix-exponential distribution.
     * @param A The matrix parameter of the matrix-exponential distribution.
     * @param K Number of moments to compute. If K=0, 2*M-1 moments are computed.
     * @return The vector of moments.
     */
    public static double[] momentsFromME(Matrix alpha, Matrix A, int K) {
        int numMoments = K == 0 ? 2 * alpha.length() - 1 : K;

        double[] moms = new double[numMoments];
        Matrix iA = A.neg().inv();

        Matrix iAPower = iA.copy();
        for (int i = 0; i < numMoments; i++) {
            long factorial = 1L;
            for (int n = 1; n <= i + 1; n++) factorial *= n;
            moms[i] = (double) factorial * alpha.mult(iAPower).elementSum();
            iAPower = iAPower.mult(iA);
        }

        return moms;
    }

    public static double[] momentsFromME(Matrix alpha, Matrix A) {
        return momentsFromME(alpha, A, 0);
    }

    /**
     * Returns the first K moments of a phase-type distribution.
     */
    public static double[] momentsFromPH(Matrix alpha, Matrix A, int K) {
        return momentsFromME(alpha, A, K);
    }

    public static double[] momentsFromPH(Matrix alpha, Matrix A) {
        return momentsFromME(alpha, A, 0);
    }

    /**
     * Overload for double[] alpha.
     */
    public static double[] momentsFromME(double[] alpha, Matrix A, int K) {
        return momentsFromME(new Matrix(alpha), A, K);
    }

    public static double[] momentsFromME(double[] alpha, Matrix A) {
        return momentsFromME(new Matrix(alpha), A, 0);
    }

    /**
     * Overload for double[] alpha.
     */
    public static double[] momentsFromPH(double[] alpha, Matrix A, int K) {
        return momentsFromME(new Matrix(alpha), A, K);
    }

    public static double[] momentsFromPH(double[] alpha, Matrix A) {
        return momentsFromME(new Matrix(alpha), A, 0);
    }
}
