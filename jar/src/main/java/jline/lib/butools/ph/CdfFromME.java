/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.util.matrix.Matrix;

public final class CdfFromME {
    private CdfFromME() {}

    /**
     * Returns the cumulative distribution function of a
     * matrix-exponential distribution.
     *
     * @param alpha The initial vector of the matrix-exponential distribution.
     * @param A The matrix parameter of the matrix-exponential distribution.
     * @param x The points at which the CDF will be computed.
     * @return The values of the CDF at the corresponding "x" values.
     */
    public static double[] cdfFromME(Matrix alpha, Matrix A, double[] x) {
        double[] cdf = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            if (x[i] < 0) {
                cdf[i] = 0.0;
            } else {
                Matrix expAt = A.scale(x[i]).expm();
                cdf[i] = 1.0 - alpha.mult(expAt).elementSum();
            }
        }
        return cdf;
    }

    /**
     * Returns the cumulative distribution function of a
     * phase-type distribution.
     */
    public static double[] cdfFromPH(Matrix alpha, Matrix A, double[] x) {
        return cdfFromME(alpha, A, x);
    }

    /**
     * Overload for double[] alpha.
     */
    public static double[] cdfFromME(double[] alpha, Matrix A, double[] x) {
        return cdfFromME(new Matrix(alpha), A, x);
    }

    /**
     * Overload for double[] alpha.
     */
    public static double[] cdfFromPH(double[] alpha, Matrix A, double[] x) {
        return cdfFromME(new Matrix(alpha), A, x);
    }
}
