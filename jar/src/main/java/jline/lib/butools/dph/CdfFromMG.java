/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import jline.util.matrix.Matrix;

public final class CdfFromMG {
    private CdfFromMG() {}

    /**
     * Returns the cumulative distribution function of a matrix-geometric distribution.
     *
     * @param alpha The initial vector of the matrix-geometric distribution.
     * @param A The matrix parameter of the matrix-geometric distribution.
     * @param x Vector of non-negative integers at which to compute the CDF.
     * @return The probabilities that the matrix-geometrically distributed random variable
     *         is less or equal to the corresponding "x" values.
     */
    public static double[] cdfFromMG(Matrix alpha, Matrix A, int[] x) {
        double[] cdf = new double[x.length];

        for (int i = 0; i < x.length; i++) {
            // CDF(x) = 1 - alpha * A^x * ones
            Matrix APower = Matrix.pow(A, x[i]);
            Matrix term = alpha.mult(APower);
            cdf[i] = 1.0 - term.elementSum();
        }

        return cdf;
    }

    /**
     * Overload for double[] alpha.
     */
    public static double[] cdfFromMG(double[] alpha, Matrix A, int[] x) {
        return cdfFromMG(new Matrix(alpha), A, x);
    }

    /**
     * Overload for single integer x.
     */
    public static double cdfFromMG(Matrix alpha, Matrix A, int x) {
        return cdfFromMG(alpha, A, new int[]{x})[0];
    }

    /**
     * Overload for double[] alpha and single integer x.
     */
    public static double cdfFromMG(double[] alpha, Matrix A, int x) {
        return cdfFromMG(new Matrix(alpha), A, new int[]{x})[0];
    }
}
