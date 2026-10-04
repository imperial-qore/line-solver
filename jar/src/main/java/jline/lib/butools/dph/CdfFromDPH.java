/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import jline.util.matrix.Matrix;

public final class CdfFromDPH {
    private CdfFromDPH() {}

    /**
     * Returns the cumulative distribution function of a discrete phase-type distribution.
     *
     * @param alpha The initial probability vector of the discrete phase-type distribution.
     * @param A The transition probability matrix of the discrete phase-type distribution.
     * @param x Vector of non-negative integers at which to compute the CDF.
     * @return The probabilities that the discrete phase type distributed random variable
     *         is less or equal to the corresponding "x" values.
     */
    public static double[] cdfFromDPH(Matrix alpha, Matrix A, int[] x) {
        return CdfFromMG.cdfFromMG(alpha, A, x);
    }

    /**
     * Overload for double[] alpha.
     */
    public static double[] cdfFromDPH(double[] alpha, Matrix A, int[] x) {
        return cdfFromDPH(new Matrix(alpha), A, x);
    }

    /**
     * Overload for single integer x.
     */
    public static double cdfFromDPH(Matrix alpha, Matrix A, int x) {
        return cdfFromDPH(alpha, A, new int[] { x })[0];
    }

    /**
     * Overload for double[] alpha and single integer x.
     */
    public static double cdfFromDPH(double[] alpha, Matrix A, int x) {
        return cdfFromDPH(new Matrix(alpha), A, new int[] { x })[0];
    }
}
