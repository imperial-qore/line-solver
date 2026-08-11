/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.dph;

import jline.util.matrix.Matrix;

public final class PmfFromDPH {
    private PmfFromDPH() {}

    /**
     * Returns the probability mass function of a discrete phase-type distribution.
     *
     * @param alpha The initial probability vector of the discrete phase-type distribution.
     *        The sum of the entries of alpha is less or equal to 1.
     * @param A The transient generator matrix of the discrete phase-type distribution.
     * @param x Vector of non-negative integers at which to compute the PMF.
     * @return The probabilities that the discrete phase type distributed random variable
     *         takes the corresponding "x" values.
     */
    public static double[] pmfFromDPH(Matrix alpha, Matrix A, int[] x) {
        return PmfFromMG.pmfFromMG(alpha, A, x);
    }

    /**
     * Overload for double[] alpha.
     */
    public static double[] pmfFromDPH(double[] alpha, Matrix A, int[] x) {
        return pmfFromDPH(new Matrix(alpha), A, x);
    }

    /**
     * Overload for single integer x.
     */
    public static double pmfFromDPH(Matrix alpha, Matrix A, int x) {
        return pmfFromDPH(alpha, A, new int[] { x })[0];
    }

    /**
     * Overload for double[] alpha and single integer x.
     */
    public static double pmfFromDPH(double[] alpha, Matrix A, int x) {
        return pmfFromDPH(new Matrix(alpha), A, new int[] { x })[0];
    }
}
