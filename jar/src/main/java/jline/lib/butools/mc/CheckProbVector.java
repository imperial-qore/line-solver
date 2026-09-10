/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.mc;

import jline.util.matrix.Matrix;

public final class CheckProbVector {
    private CheckProbVector() {}

    /**
     * Checks if the vector is a valid probability vector: the
     * vector has only non-negative elements, the sum of the
     * vector elements is 1.
     *
     * @param pi The vector to check.
     * @param sub If false, checks for stochastic; if true, checks for sub-stochastic.
     * @param prec Numerical precision.
     * @return The result of the check.
     */
    public static boolean checkProbVector(Matrix pi, boolean sub, double prec) {
        int n = pi.length();

        // Check for negative elements
        for (int i = 0; i < n; i++) {
            if (pi.get(i) < -prec) {
                return false;
            }
        }

        double sum = pi.elementSum();

        if (sub) {
            if (sum > 1 + prec * n) {
                return false;
            }
        } else {
            if (Math.abs(sum - 1.0) > prec * n) {
                return false;
            }
        }

        return true;
    }

    public static boolean checkProbVector(Matrix pi, boolean sub) {
        return checkProbVector(pi, sub, 1e-14);
    }

    public static boolean checkProbVector(Matrix pi) {
        return checkProbVector(pi, false, 1e-14);
    }

    public static boolean checkProbVector(double[] pi, boolean sub, double prec) {
        return checkProbVector(new Matrix(pi), sub, prec);
    }

    public static boolean checkProbVector(double[] pi, boolean sub) {
        return checkProbVector(new Matrix(pi), sub, 1e-14);
    }

    public static boolean checkProbVector(double[] pi) {
        return checkProbVector(new Matrix(pi), false, 1e-14);
    }
}
