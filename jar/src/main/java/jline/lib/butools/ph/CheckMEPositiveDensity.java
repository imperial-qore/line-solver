/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.util.matrix.Matrix;

public final class CheckMEPositiveDensity {
    private CheckMEPositiveDensity() {}

    /**
     * Checks if the given ME distribution has positive density.
     *
     * This is done by trying to transform the ME distribution to
     * a monocyclic PH representation. If successful, the density
     * is positive.
     *
     * @param alpha Initial vector of the matrix-exponential distribution to check (shape 1 x M)
     * @param A Matrix parameter of the matrix-exponential distribution to check (shape M x M)
     * @param maxSize The procedure tries to transform the ME distribution
     *        to phase-type up to order maxSize. The default value is 100.
     * @param prec Numerical precision. The default value is 1e-14.
     * @return True, if the given matrix-exponential distribution has a positive density.
     *
     * Notes:
     *   This procedure calls MonocyclicPHFromME, and can be time consuming.
     */
    public static boolean checkMEPositiveDensity(Matrix alpha, Matrix A, int maxSize, double prec) {
        if (!CheckMERepresentation.checkMERepresentation(alpha, A, prec)) {
            return false;
        }
        try {
            PHRepresentation result = MonocyclicPHFromME.monocyclicPHFromME(alpha, A, maxSize, prec);
            return CheckPHRepresentation.checkPHRepresentation(result.getAlpha(), result.getA(), prec);
        } catch (Exception e) {
            return false;
        }
    }

    public static boolean checkMEPositiveDensity(Matrix alpha, Matrix A, int maxSize) {
        return checkMEPositiveDensity(alpha, A, maxSize, 1e-14);
    }

    public static boolean checkMEPositiveDensity(Matrix alpha, Matrix A) {
        return checkMEPositiveDensity(alpha, A, 100, 1e-14);
    }
}
