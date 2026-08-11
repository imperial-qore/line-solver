/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.util.matrix.Matrix;

public final class PdfFromME {
    private PdfFromME() {}

    /**
     * Returns the probability density function of a matrix-exponential distribution.
     *
     * @param alpha The initial vector of the matrix-exponential distribution.
     * @param A The matrix parameter of the matrix-exponential distribution.
     * @param x The points at which the density function will be computed.
     * @return The values of the density function at the corresponding "x" values.
     */
    public static double[] pdfFromME(Matrix alpha, Matrix A, double[] x) {
        double[] pdf = new double[x.length];
        Matrix negA = A.neg();
        for (int i = 0; i < x.length; i++) {
            if (x[i] < 0) {
                pdf[i] = 0.0;
            } else {
                Matrix expAt = A.scale(x[i]).expm();
                pdf[i] = alpha.mult(expAt).mult(negA).elementSum();
            }
        }
        return pdf;
    }

    /**
     * Returns the probability density function of a phase-type distribution.
     */
    public static double[] pdfFromPH(Matrix alpha, Matrix A, double[] x) {
        return pdfFromME(alpha, A, x);
    }

    /**
     * Overload for double[] alpha.
     */
    public static double[] pdfFromME(double[] alpha, Matrix A, double[] x) {
        return pdfFromME(new Matrix(alpha), A, x);
    }

    /**
     * Overload for double[] alpha.
     */
    public static double[] pdfFromPH(double[] alpha, Matrix A, double[] x) {
        return pdfFromME(new Matrix(alpha), A, x);
    }
}
