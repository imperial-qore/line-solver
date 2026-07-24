/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class IntervalPdfFromPH {
    private IntervalPdfFromPH() {}

    /**
     * Returns the approximate probability density function of a
     * continuous phase-type distribution, based on the probability
     * of falling into intervals.
     *
     * The interval pdf is computed as the probability of falling
     * into an interval divided by the interval length.
     *
     * @param alpha Initial probability vector of the phase-type distribution (shape 1 x M)
     * @param A Transient generator matrix of the phase-type distribution (shape M x M)
     * @param intBounds The array of interval boundaries. The pdf is the
     *        probability of falling into an interval divided by
     *        the interval length. If the size of intBounds is K,
     *        the size of the result is K-1.
     * @return A pair of (x, y) where x contains the midpoints of the intervals
     *         and y contains the interval pdf values.
     */
    public static Pair<double[], double[]> intervalPdfFromPH(Matrix alpha, Matrix A, double[] intBounds) {
        int steps = intBounds.length;
        int n = steps - 1;
        double[] x = new double[n];
        double[] y = new double[n];

        for (int i = 0; i < n; i++) {
            x[i] = (intBounds[i + 1] + intBounds[i]) / 2.0;
            double width = intBounds[i + 1] - intBounds[i];
            if (width > 0.0) {
                // CDF(t) = 1 - alpha * expm(A*t) * ones
                double cdfLower;
                if (intBounds[i] <= 0.0) {
                    cdfLower = 0.0;
                } else {
                    cdfLower = 1.0 - alpha.mult(A.scale(intBounds[i]).expm()).elementSum();
                }
                double cdfUpper = 1.0 - alpha.mult(A.scale(intBounds[i + 1]).expm()).elementSum();
                y[i] = (cdfUpper - cdfLower) / width;
            } else {
                y[i] = 0.0;
            }
        }

        return new Pair<double[], double[]>(x, y);
    }

    /**
     * Returns the approximate probability density function of a
     * matrix-exponential distribution, based on the probability
     * of falling into intervals. Equivalent to intervalPdfFromPH for ME distributions.
     */
    public static Pair<double[], double[]> intervalPdfFromME(Matrix alpha, Matrix A, double[] intBounds) {
        return intervalPdfFromPH(alpha, A, intBounds);
    }
}
