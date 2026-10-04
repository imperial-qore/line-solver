/**
 * @file Shape-preserving interpolation of the flow-equivalent descriptors
 *
 * @since LINE 3.0
 */
package jline.api.fes;

/**
 * Monotone piecewise cubic Hermite interpolation.
 *
 * Interpolates the descriptors of a MAP flow-equivalent server between the populations at
 * which they were evaluated. Fritsch and Carlson slopes are used, with the noncentered
 * three-point endpoint rule of de Boor, so the interpolant never overshoots and a monotone
 * sequence of throughputs stays monotone. The algorithm is written out rather than
 * delegated to a library so that the MATLAB, Java, Python and C++ ports agree.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class Fes_map_interp {
    private Fes_map_interp() {}

    /**
     * Interpolates one series at the query abscissae.
     *
     * @param x  sample abscissae, strictly increasing
     * @param y  sample values
     * @param xq query abscissae
     * @return interpolated values
     */
    public static double[] fes_map_interp(double[] x, double[] y, double[] xq) {
        int n = x.length;
        double[] yq = new double[xq.length];
        if (n == 1) {
            for (int q = 0; q < xq.length; q++) {
                yq[q] = y[0];
            }
            return yq;
        }

        double[] h = new double[n - 1];
        double[] delta = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            h[i] = x[i + 1] - x[i];
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }

        double[] d = new double[n];
        if (n == 2) {
            d[0] = delta[0];
            d[1] = delta[0];
        } else {
            for (int i = 1; i < n - 1; i++) {
                if (delta[i - 1] * delta[i] > 0) {
                    double w1 = 2 * h[i] + h[i - 1];
                    double w2 = h[i] + 2 * h[i - 1];
                    d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
                }
            }
            d[0] = edgeSlope(h[0], h[1], delta[0], delta[1]);
            d[n - 1] = edgeSlope(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);
        }

        for (int q = 0; q < xq.length; q++) {
            double t = xq[q];
            int i;
            if (t <= x[0]) {
                i = 0;
            } else if (t >= x[n - 1]) {
                i = n - 2;
            } else {
                i = 0;
                while (i < n - 2 && x[i + 1] <= t) {
                    i++;
                }
            }
            double s = t - x[i];
            double c2 = (3 * delta[i] - 2 * d[i] - d[i + 1]) / h[i];
            double c3 = (d[i] - 2 * delta[i] + d[i + 1]) / (h[i] * h[i]);
            yq[q] = y[i] + s * (d[i] + s * (c2 + s * c3));
        }
        return yq;
    }

    /** Noncentered three-point endpoint slope with the monotonicity clamps of de Boor. */
    private static double edgeSlope(double h1, double h2, double del1, double del2) {
        double d = ((2 * h1 + h2) * del1 - h1 * del2) / (h1 + h2);
        if (Math.signum(d) != Math.signum(del1)) {
            return 0;
        }
        if (Math.signum(del1) != Math.signum(del2) && Math.abs(d) > Math.abs(3 * del1)) {
            return 3 * del1;
        }
        return d;
    }
}
