/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * Min-plus convolution of two service envelopes (tandem concatenation).
 *
 * <p>Two stations traversed in series offer the flow their min-plus
 * convolution. For independent servers with exponential-form envelopes, summing
 * the geometric series over the intermediate epoch gives</p>
 *
 * <pre>
 *   rho   = min(rho1,rho2),
 *   sigma = sigma1 + sigma2 - log(1-exp(-theta*|rho1-rho2|))/theta.
 * </pre>
 *
 * <p>This is the pay-bursts-only-once result: the end-to-end burst term grows
 * additively rather than the delay bounds of the two stations being summed. The
 * series diverges when the rates are equal, so equal rates are handled by
 * shifting the slower server down by {@code delta}, the usual regularization;
 * delta then trades rate against burst and can be optimized jointly with
 * theta.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_conv.m. Original: F. Ciucu, A. Burchard,
 * J. Liebeherr, "Scaling Properties of Statistical End-to-End Bounds in the
 * Network Calculus", IEEE Trans. Inf. Theory 52(6), 2300-2312, 2006.</p>
 */
public final class Snc_conv {
    private Snc_conv() {}

    /**
     * @param sigma1 burst term of the first station
     * @param rho1   rate term of the first station
     * @param sigma2 burst term of the second station
     * @param rho2   rate term of the second station
     * @param theta  Chernoff parameter, theta &gt; 0
     * @return {sigma, rho} of the concatenated element
     */
    public static double[] snc_conv(double sigma1, double rho1, double sigma2, double rho2,
                                    double theta) {
        return snc_conv(sigma1, rho1, sigma2, rho2, theta, 1e-2 * Math.min(rho1, rho2));
    }

    /**
     * @param sigma1 burst term of the first station
     * @param rho1   rate term of the first station
     * @param sigma2 burst term of the second station
     * @param rho2   rate term of the second station
     * @param theta  Chernoff parameter, theta &gt; 0
     * @param delta  rate separation used when the two rates coincide
     * @return {sigma, rho} of the concatenated element
     */
    public static double[] snc_conv(double sigma1, double rho1, double sigma2, double rho2,
                                    double theta, double delta) {
        if (theta <= 0) {
            throw new IllegalArgumentException("snc_conv: theta must be positive, got " + theta);
        }
        if (delta <= 0) {
            throw new IllegalArgumentException("snc_conv: delta must be positive, got " + delta);
        }
        if (!Double.isFinite(rho1) || !Double.isFinite(rho2)) {
            return new double[] {Double.POSITIVE_INFINITY, Math.min(rho1, rho2)};
        }
        double gap = Math.abs(rho1 - rho2);
        double rho;
        if (gap <= delta) {
            gap = delta; // equal rates: shift the slower server down to close the series
            rho = Math.min(rho1, rho2) - delta;
        } else {
            rho = Math.min(rho1, rho2);
        }
        if (rho <= 0) {
            return new double[] {Double.POSITIVE_INFINITY, rho};
        }
        return new double[] {sigma1 + sigma2 - Math.log(1.0 - Math.exp(-theta * gap)) / theta, rho};
    }

    /**
     * @param s1 the first service element
     * @param s2 the second service element
     * @return the concatenated element as a function of theta
     */
    public static SncEnvelope of(final SncEnvelope s1, final SncEnvelope s2) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                double[] a = s1.eval(theta);
                double[] b = s2.eval(theta);
                return snc_conv(a[0], a[1], b[0], b[1], theta);
            }
        };
    }
}
