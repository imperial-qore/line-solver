/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * Output (departure) arrival envelope of a flow leaving a server.
 *
 * <p>For independent processes and a stable station,</p>
 *
 * <pre>
 *   rho   = rhoA,
 *   sigma = sigmaA + sigmaS - log(1-exp(-theta*(rhoS-rhoA)))/theta.
 * </pre>
 *
 * <p>The rate is conserved and the server adds burstiness. This is what carries
 * a flow across a feed-forward network one hop at a time; for a tandem
 * traversed by the same flow, {@link Snc_conv} gives the tighter end-to-end
 * result.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_output.m.</p>
 */
public final class Snc_output {
    private Snc_output() {}

    /**
     * @param sigmaA burst term of the arrival envelope
     * @param rhoA   rate term of the arrival envelope
     * @param sigmaS burst term of the service envelope
     * @param rhoS   rate term of the service envelope
     * @param theta  Chernoff parameter, theta &gt; 0
     * @return {sigma, rho} of the departure envelope; sigma is Infinity when
     *         the station is unstable
     */
    public static double[] snc_output(double sigmaA, double rhoA, double sigmaS, double rhoS,
                                      double theta) {
        if (theta <= 0) {
            throw new IllegalArgumentException("snc_output: theta must be positive, got " + theta);
        }
        if (!Double.isFinite(rhoA) || !Double.isFinite(rhoS) || rhoS <= rhoA) {
            return new double[] {Double.POSITIVE_INFINITY, rhoA};
        }
        return new double[] {
                sigmaA + sigmaS - Math.log(1.0 - Math.exp(-theta * (rhoS - rhoA))) / theta, rhoA};
    }

    /**
     * @param arv the arrival envelope entering the server
     * @param srv the service envelope
     * @return the departure envelope as a function of theta
     */
    public static SncEnvelope of(final SncEnvelope arv, final SncEnvelope srv) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                double[] a = arv.eval(theta);
                double[] s = srv.eval(theta);
                return snc_output(a[0], a[1], s[0], s[1], theta);
            }
        };
    }
}
