/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * Leftover service envelope under blind (arbitrary) multiplexing.
 *
 * <p>A server with envelope (sigmaS,rhoS) shared with a cross flow of arrival
 * envelope (sigmaX,rhoX) leaves the flow of interest the service process
 * {@code S-X}, whose envelope is {@code rho = rhoS-rhoX},
 * {@code sigma = sigmaS+sigmaX}.</p>
 *
 * <p>The subtraction of the two exponential forms is exact when the processes
 * are INDEPENDENT; otherwise the pair must be split by Hoelder's inequality,
 * which this elementary version does not do. A nonpositive rho means the cross
 * traffic can exhaust the server, which the bound functions report as a
 * violation probability of 1.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_leftover.m.</p>
 */
public final class Snc_leftover {
    private Snc_leftover() {}

    /**
     * @param sigmaS burst term of the service envelope
     * @param rhoS   rate term of the service envelope
     * @param sigmaX burst term of the cross-flow arrival envelope
     * @param rhoX   rate term of the cross-flow arrival envelope
     * @return {sigma, rho} of the leftover service
     */
    public static double[] snc_leftover(double sigmaS, double rhoS, double sigmaX, double rhoX) {
        return new double[] {sigmaS + sigmaX, rhoS - rhoX};
    }

    /**
     * @param srv   the service envelope
     * @param cross the cross-flow arrival envelope
     * @return the leftover service as a function of theta
     */
    public static SncEnvelope of(final SncEnvelope srv, final SncEnvelope cross) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                double[] s = srv.eval(theta);
                double[] x = cross.eval(theta);
                return snc_leftover(s[0], s[1], x[0], x[1]);
            }
        };
    }
}
