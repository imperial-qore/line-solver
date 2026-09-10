/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * MGF service envelope of a constant-rate work-conserving server.
 *
 * <p>{@code S(s,t) = C*(t-s)}, so {@code (sigma, rho) = (0, C)} for every theta.
 * This is the elementary service element; a station shared by cross traffic is
 * obtained from it through {@link Snc_leftover} and a tandem through
 * {@link Snc_conv}.</p>
 *
 * <p>USE {@link Snc_srv_exp} INSTEAD WHENEVER THE WORK UNIT IS THE JOB: pairing
 * this element with a job-counting arrival envelope models an M/D/1 and
 * understates the delay of an exponential server.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_srv_rate.m.</p>
 */
public final class Snc_srv_rate {
    private Snc_srv_rate() {}

    /**
     * @param C server capacity, work per slot
     * @return {0, C}
     */
    public static double[] snc_srv_rate(double C) {
        if (C <= 0) {
            throw new IllegalArgumentException("snc_srv_rate: C must be positive, got " + C);
        }
        return new double[] {0.0, C};
    }

    /**
     * @param C     server capacity
     * @param theta Chernoff parameter, accepted and ignored
     * @return {0, C}
     */
    public static double[] snc_srv_rate(double C, double theta) {
        if (theta <= 0) {
            throw new IllegalArgumentException("snc_srv_rate: theta must be positive, got " + theta);
        }
        return snc_srv_rate(C);
    }

    /**
     * @param C server capacity
     * @return the envelope as a function of theta
     */
    public static SncEnvelope of(final double C) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                return snc_srv_rate(C, theta);
            }
        };
    }
}
