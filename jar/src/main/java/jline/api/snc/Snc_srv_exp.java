/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * MGF service envelope of an exponential server, in JOB units.
 *
 * <p>A single server with Exp(mu) service times completes jobs at the epochs of
 * a Poisson process of rate mu while it is busy, so its cumulative service
 * counted in JOBS is Poisson with mean mu*(t-s) and
 * {@code rho(theta) = mu*(1-exp(-theta))/theta} with a zero burst.</p>
 *
 * <p>THIS IS THE SERVICE ELEMENT TO USE WHENEVER THE WORK UNIT IS THE JOB.
 * Pairing {@link Snc_srv_rate} with a job-counting arrival envelope would model
 * a server that completes jobs at deterministic intervals, an M/D/1, and would
 * UNDERSTATE the delay of an exponential server rather than bound it. The M/M/1
 * read with this element instead reproduces both exact decay rates: the backlog
 * bound decays as (lambda/mu)^n in jobs and the delay bound as
 * exp(-(mu-lambda)*d) in time, since the optimal theta tends to
 * log(mu/lambda).</p>
 *
 * <p>Job units also compose across hops: a departure envelope from
 * {@link Snc_output} is a job count and is directly the arrival envelope of the
 * next station, whereas service-time work units differ from station to
 * station.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_srv_exp.m.</p>
 */
public final class Snc_srv_exp {
    private Snc_srv_exp() {}

    /**
     * @param mu    service rate, jobs per slot
     * @param theta Chernoff parameter, theta &gt; 0
     * @return {0, rho}
     */
    public static double[] snc_srv_exp(double mu, double theta) {
        if (mu <= 0) {
            throw new IllegalArgumentException("snc_srv_exp: mu must be positive, got " + mu);
        }
        if (theta <= 0) {
            throw new IllegalArgumentException("snc_srv_exp: theta must be positive, got " + theta);
        }
        return new double[] {0.0, mu * (1.0 - Math.exp(-theta)) / theta};
    }

    /**
     * @param mu service rate
     * @return the envelope as a function of theta
     */
    public static SncEnvelope of(final double mu) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                return snc_srv_exp(mu, theta);
            }
        };
    }
}
