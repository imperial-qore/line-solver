/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * MGF arrival envelope of a Poisson flow with unit-size jobs.
 *
 * <p>{@code log E[exp(theta*A(0,t))] = lambda*t*(exp(theta)-1)} exactly, so the
 * envelope is tight with a zero burst term:
 * {@code rho(theta) = lambda*(exp(theta)-1)/theta}.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_env_poisson.m.</p>
 */
public final class Snc_env_poisson {
    private Snc_env_poisson() {}

    /**
     * @param lambda arrival rate, jobs per slot
     * @param theta  Chernoff parameter, theta &gt; 0
     * @return {sigma, rho}, with sigma = 0
     */
    public static double[] snc_env_poisson(double lambda, double theta) {
        if (lambda < 0) {
            throw new IllegalArgumentException("snc_env_poisson: lambda must be nonnegative, got " + lambda);
        }
        if (theta <= 0) {
            throw new IllegalArgumentException("snc_env_poisson: theta must be positive, got " + theta);
        }
        return new double[] {0.0, lambda * (Math.exp(theta) - 1.0) / theta};
    }

    /**
     * @param lambda arrival rate, jobs per slot
     * @return the envelope as a function of theta
     */
    public static SncEnvelope of(final double lambda) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                return snc_env_poisson(lambda, theta);
            }
        };
    }
}
