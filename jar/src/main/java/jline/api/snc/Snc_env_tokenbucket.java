/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * Deterministic token-bucket arrival envelope.
 *
 * <p>A flow policed by a (b,r) token bucket satisfies {@code A(s,t) <= b+r*(t-s)}
 * with probability one, so the envelope is constant in theta. This is the
 * deterministic network calculus arrival curve read as a degenerate MGF
 * envelope, so it mixes freely with the stochastic ones.</p>
 *
 * <p>Port of matlab/src/api/snc/snc_env_tokenbucket.m.</p>
 */
public final class Snc_env_tokenbucket {
    private Snc_env_tokenbucket() {}

    /**
     * @param b bucket depth, units of work
     * @param r token rate, work per slot
     * @return {b, r}
     */
    public static double[] snc_env_tokenbucket(double b, double r) {
        if (b < 0 || r < 0) {
            throw new IllegalArgumentException(
                    "snc_env_tokenbucket: b and r must be nonnegative, got " + b + ", " + r);
        }
        return new double[] {b, r};
    }

    /**
     * @param b     bucket depth
     * @param r     token rate
     * @param theta Chernoff parameter, accepted and ignored
     * @return {b, r}
     */
    public static double[] snc_env_tokenbucket(double b, double r, double theta) {
        if (theta <= 0) {
            throw new IllegalArgumentException("snc_env_tokenbucket: theta must be positive, got " + theta);
        }
        return snc_env_tokenbucket(b, r);
    }

    /**
     * @param b bucket depth
     * @param r token rate
     * @return the envelope as a function of theta
     */
    public static SncEnvelope of(final double b, final double r) {
        return new SncEnvelope() {
            public double[] eval(double theta) {
                return snc_env_tokenbucket(b, r, theta);
            }
        };
    }
}
