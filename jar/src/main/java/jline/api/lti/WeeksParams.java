/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.lti;

/**
 * A Laguerre expansion: the exponential damping, the time scaling and the
 * coefficients they were computed at.
 *
 * <p>The three travel together on purpose. Evaluating a coefficient set at any
 * other (sigma, b) silently returns a different function, so
 * {@link Laplace_invert#laplace_invert_weeks} takes this object rather than
 * three loose numbers.
 */
public class WeeksParams {
    /** Exponential damping applied before the expansion. */
    public final double sigma;
    /** Time scaling applied before the expansion. */
    public final double b;
    /** Laguerre coefficients q_n, n = 0..2*p0-1. */
    public final double[] q;

    public WeeksParams(double sigma, double b, double[] q) {
        this.sigma = sigma;
        this.b = b;
        this.q = q;
    }
}
