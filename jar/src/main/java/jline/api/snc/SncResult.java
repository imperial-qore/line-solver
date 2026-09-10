/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.snc;

/**
 * The value of a Chernoff bound together with the theta that attains it.
 *
 * <p>The minimizing theta is reported rather than hidden because it is
 * diagnostic: on an M/M/1 read in job units it tends to {@code log(mu/lambda)},
 * which is exactly what makes the backlog decay rate exact, and a theta pinned
 * at the edge of the search range says the range was too narrow.</p>
 */
public final class SncResult {
    /** The bound: a violation probability, a quantile or a mean bound. */
    public final double value;
    /** The minimizing theta, NaN when no feasible theta exists. */
    public final double theta;

    public SncResult(double value, double theta) {
        this.value = value;
        this.theta = theta;
    }
}
