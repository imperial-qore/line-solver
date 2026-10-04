/**
 * @file Result of the reduced-model solve behind a MAP flow-equivalent server
 *
 * @since LINE 3.0
 */
package jline.api.fes;

/**
 * System metrics of the closed model made of a delay and a MAP flow-equivalent server.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class FesMapSolveResult {
    /** System throughput. */
    public final double X;
    /** Mean response time of the aggregated subnetwork, N/X - E[Z]. */
    public final double R;
    /** Mean number of jobs held by the flow-equivalent server. */
    public final double Q;
    /** Distribution of the jobs held by the flow-equivalent server, index k = P(k jobs). */
    public final double[] pk;

    public FesMapSolveResult(double X, double R, double Q, double[] pk) {
        this.X = X;
        this.R = R;
        this.Q = Q;
        this.pk = pk;
    }
}
