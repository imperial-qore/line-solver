/**
 * @file Per-station metrics behind a MAP flow-equivalent server
 *
 * @since LINE 3.0
 */
package jline.api.fes;

/**
 * Metrics of the stations an aggregate stands for, recovered by conditioning on its
 * population.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class FesMapDeaggregateResult {
    /** Mean queue length per station. */
    public final double[] Q;
    /** Utilization per station. */
    public final double[] U;
    /** Throughput per station. */
    public final double[] X;
    /** Mean residence time per station. */
    public final double[] R;

    public FesMapDeaggregateResult(double[] Q, double[] U, double[] X, double[] R) {
        this.Q = Q;
        this.U = U;
        this.X = X;
        this.R = R;
    }
}
