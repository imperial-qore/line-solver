/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * Sojourn-time moments of the multiclass M/M/1-PS queue.
 *
 * <p>Port of MATLAB qsys_mm1_ps.m.</p>
 */
public class QsysMm1PsResult {
    /** Per-class mean sojourn times. */
    public final double[] W;
    /** Per-class second moments of the sojourn time. */
    public final double[] W2;
    /** Unutilized fraction of the processor, 1 - sum_j lambda_j/mu_j. */
    public final double alpha;

    public QsysMm1PsResult(double[] W, double[] W2, double alpha) {
        this.W = W;
        this.W2 = W2;
        this.alpha = alpha;
    }
}
