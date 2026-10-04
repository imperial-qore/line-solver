/**
 * @file Result of a MAP flow-equivalent server aggregation
 *
 * @since LINE 3.0
 */
package jline.api.fes;

import java.util.List;

import jline.util.matrix.MatrixCell;

/**
 * The load-dependent MAP that replaces an aggregated subnetwork, with the descriptors it
 * was fitted from.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class FesMapAggregateResult {
    /** Service process of the flow-equivalent server, index k-1 holding k jobs. */
    public final List<MatrixCell> fes;
    /** Throughput of the subnetwork at each population. */
    public final double[] throughput;
    /** Descriptors e1, e2, e3 and the index of dispersion at each population. */
    public final double[][] moments;
    /** Fit status at each population, see Map2_fit_idc. */
    public final int[] status;
    /** Populations at which the inter-departure MAP was evaluated. */
    public final int[] grid;

    public FesMapAggregateResult(List<MatrixCell> fes, double[] throughput, double[][] moments,
                                 int[] status, int[] grid) {
        this.fes = fes;
        this.throughput = throughput;
        this.moments = moments;
        this.status = status;
        this.grid = grid;
    }
}
