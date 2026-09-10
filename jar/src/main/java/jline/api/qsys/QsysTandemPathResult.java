/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

/**
 * Waiting times and epochs of a tandem sample path.
 *
 * <p>Returned by {@link Qsys_tandem_lindley}. Every matrix is indexed
 * {@code [customer][station]}.
 *
 * @since LINE 3.1.0
 */
public class QsysTandemPathResult {
    /** Waiting times, {@code W[n][k]} for customer n at station k. */
    public final double[][] W;
    /**
     * Interarrival times, {@code G[n][k]} between customers n and n+1 at station
     * k, so column 0 is the external stream. The last row is NaN, there being no
     * customer N+1 to separate from.
     */
    public final double[][] G;
    /** Sojourn times, W + S. */
    public final double[][] T;
    /** Departure epochs of each customer from each station. */
    public final double[][] departure;

    public QsysTandemPathResult(double[][] W, double[][] G, double[][] T,
                                double[][] departure) {
        this.W = W;
        this.G = G;
        this.T = T;
        this.departure = departure;
    }
}
