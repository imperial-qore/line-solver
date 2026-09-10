/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.dpfqn;

/**
 * Normalizing constants of a discrete-time closed cycle of Bernoulli servers
 * with state independent service probabilities.
 *
 * <p>{@link #G} and {@link #G1} are returned in a common scale, which is unity
 * unless the recursion had to be rescaled to keep {@code (q/p)^N} inside double
 * range. Ratios such as {@code G1[k]/G} are exact in either case, and
 * {@link #lG} is always the true log constant.</p>
 */
public final class DpfqnNcResult {

    /** Logarithm of the time-stationary normalizing constant G(N,J). */
    public final double lG;

    /** Time-stationary normalizing constant G(N,J), in the common scale. */
    public final double G;

    /** Arrival normalizing constants, {@code G1[k] = G_1(k,J)} for k = 0..N. */
    public final double[] G1;

    public DpfqnNcResult(double lG, double G, double[] G1) {
        this.lG = lG;
        this.G = G;
        this.G1 = G1;
    }

    /**
     * Throughput per slot, {@code G_1(N,J)/G(N,J)}. It is the same at every
     * node of the cycle, which is the discrete-time counterpart of the
     * continuous-time identity X = G(N-1)/G(N).
     */
    public double throughput() {
        return G1[G1.length - 1] / G;
    }
}
