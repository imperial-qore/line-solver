/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.mam;

/**
 * Options for LDQBD solver.
 */
public final class LdqbdOptions {
    private final double epsilon;
    private final int maxIter;
    private final boolean verbose;

    public LdqbdOptions() {
        this(1e-10, 1000, false);
    }

    public LdqbdOptions(double epsilon, int maxIter, boolean verbose) {
        this.epsilon = epsilon;
        this.maxIter = maxIter;
        this.verbose = verbose;
    }

    public double getEpsilon() { return epsilon; }
    public int getMaxIter() { return maxIter; }
    public boolean getVerbose() { return verbose; }
}
