package jline.api.fes;

import java.util.Objects;

import jline.util.matrix.Matrix;

/**
 * Options for Flow-Equivalent Server (FES) aggregation.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class FESOptions {
    private final String solver;
    private final Matrix cutoffs;
    private final boolean verbose;

    public FESOptions(String solver, Matrix cutoffs, boolean verbose) {
        this.solver = solver;
        this.cutoffs = cutoffs;
        this.verbose = verbose;
    }

    public FESOptions() {
        this("mva", null, false);
    }

    /** Solver to use for throughput computation ('mva' default). */
    public String getSolver() {
        return solver;
    }

    /** Per-class population cutoffs (default: null, uses total jobs per class). */
    public Matrix getCutoffs() {
        return cutoffs;
    }

    /** Whether verbose output is enabled. */
    public boolean isVerbose() {
        return verbose;
    }

    /**
     * Create default options.
     */
    public static FESOptions defaults() {
        return new FESOptions();
    }

    /**
     * Create options with specified solver.
     */
    public static FESOptions withSolver(String solver) {
        return new FESOptions(solver, null, false);
    }

    /**
     * Create options with specified cutoffs.
     */
    public static FESOptions withCutoffs(Matrix cutoffs) {
        return new FESOptions("mva", cutoffs, false);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof FESOptions)) return false;
        FESOptions that = (FESOptions) o;
        return verbose == that.verbose
                && Objects.equals(solver, that.solver)
                && Objects.equals(cutoffs, that.cutoffs);
    }

    @Override
    public int hashCode() {
        return Objects.hash(solver, cutoffs, verbose);
    }

    @Override
    public String toString() {
        return "FESOptions(solver=" + solver + ", cutoffs=" + cutoffs + ", verbose=" + verbose + ")";
    }
}
