/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.analyzers;

/**
 * Transient refined-mean-field cache trajectory for a network, produced by
 * {@link FluidCacheTran}. Java analog of the outputs of the MATLAB
 * {@code solver_fld_cacheqn_tran.m}.
 */
public final class CacheTranResult {
    /** Shared time grid (nt,). */
    public final double[] t;
    /** Cache node indices (rows ordered as found by nodetype==Cache). */
    public final int[] caches;
    /** Per-cache per-class hit-probability trajectory (ncaches x K x nt). */
    public final double[][][] hitprob;
    /** Per-cache per-class miss-probability trajectory (ncaches x K x nt). */
    public final double[][][] missprob;
    /** Per-cache per-class arrival rate (ncaches x K). */
    public final double[][] arate;
    /** Per-cache full DDPP occupancy trajectory (ncaches x dim x nt). */
    public final double[][][] xocc;

    public CacheTranResult(double[] t, int[] caches, double[][][] hitprob,
                           double[][][] missprob, double[][] arate, double[][][] xocc) {
        this.t = t;
        this.caches = caches;
        this.hitprob = hitprob;
        this.missprob = missprob;
        this.arate = arate;
        this.xocc = xocc;
    }
}
