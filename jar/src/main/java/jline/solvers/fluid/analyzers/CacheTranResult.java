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
    /**
     * Per-cache isolated per-class per-item request rates (ncaches x classes),
     * each entry the list-0 column of that class's lambda_cache. They are the
     * weights that turn {@link #xocc} into a hit probability RESOLVED BY LIST,
     * which the total hit probability alone cannot supply.
     */
    public final double[][][] lambdaItem;
    /**
     * Per-cache grid-independent average of the transient against a stage
     * holding time ({@link jline.api.cache.CacheSojourn}), or null when no
     * holding time was given; an entry is null for a cache with no request rate.
     */
    public jline.api.cache.CacheSojourn.Result[] sojourn;
    /** Per-cache per-class hit probability of the sojourn-averaged occupancy (ncaches x K), or null. */
    public double[][] sojournHit;
    /** Per-cache per-class miss probability of the sojourn-averaged occupancy (ncaches x K), or null. */
    public double[][] sojournMiss;

    public CacheTranResult(double[] t, int[] caches, double[][][] hitprob,
                           double[][][] missprob, double[][] arate, double[][][] xocc,
                           double[][][] lambdaItem) {
        this.t = t;
        this.caches = caches;
        this.hitprob = hitprob;
        this.missprob = missprob;
        this.arate = arate;
        this.xocc = xocc;
        this.lambdaItem = lambdaItem;
    }

    public CacheTranResult(double[] t, int[] caches, double[][][] hitprob,
                           double[][][] missprob, double[][] arate, double[][][] xocc) {
        this(t, caches, hitprob, missprob, arate, xocc, null);
    }
}
