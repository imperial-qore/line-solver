/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import jline.lang.layered.LayeredNetworkStruct;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

/**
 * Result container for LayeredNetwork LDES simulation.
 *
 * Metrics are indexed by LQN element type:
 * - Hosts: indices 1 to nhosts
 * - Tasks: indices tshift+1 to tshift+ntasks
 * - Entries: indices eshift+1 to eshift+nentries
 * - Activities: indices ashift+1 to ashift+nacts
 */
public class LNLDESResult extends SolverResult {

    /** LQN structure reference */
    public LayeredNetworkStruct lsn;

    /** Queue length per LQN element [1 x nidx] */
    public Matrix QLN;

    /** Utilization per LQN element [1 x nidx] */
    public Matrix ULN;

    /** Response time per LQN element [1 x nidx] */
    public Matrix RLN;

    /** Residence/waiting time per LQN element [1 x nidx] */
    public Matrix WLN;

    /** Throughput per LQN element [1 x nidx] */
    public Matrix TLN;

    /** Arrival rate per LQN element [1 x nidx] */
    public Matrix ALN;

    /** Queue length confidence interval half-widths [1 x nidx] */
    public Matrix QLNCI;

    /** Utilization confidence interval half-widths [1 x nidx] */
    public Matrix ULNCI;

    /** Response time confidence interval half-widths [1 x nidx] */
    public Matrix RLNCI;

    /** Throughput confidence interval half-widths [1 x nidx] */
    public Matrix TLNCI;

    // --- Per-CacheTask metrics for LCQ models (getAvgCacheTable). Parallel lists,
    // one entry per (CacheTask, ItemEntry) read. Empty when the LQN has no cache. ---
    /** CacheTask absolute LQN index for each cache read. */
    public java.util.List<Integer> cacheTaskIdx = new java.util.ArrayList<Integer>();
    /** ItemEntry absolute LQN index for each cache read. */
    public java.util.List<Integer> cacheItemEntryIdx = new java.util.ArrayList<Integer>();
    /** True-hit probability (excludes delayed hits). */
    public java.util.List<Double> cacheHitProb = new java.util.ArrayList<Double>();
    /** Miss probability. */
    public java.util.List<Double> cacheMissProb = new java.util.ArrayList<Double>();
    /** Delayed-hit probability (0 when no retrieval system). */
    public java.util.List<Double> cacheDelayedProb = new java.util.ArrayList<Double>();
    /** Read throughput (reads per unit time) into this cache. */
    public java.util.List<Double> cacheReadRate = new java.util.ArrayList<Double>();

    public LNLDESResult() {
        super();
    }
}
