/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.des;

import jline.lang.layered.LayeredNetworkStruct;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

/**
 * Result container for LayeredNetwork DES simulation.
 *
 * Metrics are indexed by LQN element type:
 * - Hosts: indices 1 to nhosts
 * - Tasks: indices tshift+1 to tshift+ntasks
 * - Entries: indices eshift+1 to eshift+nentries
 * - Activities: indices ashift+1 to ashift+nacts
 */
public class LNDESResult extends SolverResult {

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

    public LNDESResult() {
        super();
    }
}
