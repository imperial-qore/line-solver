/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.util.matrix.Matrix;

/**
 * One cache node's average results, as a portable RECORD rather than a write on
 * the node object.
 * <p>
 * A Cache does not report in Q/U/R/T: it carries its result on the node, which
 * is why a solver that hands its work to another model (SolverENV under
 * mapEnvApprox, say) writes the split onto a model the caller never sees. This
 * record is the hand-back: keyed by node NAME, never by position, so it survives
 * a model whose node order differs.
 * <p>
 * EVERY FIELD IS OPTIONAL AND ABSENT (null) MEANS NOT COMPUTED, never zero. A
 * recombination that cannot reach a quantity must leave it null so the caller
 * reports it missing rather than reporting a fabricated zero.
 *
 * @see jline.solvers.NetworkSolver#applyCacheResults
 */
public class CacheResultRecord {

    /** Name of the cache node these results belong to. */
    public String name;

    /** Node index in the model that produced the record, for diagnostics only. */
    public int node;

    /** Actual hit probability per read class [1 x classes], or null. */
    public Matrix hitprob;

    /** Actual miss probability per read class [1 x classes], or null. */
    public Matrix missprob;

    /** Actual delayed-hit probability per read class [1 x classes], or null. */
    public Matrix delayedhitprob;

    /** Actual hit probability by cache list [classes x lists], or null. */
    public Matrix hitproblist;

    public CacheResultRecord() {
    }

    public CacheResultRecord(String name, int node) {
        this.name = name;
        this.node = node;
    }
}
