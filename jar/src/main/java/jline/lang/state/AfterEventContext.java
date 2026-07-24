/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.state;

import java.io.Serializable;
import java.util.Map;

import jline.lang.nodes.Station;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

/**
 * Loop-invariant context precomputed by {@link State#afterEventInit} for
 * {@link State#afterEvent}. Hot callers (the SSA Gillespie loop re-evaluates
 * every synchronization at every step) pass it to skip the per-call
 * derivation of ismkvmodclass, lldscaling and cdscaling; semantics are
 * identical to the per-call path.
 *
 * Build the context from the SAME NetworkStruct instance later passed to
 * afterEvent, after any caller-side rewrite of its fields (see the
 * Solver_ssa preamble that rewrites nservers/cap/classcap in place).
 */
public class AfterEventContext implements Serializable {

    /** Load-dependent scaling matrix (defaulted to ones if sn.lldscaling is empty). */
    public final Matrix lldscaling;

    /** Number of columns of lldscaling. */
    public final int lldlimit;

    /** Class-dependent scaling functions (defaulted to the constant 1 map). */
    public final Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling;

    /** Per-node (R x 1) indicator of MAP/MMPP2/BMAP service; stations only. */
    public final Map<Integer, Matrix> ismkvmodclass;

    public AfterEventContext(Matrix lldscaling, int lldlimit,
                             Map<Station, SerializableFunction<Matrix, Matrix>> cdscaling,
                             Map<Integer, Matrix> ismkvmodclass) {
        this.lldscaling = lldscaling;
        this.lldlimit = lldlimit;
        this.cdscaling = cdscaling;
        this.ismkvmodclass = ismkvmodclass;
    }
}
