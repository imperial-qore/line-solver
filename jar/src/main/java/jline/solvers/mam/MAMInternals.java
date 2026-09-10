/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.api.qsys.QsysBmapM1Result;
import jline.api.qsys.QsysRetrialResult;

/**
 * Intermediate quantities of the matrix-analytic analysis of a single-queue
 * model, in addition to the mean performance measures returned by getAvg.
 *
 * <p>Mean values alone hide the objects the method is actually built on, so a
 * matrix-analytic result cannot be inspected, taught, or checked against a
 * published derivation. This carries them.</p>
 *
 * <p>Exactly one of the two fields is populated: for a BMAP (or MAP) arrival
 * stream feeding an exponential single server the result is that of
 * {@code qsys_bmapm1} and carries the M/G/1-type quantities; for a retrial
 * station it is that of {@code qsys_bmapphnn_retrial} and carries the
 * orbit-level stationary distribution together with the truncation level and its
 * residual.</p>
 */
public class MAMInternals {
    /** M/G/1-type internals of a BMAP/M/1 queue, or null for a retrial station. */
    public final QsysBmapM1Result bmapm1;
    /** Orbit-level internals of a retrial station, or null otherwise. */
    public final QsysRetrialResult retrial;

    public MAMInternals(QsysBmapM1Result bmapm1, QsysRetrialResult retrial) {
        this.bmapm1 = bmapm1;
        this.retrial = retrial;
    }
}
