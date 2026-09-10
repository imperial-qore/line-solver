/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ln;

import jline.util.matrix.Matrix;

/**
 * Transient average metrics of a layered network, assembled block-diagonally
 * across the ensemble layers by {@link SolverLN#getTranAvg()}.
 *
 * <p>Each layer e occupies a disjoint block of rows (its stations) and columns
 * (its classes), mirroring the MATLAB SolverLN.getTranAvg cell-array layout.
 * Off-block cells and disabled metrics are left null.</p>
 *
 * <p>Each populated cell {@code QNt[i][r]} (and {@code UNt}, {@code TNt}) is a
 * column vector of length Tmax holding the metric time series; {@code t[i][r]}
 * is the matching column vector of time instants for that cell.</p>
 */
public class LNTranAvgResult {

    /** Transient queue lengths, block-diagonal [rows][cols]; each cell Tmax x 1 or null. */
    public Matrix[][] QNt;

    /** Transient utilizations, block-diagonal [rows][cols]; each cell Tmax x 1 or null. */
    public Matrix[][] UNt;

    /** Transient throughputs, block-diagonal [rows][cols]; each cell Tmax x 1 or null. */
    public Matrix[][] TNt;

    /** Time instants per cell, block-diagonal [rows][cols]; each cell Tmax x 1 or null. */
    public Matrix[][] t;

    public LNTranAvgResult(Matrix[][] QNt, Matrix[][] UNt, Matrix[][] TNt, Matrix[][] t) {
        this.QNt = QNt;
        this.UNt = UNt;
        this.TNt = TNt;
        this.t = t;
    }
}
