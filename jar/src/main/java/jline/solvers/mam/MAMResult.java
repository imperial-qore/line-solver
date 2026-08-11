/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import java.util.List;

public class MAMResult extends SolverResult {
    public int iter;
    public double lG;
    public Matrix actionRates;
    public List<Matrix> equilibrium;
    public List<Matrix> generators;
    /**
     * Orbit-level internals of the matrix-analytic retrial engine, when the model
     * was solved as a retrial station. Exposed through SolverMAM.getMAMResult so
     * that the orbit-level stationary distribution, the truncation level and its
     * residual can be inspected rather than only the mean measures.
     */
    public jline.api.qsys.QsysRetrialResult retrialInternals;
}
