/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.analyzers;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

public interface FluidAnalyzer {

    void analyze(NetworkStruct sn, SolverOptions options, SolverResult result);

    Matrix getXVecIt();
}
