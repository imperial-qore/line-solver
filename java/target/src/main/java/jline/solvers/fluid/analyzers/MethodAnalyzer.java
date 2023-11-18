// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers.fluid.analyzers;

import jline.util.Matrix;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;

public interface MethodAnalyzer {

  void analyze(NetworkStruct sn, SolverOptions options, SolverResult result);

  Matrix getXVecIt();
}
