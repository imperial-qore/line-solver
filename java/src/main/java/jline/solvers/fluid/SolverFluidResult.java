// Copyright (c) 2012-2024, Imperial College London
// All rights reserved.

package jline.solvers.fluid;

import jline.solvers.SolverResult;
import jline.util.Matrix;

// Class to store results specific to SolverFluid
public class SolverFluidResult extends SolverResult {

  // Results from getAvg
  public Matrix odeStateVec;

  // Results from getCdfRespT and getTranCdfPassT
  public Matrix[][] distribC;
  public double distribRuntime;

  // Results from getProbAggr
  public double Pnir;
  public double logPnir;
}
