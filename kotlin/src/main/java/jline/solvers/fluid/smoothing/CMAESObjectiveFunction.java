// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers.fluid.smoothing;

import jline.util.Matrix;
import jline.lang.Network;
import jline.solvers.fluid.SolverFluid;
import org.apache.commons.math3.analysis.MultivariateFunction;

public class CMAESObjectiveFunction implements MultivariateFunction {

  private final Matrix targetQueueLengths;
  private final Network model;
  private final boolean stiff;
  private int evaluation;

  public CMAESObjectiveFunction(Matrix targetQueueLengths, Network model, boolean stiff) {
    this.targetQueueLengths = targetQueueLengths;
    this.model = model;
    this.stiff = stiff;
    this.evaluation = 0;
  }

  private double computeErrorValue(Matrix QNFluid) {
    // Compute Error (L2 Norm - (Squared) Euclidean)
    double errorValue = 0;
    int M = targetQueueLengths.getNumRows();
    int K = targetQueueLengths.getNumCols();
    for (int i = 0; i < M; i++) {
      for (int k = 0; k < K; k++) {
        errorValue += Math.pow(QNFluid.get(i, k) - targetQueueLengths.get(i, k), 2);
      }
    }
    return errorValue;
  }

  @Override
  public double value(double[] doubles) {

    /*    System.out.format("Evaluation %d\n", this.evaluation);
    for (int i = 0; i < doubles.length; i++) {
      System.out.format("pStar for Station %d: %f\n", i, doubles[i]);
    }*/

    SolverFluid solverFluid = new SolverFluid(this.model);
    solverFluid.options.stiff = this.stiff;
    for (int i = 0; i < doubles.length; i++) {
      solverFluid.options.config.pstar.add(i, doubles[i]);
    }
    solverFluid.runAnalyzer();
    Matrix QNFluid = solverFluid.result.QN;
    double errorValue = computeErrorValue(QNFluid);

    // System.out.format("Error Value: %f\n\n", errorValue);

    evaluation++;
    return errorValue;
  }
}
