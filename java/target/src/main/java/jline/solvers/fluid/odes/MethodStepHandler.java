// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers.fluid.odes;

import jline.util.Matrix;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

public class MethodStepHandler implements StepHandler {

  public Matrix tVec;
  public Matrix xVec;
  private int maxSteps;
  private int stepCount;

  public MethodStepHandler(int numDimensions) {

    maxSteps = 2; // Can be increased but must be 2 or more
    tVec = new Matrix(maxSteps, 1);
    xVec = new Matrix(maxSteps, numDimensions);
    stepCount = 0;
  }

  @Override
  public void init(double t0, double[] x0, double t) {

    tVec.set(stepCount, 0, t0);
    for (int i = 0; i < x0.length; i++) {
      xVec.set(stepCount, i, x0[i]);
    }
    stepCount++;
  }

  @Override
  public void handleStep(StepInterpolator interpolator, boolean isLast) {

    tVec.set(stepCount, 0, interpolator.getCurrentTime());

    double[] x = interpolator.getInterpolatedState();
    for (int i = 0; i < x.length; i++) {
      xVec.set(stepCount, i, Math.max(0, x[i])); // Equivalent of NonNegative odeset in LINE
    }
    stepCount++;

    if (isLast) {
      // Shrinking arrays to minimum number of rows
      tVec.numRows = stepCount;
      xVec.numRows = stepCount;
      return;
    }

    // Dynamic resizing of JLineMatrix objects per Dynamic Array logic
    if (stepCount == maxSteps) {
      maxSteps *= 2;
      tVec.expandMatrix(maxSteps, tVec.getNumCols(), maxSteps * tVec.getNumCols());
      xVec.expandMatrix(maxSteps, xVec.getNumCols(), maxSteps * xVec.getNumCols());
    }
  }
}
