// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers.fluid.odes;

import jline.lang.constant.GlobalConstants;
import jline.util.Matrix;
import jline.lang.NetworkStruct;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.equation.Equation;
import org.ejml.simple.SimpleMatrix;

import java.util.List;

import static java.lang.Math.min;

public class MatrixMethodODE implements FirstOrderDifferentialEquations {

  private final DMatrixSparseCSC W;
  private final DMatrixSparseCSC SQ;
  private final DMatrixSparseCSC S;
  private final DMatrixSparseCSC Qa;
  private final DMatrixSparseCSC ALambda;
  private final int numDimensions;
  private DMatrixSparseCSC pQa;

  public MatrixMethodODE(
          Matrix W, Matrix SQ, Matrix S, Matrix Qa, Matrix ALambda, int numDimensions) {
    this.W = W.toDMatrixSparseCSC();
    this.SQ = SQ.toDMatrixSparseCSC();
    this.S = S.toDMatrixSparseCSC();
    this.Qa = Qa.toDMatrixSparseCSC();
    this.ALambda = ALambda.toDMatrixSparseCSC();
    this.numDimensions = numDimensions;
    this.pQa = new DMatrixSparseCSC(0, 0);
  }

  public MatrixMethodODE(
      Matrix W,
      Matrix SQ,
      Matrix S,
      Matrix Qa,
      Matrix ALambda,
      int numDimensions,
      NetworkStruct sn,
      List<Double> pStarValues) {

    this(W, SQ, S, Qa, ALambda, numDimensions);

    this.pQa = new DMatrixSparseCSC(SQ.getNumRows(), 1);
    int row = 0;
    for (int i = 0; i < sn.nstations; i++) {
      double pStarValue = pStarValues.get(i);
      for (int j = 0; j < sn.nclasses; j++) {
        int nPhases = (int) sn.phases.get(i, j);
        for (int k = 0; k < nPhases; k++) {
          pQa.set(row, 0, pStarValue);
          row++;
        }
      }
    }
  }

  @Override
  public int getDimension() {
    return numDimensions;
  }

  @Override
  public void computeDerivatives(double t, double[] x, double[] dxdt)
      throws MaxCountExceededException, DimensionMismatchException {

    DMatrixSparseCSC xDMS = new DMatrixSparseCSC(x.length, 1);
    for (int i = 0; i < x.length; i++) {
      xDMS.set(i, 0, x[i]);
    }

    Equation calculateSumXQa = new Equation();
    calculateSumXQa.alias(xDMS, "x", SQ, "SQ", GlobalConstants.FineTol, "distribZero");
    calculateSumXQa.process("sumXQa = distribZero + SQ * x");
    SimpleMatrix sumXQa = calculateSumXQa.lookupSimple("sumXQa");

    int QaCols = this.Qa.getNumCols();
    DMatrixSparseCSC SQa = new DMatrixSparseCSC(QaCols, 1);
    for (int i = 0; i < QaCols; i++) {
      SQa.set(i, 0, S.get((int) Qa.get(0, i), 0));
    }

    SimpleMatrix dxdtTmp;
    if (this.pQa.getNumRows() == 0) { // If no pStar values have been specified
      dxdtTmp = computeDerivativesWithoutSmoothing(xDMS, sumXQa, SQa);
    } else {
      dxdtTmp = computeDerivativesUsingPNormSmoothing(xDMS, sumXQa, SQa);
    }

    for (int i = 0; i < dxdt.length; i++) {
      dxdt[i] = dxdtTmp.get(i);
    }
  }

  private SimpleMatrix computeDerivativesWithoutSmoothing(
      DMatrixSparseCSC x, SimpleMatrix sumXQa, DMatrixSparseCSC SQa) {

    int sumXQaRows = sumXQa.numRows();
    int sumXQaCols = sumXQa.numCols();
    DMatrixSparseCSC minOfSumXQaAndSQa = new DMatrixSparseCSC(sumXQaRows, sumXQaCols);
    for (int i = 0; i < sumXQaRows; i++) {
      for (int j = 0; j < sumXQaCols; j++) {
        minOfSumXQaAndSQa.set(i, j, min(sumXQa.get(i, j), SQa.get(i, j)));
      }
    }

    Equation computeDerivatives = new Equation();
    computeDerivatives.alias(
        W, "W", x, "x", sumXQa, "sumXQa", minOfSumXQaAndSQa, "minOfSumXQaAndSQa", ALambda, "ALambda");
    computeDerivatives.process("dxdt = W' * (x ./ sumXQa .* minOfSumXQaAndSQa) + ALambda");
    return computeDerivatives.lookupSimple("dxdt");
  }

  private SimpleMatrix computeDerivativesUsingPNormSmoothing(
      DMatrixSparseCSC x, SimpleMatrix sumXQa, DMatrixSparseCSC SQa) {

    // ghat = smoothed processor-share constraint approximation, per Ruuskanen et. al
    DMatrixSparseCSC ghat = x.createLike();
    for (int i = 0; i < x.getNumRows(); i++) {
      // x, c and p as per Ruuskanen's Julia implementation
      double xVal = sumXQa.get(i, 0);
      double cVal = SQa.get(i, 0);
      double pVal = pQa.get(i, 0);
      double ghatVal = 1 / Math.pow(1 + Math.pow(xVal / cVal, pVal), 1 / pVal);
      if (Double.isNaN(ghatVal)) {
        ghat.set(i, 0, 0);
      } else {
        ghat.set(i, 0, ghatVal);
      }
    }

    Equation computeDerivatives = new Equation();
    computeDerivatives.alias(W, "W", x, "x", ghat, "ghat", ALambda, "ALambda");
    computeDerivatives.process("dxdt = W' * (x .* ghat) + ALambda");
    return computeDerivatives.lookupSimple("dxdt");
  }
}
