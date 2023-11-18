// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers.fluid.analyzers;

import jline.lang.constant.GlobalConstants;
import jline.lang.constant.VerboseLevel;
import jline.util.Matrix;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.JobClassType;
import jline.lang.constant.NodeType;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.fluid.odes.MatrixMethodODE;
import jline.solvers.fluid.odes.TransientDataHandler;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.equation.Equation;

import java.util.*;

import static java.lang.Double.*;
import static java.lang.Math.abs;
import static java.lang.Math.min;

import odesolver.LSODA;

public class MatrixMethodAnalyzer implements MethodAnalyzer {

  public Matrix xvec_t;
  public Matrix xvec_it;

  @Override
  public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {

    int M = sn.nstations; // Number of stations
    int K = sn.nclasses; // Number of classes

    Matrix S = sn.nservers.clone();
    double initialPopulation = sn.njobs.elementSum();
    int SRows = S.getNumRows();
    for (int i = 0; i < SRows; i++) {
      if (isInfinite(S.get(i, 0))) {
        S.set(i, 0, initialPopulation);
      }
    }

    // W (complete graph of transition rates) as per Ruuskanen et al., PEVA 151 (2021)
    Matrix psi = new Matrix(0, 0);
    Matrix A = new Matrix(0, 0);
    Matrix B = new Matrix(0, 0);
    Matrix lambda = new Matrix(M * K, 1);

    for (int i = 0; i < M; i++) {
      if (sn.nodetypes.get((int) sn.stationToNode.get(i)) == NodeType.Source) {
        for (int k = 0; k < K; k++) {
          if (sn.jobclasses.get(k).getJobClassType() == JobClassType.Open) {
            lambda.set(i * K + k, 0, sn.rates.get(i, k));
          }
        }
      }
    }

    for (int i = 0; i < M; i++) {
      Station station = sn.stations.get(i);
      for (int r = 0; r < K; r++) {
        JobClass jobClass = sn.jobclasses.get(r);
        if (sn.phases.get(i, r) == 0) {
          Matrix zeroMatrix = new Matrix(1, 1, 1);
          Matrix nanMatrix = new Matrix(1, 1, 1);
          nanMatrix.set(0, 0, Double.NaN);
          psi = psi.createBlockDiagonal(zeroMatrix);
          A = A.createBlockDiagonal(nanMatrix);
          B = B.createBlockDiagonal(zeroMatrix);
        } else {
          psi = psi.createBlockDiagonal(sn.proc.get(station).get(jobClass).get(0));
          A = A.createBlockDiagonal(sn.pie.get(station).get(jobClass).transpose());
          B = B.createBlockDiagonal(sn.proc.get(station).get(jobClass).get(1).sumRows());
        }
      }
    }

    // W (complete graph of transition rates) as per Ruuskanen et al., PEVA 151 (2021)
    Matrix W = calculateW(sn, psi, A, B);

    Matrix ALambda = A.mult(lambda, null);

    int totalNumPhases = (int) sn.phases.elementSum();
    // State mapping to queues (called Q(a) in Ruuskanen et al.)
    Matrix Qa = new Matrix(1, totalNumPhases);
    // To compute per-class queue length, utilisation and throughput at the end
    Matrix SQC = new Matrix(M * K, totalNumPhases);
    Matrix SUC = new Matrix(M * K, totalNumPhases);
    Matrix STC = new Matrix(M * K, totalNumPhases);
    // To compute total queue length in ODEs
    Matrix SQ = new Matrix(totalNumPhases, totalNumPhases);

    int state = 0;
    for (int i = 0; i < M; i++) {
      Station station = sn.stations.get(i);
      for (int r = 0; r < K; r++) {
        JobClass jobClass = sn.jobclasses.get(r);
        int nPhases = (int) sn.phases.get(i, r);
        for (int k = 0; k < nPhases; k++) {
          Qa.set(0, state, i);
          SQC.set(i * K + r, state, 1);
          SUC.set(i * K + r, state, 1 / S.get(i, 0));
          STC.set(i * K + r, state, sn.proc.get(station).get(jobClass).get(1).sumRows(k));
          state++;
        }
      }
    }

    int nextSQRow = 0;
    for (int i = 0; i < M; i++) {
      for (int r = 0; r < K; r++) {
        int nPhases = (int) sn.phases.get(i, r);
        for (int k = 0; k < nPhases; k++) {
          for (int col = 0; col < totalNumPhases; col++) {
            if (Qa.get(0, col) == i) {
              SQ.set(nextSQRow, col, 1); // Setting weights
            }
          }
          nextSQRow++;
        }
      }
    }

    Matrix Sa = new Matrix(totalNumPhases, 1);
    for (int i = 0; i < totalNumPhases; i++) {
      Sa.set(i, 0, S.get((int) Qa.get(0, i), 0));
    }

    double minNonZeroRate = POSITIVE_INFINITY;
    int WRows = W.getNumRows();
    int WCols = W.getNumCols();
    for (int i = 0; i < WRows; i++) {
      for (int j = 0; j < WCols; j++) {
        double tmpRate = abs(W.get(i, j));
        if (tmpRate < minNonZeroRate && tmpRate > 0) {
          minNonZeroRate = tmpRate;
        }
      }
    }
    double[] tRange = {
      options.timespan[0], min(options.timespan[1], abs(10 * options.iter_max / minNonZeroRate))
    };

    int initSolLength = options.init_sol.length();
    double[] initialState = new double[initSolLength];
    double[] nextState = new double[initSolLength];
    for (int i = 0; i < initSolLength; i++) {
      initialState[i] = options.init_sol.get(0, i);
      nextState[i] = 0;
    }

    // Choose between original compact matrix form representation, and p-norm smoothed
    // representation as per Ruuskanen et al., PEVA 151 (2021).
    FirstOrderDifferentialEquations ode;

    if (options.config.pstar.size() == 0) { // If pStar values do not exist
      ode = new MatrixMethodODE(W, SQ, S, Qa, ALambda, initSolLength);
    } else {
      ode = new MatrixMethodODE(W, SQ, S, Qa, ALambda, initSolLength, sn, options.config.pstar);
    }

    int Tmax;
    if (options.stiff) {
      LSODA odeSolver;
      if (options.tol > GlobalConstants.CoarseTol){
        odeSolver = options.odesolvers.fastStiffODESolver;
      } else {
        odeSolver = options.odesolvers.accurateStiffODESolver;
      }
      //System.out.print("Starting ODE integration cycle...");
      odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState);
      //System.out.println("done.");
      Tmax = odeSolver.getStepsTaken()+1;
      Matrix tVec = new Matrix(Tmax, 1);
      Matrix xVec = new Matrix(Tmax,ode.getDimension());
      ArrayList<Double> tHistory = odeSolver.getTvec();
      ArrayList<Double[]> yHistory = odeSolver.getYvec();
      for (int i = 0; i < Tmax; i++) {
        tVec.set(i,0,tHistory.get(i));
        for (int j = 0; j < ode.getDimension(); j++)
          xVec.set(i,j,Math.max(0,yHistory.get(i)[j]));
      }
      result.t = tVec;
      this.xvec_t = xVec;
    }
    else {
      if (options.tol > GlobalConstants.CoarseTol && (options.verbose == VerboseLevel.DEBUG)) {
        System.err.println(
                "Fast, non-stiff ODE solver is not yet available in JLINE. Using accurate non-stiff ODE solver instead.");
      }
      FirstOrderIntegrator odeSolver = options.odesolvers.accurateODESolver;
      odeSolver.clearStepHandlers();
      TransientDataHandler stepHandler = new TransientDataHandler(initSolLength);
      odeSolver.addStepHandler(stepHandler);
      //System.out.print("Starting ODE integration cycle...");
      odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState);
      //System.out.println("done.");

      // Retrieve Transient Data
      result.t = stepHandler.tVec;
      Tmax = result.t.getNumRows();
      this.xvec_t = stepHandler.xVec;
    }
    this.xvec_it = Matrix.extractRows(xvec_t, Tmax - 1, Tmax, null);

    // Use Transient Data to Store Results
    result.QNt = new Matrix[M][K];
    result.UNt = new Matrix[M][K];
    result.TNt = new Matrix[M][K];
    for (int i = 0; i < M; i++) {
      for (int k = 0; k < K; k++) {
        result.QNt[i][k] = new Matrix(Tmax, 1);
        result.UNt[i][k] = new Matrix(Tmax, 1);
        result.TNt[i][k] = new Matrix(Tmax, 1);
      }
    }

    for (int step = 0; step < Tmax; step++) {

      Matrix x = Matrix.extractRows(xvec_t, step, step + 1, null);
      x = x.transpose();
      Matrix theta = x.clone(); // Theta per Ruuskanen et al., PEVA 151 (2021).
      Matrix SQx = SQ.mult(x, null);
      for (int phase = 0; phase < totalNumPhases; phase++) {
        double valSQx = SQx.get(phase, 0) + GlobalConstants.Zero;
        double valSa = Sa.get(phase, 0);
        if (valSQx > valSa) {
          theta.set(phase, 0, x.get(phase, 0) / valSQx * valSa);
        }
      }

      Matrix QNtmp = SQC.mult(x, null);
      Matrix UNtmp = SUC.mult(theta, null);
      Matrix TNtmp = STC.mult(theta, null);

      for (int i = 0; i < M; i++) {
        for (int k = 0; k < K; k++) {
          result.QNt[i][k].set(step, 0, QNtmp.get(i * K + k, 0));
          result.UNt[i][k].set(step, 0, UNtmp.get(i * K + k, 0));
          result.TNt[i][k].set(step, 0, TNtmp.get(i * K + k, 0));
        }
      }
    }

    result.QN = new Matrix(M, K);
    result.UN = new Matrix(M, K);
    result.RN = new Matrix(M, K);
    result.TN = new Matrix(M, K);
    result.WN = new Matrix(M, K); // TODO
    result.AN = new Matrix(M, K); // TODO
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < K; j++) {
        result.QN.set(i, j, result.QNt[i][j].get(Tmax - 1, 0));
        result.UN.set(i, j, result.UNt[i][j].get(Tmax - 1, 0));
        result.TN.set(i, j, result.TNt[i][j].get(Tmax - 1, 0));
        result.RN.set(i, j, result.QN.get(i, j) / result.TN.get(i, j));
      }
    }
  }

  @Override
  public Matrix getXVecIt() {
    return this.xvec_it;
  }

  private Matrix calculateW(NetworkStruct sn, Matrix psi, Matrix A, Matrix B) {

    // ODE building as per Ruuskanen et al., PEVA 151 (2021).
    DMatrixSparseCSC psiDMatrix = psi.toDMatrixSparseCSC();
    DMatrixSparseCSC ADMatrix = A.toDMatrixSparseCSC();
    DMatrixSparseCSC BDMatrix = B.toDMatrixSparseCSC();
    DMatrixSparseCSC P = sn.rt.toDMatrixSparseCSC();

    Equation calculateW = new Equation();
    calculateW.alias(psiDMatrix, "psi", ADMatrix, "A", P, "P", BDMatrix, "B");
    calculateW.process("W = psi + B*P*A'");
    Matrix W = new Matrix(calculateW.lookupSimple("W"));

    // Remove disabled transitions
    if (W.hasNaN()) {
      int WRows = W.getNumRows();
      int WCols = W.getNumCols();
      List<Integer> colsWithNans = new LinkedList<>();
      Matrix sumWCols = W.sumCols();
      for (int i = 0; i < WCols; i++) {
        if (isNaN(sumWCols.get(0, i))) {
          colsWithNans.add(i);
        }
      }

      Matrix tmpW =
          new Matrix(
              W.getNumRows() - colsWithNans.size(), W.getNumCols() - colsWithNans.size());
      int tmpWRow = 0;
      for (int i = 0; i < WRows; i++) {
        while (colsWithNans.contains(i)) {
          i++;
        }
        int tmpWCol = 0;
        for (int j = 0; j < WCols; j++) {
          while (colsWithNans.contains(j)) {
            j++;
          }
          tmpW.set(tmpWRow, tmpWCol, W.get(i, j));
          tmpWCol++;
        }
        tmpWRow++;
      }

      W = tmpW;
    }

    return W;
  }
}
