// Copyright (c) 2012-2022, Imperial College London
// All rights reserved.

package jline.solvers.env;


import jline.util.Maths;
import jline.lang.*;
import jline.lang.constant.SolverType;
import jline.solvers.EnsembleSolver;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.fluid.smoothing.PStarSearcher;
import jline.util.Matrix;
import jline.util.Utils;
import org.apache.commons.math3.optim.PointValuePair;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import static jline.lib.KPCToolbox.*;
// Solver for models immersed in a random environment
public class SolverEnv extends EnsembleSolver {

  // User-supplied representation of each stage transition
  private final Env envObj;
  private final NetworkStruct[] sn;
  private final Env.ResetQueueLengthsFunction[][] resetFromMarginal;
  private final Env.ResetEnvRatesFunction[][] resetEnvRates;

  public SolverEnv(Env renv, NetworkSolver[] solvers, SolverOptions options) {
    super(renv, "SolverEnv", options);
    int E = getNumberOfModels();
    this.envObj = renv;
    this.ensemble = renv.getEnsemble().toArray(new Network[0]);
    this.solvers = solvers;
    this.sn = new NetworkStruct[E];
    this.resetFromMarginal = new Env.ResetQueueLengthsFunction[E][E];
    this.resetEnvRates = new Env.ResetEnvRatesFunction[E][E];

    for (int e = 0; e < E; e++) {
      this.sn[e] = this.ensemble[e].getStruct(true);
      // TODO: implement the below method and then uncomment
      // if (!solvers[e].supports(ensemble[e])) {
      //   throw new RuntimeException("Model is not supported by the solver.");
      // }
      System.arraycopy(renv.resetQLFun[e], 0, resetFromMarginal[e], 0, E);
      System.arraycopy(renv.resetEnvRatesFun[e], 0, resetEnvRates[e], 0, E);
    }
  }

  // Convergence test at iteration it
  @Override
  protected boolean converged(int it) {

    int M = sn[0].nstations;
    int K = sn[0].nclasses;
    int E = getNumberOfModels();

    if (it <= 1) {
      return false;
    }
    Matrix mapes = new Matrix(1, E);
    for (int e = 0; e < E; e++) {
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < K; j++) {
          // Error is calculated only on entry value (t = 0)
          mapes.set(
              0,
              e,
              Math.max(
                  mapes.get(0, e),
                  Utils.mape(
                      Matrix.extractRows(results.get(it).get(e).QNt[i][j], 0, 1, null),
                      Matrix.extractRows(results.get(it - 1).get(e).QNt[i][j], 0, 1, null))));
        }
      }
    }
    return mapes.elementMax() < options.iter_tol;
  }

  @Override
  protected void init() {
    envObj.init();
  }

  @Override
  protected void pre(int it) {

    int E = getNumberOfModels();

    if (it == 1) {
      for (int e = 0; e < E; e++) {
        if (Double.isInfinite(solvers[e].options.timespan[1])) {
          this.solvers[e].getAvg();
        } else {
          this.solvers[e].getTranAvg();
        }
        this.ensemble[e].initFromMarginal(solvers[e].result.QN);
      }
    }
  }

  // Solves model in stage e
  @Override
  protected SolverResult analyze(int it, int e) {
    // TODO: [Qt,Ut,Tt] = self.ensemble{e}.getTranHandles;
    this.solvers[e].resetResults();
    // If pStar values exist, implement p-norm smoothing
    if (this.solvers[e].options.config.pstar.size() != 0) {
      PStarSearcher searcher = new PStarSearcher();
      Matrix targetQueueLengths = searcher.generateTargetQueueLengths(this.solvers[e].model);
      PointValuePair pStarValues =
          searcher.findPStarValues(this.solvers[e].model, targetQueueLengths);
      solvers[e].options.config.pstar.clear();
      for (int i = 0; i < this.solvers[e].model.getNumberOfNodes(); i++) {
        solvers[e].options.config.pstar.add(i, pStarValues.getPoint()[i]);
      }
    }
    this.solvers[e].getTranAvg();
    return this.solvers[e].result;
  }

  @Override
  protected void post(int it) {

    int M = sn[0].nstations;
    int K = sn[0].nclasses;
    int E = getNumberOfModels();
    Matrix[][] QExit = new Matrix[E][E];
    Matrix[][] UExit = new Matrix[E][E];
    Matrix[][] TExit = new Matrix[E][E];
    Matrix[][] w = new Matrix[E][E];
    Matrix[] QEntry = new Matrix[E]; // Average entry queue-length

    for (int e = 0; e < E; e++) {
      for (int h = 0; h < E; h++) {
        QExit[e][h] = new Matrix(M, K);
        UExit[e][h] = new Matrix(M, K);
        TExit[e][h] = new Matrix(M, K);
        Matrix D0 = envObj.proc[e][h].get(0);
        Matrix D1 = envObj.proc[e][h].get(1);
        map_normalize(D0,D1);
        for (int i = 0; i < M; i++) {
          for (int r = 0; r < K; r++) {
            w[e][h] = new Matrix(1, 1);
            Matrix cdf1 =
                map_cdf(
                        D0,D1,
                    Matrix.extractRows(
                        results.get(it).get(e).t, 1, results.get(it).get(e).t.getNumRows(), null));

            Matrix cdf2 =
                map_cdf(
                        D0,D1,
                    Matrix.extractRows(
                        results.get(it).get(e).t,
                        0,
                        results.get(it).get(e).t.getNumRows() - 1,
                        null));

            Matrix cdfDiff = cdf1.sub(1, cdf2);
            cdfDiff = cdfDiff.transpose();
            w[e][h] = Matrix.concatRows(w[e][h], cdfDiff, null);

            if (!w[e][h].hasNaN()) {
              QExit[e][h].set(
                  i,
                  r,
                  results.get(it).get(e).QNt[i][r].transpose().mult(w[e][h], null).get(0, 0)
                      / w[e][h].elementSum());
              UExit[e][h].set(
                  i,
                  r,
                  results.get(it).get(e).UNt[i][r].transpose().mult(w[e][h], null).get(0, 0)
                      / w[e][h].elementSum());
              TExit[e][h].set(
                  i,
                  r,
                  results.get(it).get(e).TNt[i][r].transpose().mult(w[e][h], null).get(0, 0)
                      / w[e][h].elementSum());
            }
          }
        }
      }
    }

    for (int e = 0; e < E; e++) {
      QEntry[e] = new Matrix(M, K);
      for (int h = 0; h < E; h++) {
        // Probability of coming from h to e \times resetFun(Qexit from h to e
        if (envObj.probOrig.get(h, e) > 0) {
          Matrix partialQEntry = new Matrix(0, 0);
          resetFromMarginal[h][e]
              .reset(QExit[h][e])
              .scale(envObj.probOrig.get(h, e), partialQEntry);
          QEntry[e] = QEntry[e].add(1, partialQEntry);
        }
      }
      solvers[e].resetResults();

      ensemble[e].initFromMarginal(QEntry[e]);
    }

    // Update transition rates between stages if State Dependent
    if (Objects.equals(options.method, "statedep")) {
      for (int e = 0; e < E; e++) {
        for (int h = 0; h < E; h++) {
          if (envObj.env[e][h] != null) {
            // If not defined, rates are left unchanged
            if (resetEnvRates[e][h] != null) {
              envObj.env[e][h] =
                  resetEnvRates[e][h].reset(
                      envObj.env[e][h], QExit[e][h], UExit[e][h], TExit[e][h]);
            }
          }
        }
      }
      // Reinitialise
      envObj.init();
    }
  }

  @Override
  protected void finish() {

    // Use last iteration
    int it = results.size();
    int M = sn[0].nstations;
    int K = sn[0].nclasses;
    int E = getNumberOfModels();
    Matrix[] QExit = new Matrix[E];
    Matrix[] UExit = new Matrix[E];
    Matrix[] TExit = new Matrix[E];
    Matrix[] w = new Matrix[E];

    for (int e = 0; e < E; e++) {
      QExit[e] = new Matrix(M, K);
      UExit[e] = new Matrix(M, K);
      TExit[e] = new Matrix(M, K);
      if (it > 0) {
        Matrix D0 = envObj.holdTime[e].get(0);
        Matrix D1 = envObj.holdTime[e].get(1);
        map_normalize(D0,D1);
        for (int i = 0; i < M; i++) {
          for (int r = 0; r < K; r++) {
            w[e] = new Matrix(1, 1);
            w[e].set(0, 0, 0);
            Matrix cdf1 =
                map_cdf(
                        D0, D1,
                    Matrix.extractRows(
                        results.get(it).get(e).t, 1, results.get(it).get(e).t.getNumRows(), null));

            Matrix cdf2 =
                map_cdf(
                        D0, D1,
                    Matrix.extractRows(
                        results.get(it).get(e).t,
                        0,
                        results.get(it).get(e).t.getNumRows() - 1,
                        null));

            Matrix cdfDiff = cdf1.sub(1, cdf2);
            cdfDiff = cdfDiff.transpose();
            w[e] = Matrix.concatRows(w[e], cdfDiff, new Matrix(0, 0));

            QExit[e].set(
                i,
                r,
                results
                        .get(it)
                        .get(e)
                        .QNt[i][r]
                        .transpose()
                        .mult(w[e], new Matrix(0, 0))
                        .get(0, 0)
                    / w[e].elementSum());
            UExit[e].set(
                i,
                r,
                results
                        .get(it)
                        .get(e)
                        .UNt[i][r]
                        .transpose()
                        .mult(w[e], new Matrix(0, 0))
                        .get(0, 0)
                    / w[e].elementSum());
            TExit[e].set(
                i,
                r,
                results
                        .get(it)
                        .get(e)
                        .TNt[i][r]
                        .transpose()
                        .mult(w[e], new Matrix(0, 0))
                        .get(0, 0)
                    / w[e].elementSum());
          }
        }
      }
    }

    Matrix QVal = QExit[0].clone();
    Matrix UVal = UExit[0].clone();
    Matrix TVal = TExit[0].clone();
    QVal.zero();
    UVal.zero();
    TVal.zero();
    for (int e = 0; e < E; e++) {
      Matrix partialResult = new Matrix(0, 0);
      QExit[e].scale(envObj.probEnv.get(e), partialResult);
      QVal = partialResult.add(1, QVal);
      UExit[e].scale(envObj.probEnv.get(e), partialResult);
      UVal = partialResult.add(1, UVal);
      TExit[e].scale(envObj.probEnv.get(e), partialResult);
      TVal = partialResult.add(1, TVal);
    }

    result.QN = QVal;
    result.UN = UVal;
    result.TN = TVal;
  }

  public String getName() {
    return "SolverEnv";
  }

  public void getGenerator() {
    // TODO: implementation - note return type should likely not be void
    throw new RuntimeException(
        "getGenerator() has not yet been implemented in JLINE - requires CTMC.");
  }

  public void getAvg() {
    getEnsembleAvg();
  }

  @Override
  public void getEnsembleAvg() {
    if (this.result.QN == null || this.result.QN.isEmpty() || this.options.force) {
      iterate();
    }
  }

  // Return table of average station metrics
  public void printAvgTable() {

    // TODO: add polymorphic version where 'keepDisabled' is a parameter
    boolean keepDisabled = false;

    this.getAvg();
    Matrix QN = this.result.QN;
    Matrix UN = this.result.UN;
    Matrix TN = this.result.TN;

    int M = QN.getNumRows();
    int K = QN.getNumCols();

    if (QN.isEmpty()) {
      throw new RuntimeException(
          "Unable to compute results and therefore unable to print AvgTable.");
    }

    if (!keepDisabled) {
      List<Double> Qval = new ArrayList<>();
      List<Double> Uval = new ArrayList<>();
      List<Double> Tval = new ArrayList<>();
      List<Double> respTVal = new ArrayList<>();
      List<String> className = new ArrayList<>();
      List<String> stationName = new ArrayList<>();
      for (int i = 0; i < M; i++) {
        for (int k = 0; k < K; k++) {
          if (QN.get(i, k) + UN.get(i, k) + TN.get(i, k) > 0) {
            Qval.add(QN.get(i, k));
            Uval.add(UN.get(i, k));
            Tval.add(TN.get(i, k));
            respTVal.add(QN.get(i, k) / TN.get(i, k));
            className.add(sn[0].jobclasses.get(k).getName());
            stationName.add(sn[0].stations.get(i).getName());
          }
        }
      }

      System.out.printf(
              "\n%-12s\t %-12s\t %-10s\t %-10s\t %-10s\t %-10s",
              "Station", "JobClass", "QLen", "Util", "RespT", "Tput");
      System.out.println(
              "\n------------------------------------------------------------------------------");
      NumberFormat nf = NumberFormat.getNumberInstance();
      nf.setMinimumFractionDigits(5);
      for (int i = 0; i < stationName.size(); i++) {
        System.out.format(
                "%-12s\t %-12s\t %-10s\t %-10s\t %-10s\t %-10s\n",
                stationName.get(i),
                className.get(i),
                nf.format(Qval.get(i)),
                nf.format(Uval.get(i)),
                nf.format(respTVal.get(i)),
                nf.format(Tval.get(i)));
      }
      System.out.println(
              "------------------------------------------------------------------------------");
    } else {
      // TODO: implementation if keepDisabled is set to true
      System.out.println("Warning: unimplemented code reached in SolverEnv.printAvgTable.");
    }
  }

  @Override
  public boolean supports(Ensemble model) {
    // TODO: implementation
    throw new RuntimeException("supports() has not yet been implemented in JLINE.");
  }

  public static boolean supports(Model model) {
    // TODO: implementation
    throw new RuntimeException("supports() has not yet been implemented in JLINE.");
  }

  @Override
  public void runAnalyzer() {
    iterate();
  }

  public static SolverOptions defaultOptions() {
    return new SolverOptions(SolverType.Env);
  }
}
