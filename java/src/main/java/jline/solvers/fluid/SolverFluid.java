// Copyright (c) 2012-2024, Imperial College London
// All rights reserved.

package jline.solvers.fluid;

import jline.lang.constant.*;
import jline.lang.state.State;
import jline.util.PopulationLattice;
import jline.util.Matrix;
import jline.lang.JobClass;
import jline.lang.Model;
import jline.lang.Network;
import jline.lang.distributions.Coxian;
import jline.solvers.NetworkSolver;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.fluid.analyzers.ClosingAndStateDepMethodsAnalyzer;
import jline.solvers.fluid.analyzers.MatrixMethodAnalyzer;
import jline.solvers.fluid.analyzers.MethodAnalyzer;
import jline.solvers.fluid.odes.ClosingAndStateDepMethodsODE;
import jline.solvers.fluid.odes.TransientDataHandler;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;

import java.util.*;

import static java.lang.Double.*;
import static java.lang.Math.abs;
import static java.lang.Math.min;

// SolverFluid is based on fluid and mean-field approximation methods.
public class SolverFluid extends NetworkSolver {

  public SolverFluid(Network model) {
    this(model, SolverFluid.defaultOptions());
    this.result = new SolverFluidResult();
  }

  public SolverFluid(Network model, String method) {
    super(model, "SolverFluid", SolverFluid.defaultOptions().method(method));
    this.result = new SolverFluidResult();
  }

  public SolverFluid(Network model, SolverOptions options) {
    super(model, "SolverFluid", options);
    // TODO: self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
    this.result = new SolverFluidResult();
  }

  public void initSol() {

    Matrix initSol = new Matrix(1, 0);
    for (int ind = 0; ind < sn.nnodes; ind++) {
      if (sn.isstateful.get(ind, 0) == 1) {
        int isf = (int) sn.nodeToStateful.get(0, ind);
        int ist = (int) sn.nodeToStation.get(0, ind);
        Matrix state_i = new Matrix(1, 0);
        // Compared to state_i, initSol_i does not track disabled classes
        // and removes Inf entries in the Sources
        Matrix initSol_i = new Matrix(1, 0);

        State.StateMarginalStatistics stats =
            State.toMarginal(
                sn, ind, sn.state.get(sn.stations.get(isf)), null, null, null, null, null);
        Matrix nir = stats.nir;
        List<Matrix> kir_i = stats.kir;

        switch (sn.sched.get(sn.stations.get(ist))) {
          case EXT:
            state_i.expandMatrix(1, state_i.getNumCols() + 1, state_i.getNumElements() + 1);
            state_i.set(0, 0, POSITIVE_INFINITY); // Fluid does not model infinite buffer?
            int rMax = kir_i.get(0).getNumCols();
            for (int r = 0; r < rMax; r++) {
              int kMax = sn.mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).length();
              for (int k = 0; k < kMax; k++) {
                state_i.expandMatrix(1, state_i.getNumCols() + 1, state_i.getNumElements() + 1);
                state_i.set(0, state_i.getNumCols() - 1, kir_i.get(k).get(0, r));
                if (!isNaN(sn.rates.get(ist, r))) {
                  initSol_i.expandMatrix(
                      1, initSol_i.getNumCols() + 1, initSol_i.getNumElements() + 1);
                  initSol_i.set(0, initSol_i.getNumCols() - 1, kir_i.get(k).get(0, r));
                }
              }
            }
            break;

          case FCFS:
          case SIRO:
          case PS:
          case INF:
          case DPS:
          case HOL:
            rMax = kir_i.get(0).getNumCols();
            for (int r = 0; r < rMax; r++) {
              int kMax = sn.mu.get(sn.stations.get(ist)).get(sn.jobclasses.get(r)).length();
              for (int k = 0; k < kMax; k++) {
                state_i.expandMatrix(1, state_i.getNumCols() + 1, state_i.getNumElements() + 1);
                if (!isNaN(sn.rates.get(ist, r))) {
                  initSol_i.expandMatrix(
                      1, initSol_i.getNumCols() + 1, initSol_i.getNumElements() + 1);
                }
                if (k == 0) {
                  double sumKir_i = 0;
                  for (int m = 1; m < kir_i.size(); m++) {
                    sumKir_i = kir_i.get(m).get(0, r);
                  }
                  // Jobs in waiting buffer are re-started phase 1
                  state_i.set(0, state_i.getNumCols() - 1, nir.get(0, r) - sumKir_i);
                  if (!isNaN(sn.rates.get(ist, r))) {
                    initSol_i.set(0, initSol_i.getNumCols() - 1, nir.get(0, r) - sumKir_i);
                  }
                } else {
                  state_i.set(0, state_i.getNumCols() - 1, kir_i.get(k).get(0, r));
                  if (!isNaN(sn.rates.get(ist, r))) {
                    initSol_i.set(0, initSol_i.getNumCols() - 1, kir_i.get(k).get(0, r));
                  }
                }
              }
            }
            break;
          default:
            System.err.format("Unsupported scheduling policy at station %d.", ist);
            return;
        }

        initSol = Matrix.concatColumns(initSol, initSol_i, null);
        sn.state.put(sn.stations.get(isf), state_i);
      }
    }
    options.init_sol = initSol;
  }

  private void runMethodSpecificAnalyzer() {

    int M = sn.nstations;
    int K = sn.nclasses;
    Matrix phases = sn.phases.clone();
    Matrix phasesLast = sn.phases.clone();
    Matrix rates0 = sn.rates.clone();

    Matrix V = new Matrix(M, K);
    for (int i = 0; i < sn.visits.size(); i++) {
      V = V.add(1, sn.visits.get(i));
    }

    MethodAnalyzer analyzer;
    switch (options.method) {
      case "statedep":
      case "closing":
        analyzer = new ClosingAndStateDepMethodsAnalyzer();
        break;
      case "matrix":
        analyzer = new MatrixMethodAnalyzer();
        break;
      default:
        // TODO: not implemented
        throw new RuntimeException("Unsupported method.");
    }

    if (options.init_sol.isEmpty()) {
      initSol();
    }
    analyzer.analyze(sn, options, result);

    // Single iteration is sufficient for statedep, so do nothing unless matrix or closing
    if ((Objects.equals(options.method, "matrix")) || (Objects.equals(options.method, "closing"))) {
      if (sn.sched.containsValue(SchedStrategy.FCFS)) {
        int iter = 0;
        Matrix eta_1 = new Matrix(M, 1);
        Matrix eta = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
          eta.set(i, 0, POSITIVE_INFINITY);
        }
        double tol = 0.01;

        while (iter <= options.iter_max) {
          double maxAbsEtaEta_1 = NEGATIVE_INFINITY;
          for (int i = 0; i < M; i++) {
            double val = Math.abs(1 - (eta.get(i, 0) / eta_1.get(i, 0)));
            if (val > maxAbsEtaEta_1) {
              maxAbsEtaEta_1 = val;
            }
          }
          if (maxAbsEtaEta_1 <= tol) {
            break;
          }

          iter++;
          eta_1 = eta.clone();
          for (int i = 0; i < M; i++) {
            LinkedList<Integer> sdCols = new LinkedList<>();
            for (int k = 0; k < K; k++) {
              if (rates0.get(i, k) > 0) {
                sdCols.add(k);
              }
            }
            for (int k : sdCols) {
              result.UN.set(i, k, result.TN.get(i, k) / rates0.get(i, k));
            }
          }

          Matrix ST0 = new Matrix(M, K);
          for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
              ST0.set(i, j, 1 / rates0.get(i, j));
              if (isInfinite(ST0.get(i, j))) {
                ST0.set(i, j, 1/GlobalConstants.Zero);
              } else if (isNaN(ST0.get(i, j))) {
                ST0.set(i, j, GlobalConstants.Zero);
              }
            }
          }

          Matrix XN = new Matrix(1, K);
          for (int k = 0; k < K; k++) {
            if (sn.refstat.get(k, 0) >= 0) { // Ignore artificial classes
              XN.set(0, k, result.TN.get((int) sn.refstat.get(k, 0), k));
            }
          }

          if (Objects.equals(options.config.highvar, "interp")) {
            // TODO: implement npfqn_nonexp_approx
            throw new RuntimeException("Warning: unimplemented code reached in SF.rMSA 1");
          }
          // DELETE THIS SECTION ONCE YOU'VE IMPLEMENTED NPFQN //
          eta = new Matrix(M, 1);
          eta.ones();
          Matrix ST = ST0.clone();
          // DELETE DOWN TO HERE

          Matrix rates = new Matrix(M, K);
          for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
              rates.set(i, j, 1 / ST.get(i, j));
              if (isInfinite(rates.get(i, j))) {
                rates.set(i, j, 1/GlobalConstants.Zero);
              } else if (isNaN(rates.get(i, j))) {
                rates.set(i, j, GlobalConstants.Zero);
              }
            }
          }

          for (int i = 0; i < M; i++) {
            if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.FCFS) {
              for (int k = 0; k < K; k++) {
                if (rates.get(i, k) > 0 && sn.scv.get(i, k) > 0) {
                  Coxian cx = Coxian.fitMeanAndSCV(1 / rates.get(i, k), sn.scv.get(i, k));
                  Matrix muik = cx.getMu();
                  Matrix phiik = cx.getPhi();
                  // we now handle the case that due to either numerical issues or different
                  // relationship between scv and mean if the size of the phase-type representation
                  // has changed
                  phases.set(i, k, muik.length());
                  State.StateMarginalStatistics stats = null;
                  if (phases.get(i, k) != phasesLast.get(i, k)) { // If number of phases changed
                    // Before we update sn we adjust the initial state
                    int isf = (int) sn.stationToStateful.get(0, i);
                    stats =
                        State.toMarginal(
                            sn,
                            i,
                            sn.state.get(sn.stations.get(isf)),
                            null,
                            null,
                            null,
                            null,
                            null);
                  }
                  sn.proc.get(sn.stations.get(i)).put(sn.jobclasses.get(k), cx.getRepres());
                  sn.mu.get(sn.stations.get(i)).put(sn.jobclasses.get(k), muik);
                  sn.phi.get(sn.stations.get(i)).put(sn.jobclasses.get(k), phiik);
                  sn.phases = phases.clone();

                  sn.phasessz = sn.phases.clone();
                  for (int row = 0; row < M; row++) {
                    for (int col = 0; col < K; col++) {
                      if (sn.phasessz.get(row, col) < 1) {
                        sn.phasessz.set(row, col, 1);
                      }
                    }
                  }

                  sn.phaseshift = new Matrix(sn.phases.getNumRows(), 1);
                  sn.phaseshift =
                      Matrix.concatColumns(sn.phaseshift, sn.phasessz.cumsumViaRow(), null);

                  if (phases.get(i, k) != phasesLast.get(i, k)) {
                    int isf = (int) sn.stationToStateful.get(0, i);
                    // We now initialise the new service process
                    assert stats != null;
                    sn.state.put(
                        sn.stations.get(isf),
                        State.fromMarginalAndStarted(sn, i, stats.nir, stats.sir));
                    // Pick one as the marginals won't change
                    sn.state.put(
                        sn.stations.get(isf),
                        Matrix.extractRows(sn.state.get(sn.stations.get(isf)), 0, 1, null));
                  }
                }
              }
            }

            options.init_sol = analyzer.getXVecIt();
            options.init_sol.transpose();
            // If there is a change of phases reset
            boolean callInitSol = false;
            for (int row = 0; row < M; row++) {
              for (int col = 0; col < K; col++) {
                if (phasesLast.get(row, col) - phases.get(row, col) != 0) {
                  callInitSol = true;
                }
              }
            }
            if (callInitSol) {
              initSol();
            }
          }
          sn.phases = phases.clone();
          analyzer.analyze(sn, options, result);
          phasesLast = phases.clone();
        } // FCFS iteration ends here

        // The FCFS iteration reinitializes at the solution of the last iterative step. We now
        // have converged in the substitution of the model parameters and we rerun everything from
        // the true initial point so that we get the correct transient.
        initSol();
        analyzer.analyze(sn, options, result);
      }
    }

    if (result.t.get(0, 0) == 0) {
      result.t.set(0, 0, 0.00000001);
    }

    Matrix Ufull0 = result.UN.clone();
    for (int i = 0; i < M; i++) {
      List<Integer> sdCols = new LinkedList<>();
      for (int k = 0; k < K; k++) {
        if (result.QN.get(i, k) > 0) {
          sdCols.add(k);
        }
        if (result.QN.get(i, k) == 0) {
          result.UN.set(i, k, 0);
          result.RN.set(i, k, 0);
        }
      }

      if (sn.sched.get(sn.stations.get(i)) == SchedStrategy.INF) {
        for (int k : sdCols) {
          result.UN.set(i, k, result.QN.get(i, k));
          result.UNt[i][k] = result.QNt[i][k].clone();
          result.UNt[i][k].scale(sn.rates.get(i, k), result.TNt[i][k]);
        }
      } else {
        double sumUfull0sd = 0;
        double sumTNDivRates0sd = 0;
        for (int k : sdCols) {
          sumUfull0sd += Ufull0.get(i, k);
          sumTNDivRates0sd += result.TN.get(i, k) / rates0.get(i, k);
        }
        for (int k : sdCols) {
          // correct for the real rates, instead of the diffusion approximation rates
          result.UN.set(
              i,
              k,
              min(
                  1,
                  min(
                      result.QN.get(i, k) / sn.nservers.get(i, 0),
                      sumUfull0sd * result.TN.get(i, k) / rates0.get(i, k) / sumTNDivRates0sd)));
          result.UNt[i][k].scale(sn.rates.get(i, k) * sn.nservers.get(i, 0), result.TNt[i][k]);
          result.RN.set(i, k, (result.QN.get(i, k) / result.TN.get(i, k)));
        }
      }
    }

    for (int i = 0; i < M; i++) {
      for (int j = 0; j < K; j++) {
        if (Double.isNaN(result.UN.get(i, j))) {
          result.UN.set(i, j, 0);
        }
        if (Double.isNaN(result.RN.get(i, j))) {
          result.RN.set(i, j, 0);
        }
      }
    }

    result.XN = new Matrix(1, K);
    result.CN = new Matrix(1, K);
    for (int k = 0; k < K; k++) {
      if (sn.refstat.get(k, 0) >= 0) { // Ignore artificial classes
        result.XN.set(0, k, result.TN.get((int) sn.refstat.get(k, 0), k));
        result.CN.set(0, k, sn.njobs.get(0, k) / result.XN.get(0, k));
      }
    }

    ((SolverFluidResult) this.result).odeStateVec = analyzer.getXVecIt();
    // In LINE, sn is stored within result too - feels unnecessary, so I've not added
  }

  public SolverResult runMethodSpecificAnalyzerViaLINE() {
    runMethodSpecificAnalyzer();
    return this.result;
  }

  // Run the solver
  @Override
  public void runAnalyzer() {

    long startTime = System.nanoTime();
    boolean hasOpenClasses = false;
    for (NodeType nodetype : sn.nodetypes) {
      if (nodetype == NodeType.Source) {
        hasOpenClasses = true;
        break;
      }
    }

    switch (options.method) {
      case "matrix":
        if (hasOpenClasses && (options.verbose != VerboseLevel.SILENT)) {
          System.err.println(
              "The matrix solver does not support open arrivals. Using options.method = 'closing' instead.");
          options.method = "closing";
        }
        break;
      case "default":
        if (hasOpenClasses) {
          options.method = "closing";
        } else {
          options.method = "matrix";
        }
        break;
      case "closing":
      case "statedep":
        break;
      default:
        System.err.println(
            "SolverFluid does not support the specified method. Using options.method = 'default'.");
        options.method = "default";
    }
    result.method = options.method;

    if (isInfinite(options.timespan[0])) {
      if (options.verbose == VerboseLevel.DEBUG) {
        System.err.println(
            "SolverFluid requires options.timespan[0] to be finite. Setting it to 0.");
      }
      options.timespan[0] = 0;
    }
    if (options.timespan[0] == options.timespan[1]) {
      System.err.println(
          "SolverFluid does not support a timespan that is a single point. Setting options.timespan[0] to 0.");
      options.timespan[0] = 0;
    }
    // TODO: implement the below method and then uncomment
    /*if (enableChecks && !supports(model)) {
      throw new RuntimeException("This model contains features not supported by the solver.");
    }*/

    int M = sn.nstations;
    int K = sn.nclasses;

    Matrix Q = new Matrix(M, K);
    Matrix U = Q.clone();
    Matrix R = Q.clone();
    Matrix T = Q.clone();
    Matrix C = new Matrix(1, K);
    Matrix X = C.clone();
    Matrix[][] Qt = new Matrix[M][K];
    Matrix[][] Ut = new Matrix[M][K];
    Matrix[][] Tt = new Matrix[M][K];
    for (int i = 0; i < M; i++) {
      for (int k = 0; k < K; k++) {
        Qt[i][k] = new Matrix(0, 0);
        Ut[i][k] = new Matrix(0, 0);
        Tt[i][k] = new Matrix(0, 0);
      }
    }

    Matrix s0_sz = new Matrix(1, sn.state.size());
    Matrix s0_id = s0_sz.clone();
    int i = 0;
    for (Station station : sn.stations) {
      s0_sz.set(0, i, sn.state.get(station).getNumRows());
      s0_id.set(0, i, s0_sz.get(0, i) - 1);
      i++;
    }
    s0_id = PopulationLattice.pprod(s0_id);

    while (s0_id.elementMin() >= 0) { // For all possible initial states
      double s0prior_val = 1;
      for (int ind = 0; ind < sn.nnodes; ind++) {
        if (sn.isstateful.get(ind, 0) == 1) {
          int isf = (int) sn.nodeToStateful.get(0, ind);
          // Update prior
          s0prior_val *= sn.stateprior.get(sn.stations.get(isf)).get(0, (int) s0_id.get(0, isf));
          Matrix newState =
              Matrix.extractRows(
                  sn.state.get(sn.stations.get(isf)),
                  (int) s0_id.get(0, isf),
                  (int) s0_id.get(0, isf) + 1,
                  null);
          this.model.getStations().get(isf).setState(newState);
        }
      }
      this.sn = this.model.getStruct(true);

      if (s0prior_val > 0) {
        runMethodSpecificAnalyzer();

        // Note: in LINE, only unique time-step values (and their associated metrics) are stored.
        // The time inefficiency in determining the unique values is, I believe, worse than the
        // space inefficiency in storing larger arrays than is necessary. For that reason I've not
        // transferred the 'unique' functionality across to JLINE

        if (((SolverFluidResult) this.result).odeStateVec.isEmpty()) { // If the solution has failed
          for (int k = 0; k < K; k++) {
            for (int j = 0; j < M; j++) {
              Q.set(j, k, NaN);
              U.set(j, k, NaN);
              R.set(j, k, NaN);
              T.set(j, k, NaN);
            }
            C.set(1, k, NaN);
            X.set(1, k, NaN);
          }

          Matrix nanMatrix = new Matrix(1, 1);
          nanMatrix.set(0, 0, NaN);
          for (int ist = 0; ist < M; ist++) {
            for (int r = 0; r < K; r++) {
              Qt[ist][r] = nanMatrix.clone();
              Ut[ist][r] = nanMatrix.clone();
              Tt[ist][r] = nanMatrix.clone();
            }
          }
        } else {
          Q = Q.add(s0prior_val, result.QN);
          U = U.add(s0prior_val, result.UN);
          R = R.add(s0prior_val, result.RN);
          T = T.add(s0prior_val, result.TN);
          C = C.add(s0prior_val, result.CN);
          X = X.add(s0prior_val, result.XN);
          // Check if this is the first time adding results or not
          if (Qt[0][0].isEmpty()) {
            for (int ist = 0; ist < M; ist++) {
              for (int r = 0; r < K; r++) {
                result.QNt[ist][r].scale(s0prior_val, Qt[ist][r]);
                result.UNt[ist][r].scale(s0prior_val, Ut[ist][r]);
                result.TNt[ist][r].scale(s0prior_val, Tt[ist][r]);
              }
            }
          } else {
            for (int ist = 0; ist < M; ist++) {
              for (int r = 0; r < K; r++) {
                // TODO: implement lines 124 to 136
              }
            }
          }
        }
      }

      // Update initial state
      Matrix pprodInput = s0_sz.clone();
      for (int j = 0; j < s0_sz.length(); j++) {
        pprodInput.set(0, j, pprodInput.get(0, j) - 1);
      }
      s0_id = PopulationLattice.pprod(s0_id, pprodInput);
    }

    result.runtime = (System.nanoTime() - startTime) / 1000000000.0;
    // TODO: replace with setAvgResults method once implemented
    result.QN = Q;
    result.UN = U;
    result.RN = R;
    result.TN = T;
    result.CN = C;
    result.XN = X;
    // TODO: replace with setTranAvgResults method once implemented
    result.QNt = Qt;
    result.UNt = Ut;
    result.TNt = Tt;
  }

  public void getCdfRespT() {
    long startTime = System.nanoTime();
    this.sn = model.getStruct(true);
    this.getAvg(); // Get steady-state solution
    this.options.init_sol = ((SolverFluidResult) this.result).odeStateVec;
    ((SolverFluidResult) this.result).distribC = passageTime();
    ((SolverFluidResult) this.result).distribRuntime = (System.nanoTime() - startTime) / 1000000000.0;
    // TODO: add in setDistribResults method once implemented
  }

  public void getTranCdfPassT() {
    long startTime = System.nanoTime();
    this.sn = model.getStruct(true);
    for (int ind = 0; ind < sn.nnodes; ind++) {
      if (sn.isstateful.get(ind, 0) == 1) {
        int isf = (int) sn.nodeToStateful.get(ind, 0);
        Matrix statePrior = sn.stateprior.get(sn.stations.get(isf));
        if (statePrior.elementSum() != statePrior.elementMax()) {
          throw new RuntimeException(
              "getTranCdfPassT: multiple initial states have non-zero prior - unsupported.");
        }
        // Assign initial state to network
        sn.state.put(
            sn.stations.get(isf),
            Matrix.extractRows(sn.state.get(sn.stations.get(isf)), 0, 1, null));
      }
    }

    initSol();
    ((SolverFluidResult) this.result).distribC = passageTime();
    ((SolverFluidResult) this.result).distribRuntime = (System.nanoTime() - startTime) / 1000000000.0;
    // TODO: add in setDistribResults method once implemented
  }

  public void getProbAggr(int ist) {

    if (ist > sn.nstations) {
      throw new RuntimeException("Station number exceeds the number of stations in the model.");
    }

    if (!this.hasAvgResults()) {
      this.getAvg();
    }

    boolean allAreFinite = true;
    for (int i = 0; i < sn.njobs.getNumCols(); i++) {
      if (isInfinite(sn.njobs.get(0, i))) {
        allAreFinite = false;
        break;
      }
    }

    if (allAreFinite) {
      State.StateMarginalStatistics stats =
          State.toMarginal(
              this.sn,
              ist,
              sn.state.get(sn.stations.get((int) sn.stationToStateful.get(0, ist))),
              null,
              null,
              null,
              null,
              null);
      // Binomial approximation with mean fitted to queue-lengths.
      // Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997.
      ((SolverFluidResult) this.result).logPnir = 0;
      for (int r = 0; r < stats.nir.getNumCols(); r++) {
        // TODO x 3: logPnir = logPnir + nchoosekln(N(r),nir(r));
        // logPnir = logPnir + nir(r)*log(Q(ist,r)/N(r));
        // logPnir = logPnir + (N(r)-nir(r))*log(1-Q(ist,r)/N(r));
      }
      // TODO: return real part of the array only
      ((SolverFluidResult) this.result).Pnir = Math.exp(((SolverFluidResult) this.result).logPnir);
    } else {
      throw new RuntimeException("getProbAggr not yet implemented for models with open classes.");
    }
  }

  // TODO: passage time analysis not tested yet
  private Matrix[][] passageTime() {

    int M = sn.nstations; // Number of Stations
    int K = sn.nclasses; // Number of Classes
    double N = sn.nclosedjobs; // Population
    Matrix S = sn.nservers.clone();
    for (int i = 0; i < M; i++) {
      // Set number of servers in delay station = population
      if (Double.isInfinite(S.get(i, 0))) {
        S.set(i, 0, N);
      }
    }
    Matrix[][][] tmpRT = new Matrix[M][K][2];

    // Initialisation
    Matrix slowrate = new Matrix(M, K);
    for (int i = 0; i < M; i++) {
      for (int k = 0; k < K; k++) {
        // Service completion (exit) rates in each phase
        slowrate.set(
            i,
            k,
            Math.min(
                POSITIVE_INFINITY,
                sn.mu.get(sn.stations.get(i)).get(sn.jobclasses.get(k)).elementMin()));
      }
    }

    // Response time analysis - starting from fixed point found
    JobClass KcClass = null;
    for (int i = 0; i < M; i++) {
      if (sn.nodetypes.get((int) sn.stationToNode.get(i)) != NodeType.Source) {
        for (int k = 0; k < sn.nchains; k++) { // once for each chain
          LinkedList<Integer> idxClassesInChain = new LinkedList<>();
          for (int c = 0; c < sn.chains.getNumCols(); c++) {
            if (sn.chains.get(k, c) == 1) {
              idxClassesInChain.add(c);
            }
          }
          for (Integer c : idxClassesInChain) {
            if (sn.phases.get(i, c) > 0) {

              // Generate ODEs Passage Time
              Matrix phases_c;
              int Kc = K + 1; // Add a single new class
              // Single new class is arbitrarily open
              if (Kc > sn.jobclasses.size()) {
                KcClass = new JobClass(JobClassType.Open, "singleNewClass");
                sn.jobclasses.add(KcClass);
              }
              // Indices of the transient class corresponding to each class in the original model
              // for
              // chain k
              Matrix idxTranCl = new Matrix(1, K);
              for (Integer c2 : idxClassesInChain) {
                idxTranCl.set(0, c2, K);
              }

              Matrix newRT = new Matrix(M * Kc, M * Kc); // New routing table
              Map<Station, Map<JobClass, Map<Integer, Matrix>>> newProc = new HashMap<>();
              Map<Station, Map<JobClass, Matrix>> newMu = new HashMap<>();
              Map<Station, Map<JobClass, Matrix>> newPi = new HashMap<>();

              for (int j = 0; j < M; j++) {
                Station station = sn.stations.get(j);
                newMu.put(station, new HashMap<>());
                newPi.put(station, new HashMap<>());
                newProc.put(station, new HashMap<>());
                for (int r = 0; r < K; r++) {
                  JobClass jobClass = sn.jobclasses.get(r);
                  newProc.get(station).put(jobClass, new HashMap<>());
                  int procSize = sn.proc.get(station).get(jobClass).size();
                  for (int q = 0; q < procSize; q++) {
                    newProc
                        .get(station)
                        .get(jobClass)
                        .put(q, sn.proc.get(station).get(jobClass).get(q).clone());
                  }
                  // Service rates
                  newMu.get(station).put(jobClass, sn.mu.get(station).get(jobClass).clone());
                  // Completion Probabilities
                  newPi.get(station).put(jobClass, sn.phi.get(station).get(jobClass).clone());
                }
                newMu
                    .get(station)
                    .put(KcClass, sn.mu.get(station).get(sn.jobclasses.get(c)).clone());
                newPi
                    .get(station)
                    .put(KcClass, sn.phi.get(station).get(sn.jobclasses.get(c)).clone());

                // PHD Distribution
                for (int r = 0; r < sn.nchains; r++) {
                  JobClass jobClass = sn.jobclasses.get(r);
                  int procSize = sn.proc.get(station).get(jobClass).size();
                  for (int q = 0; q < procSize; q++) {
                    newProc
                        .get(station)
                        .get(jobClass)
                        .put(q, sn.proc.get(station).get(jobClass).get(q).clone());
                  }
                }
                newProc.get(station).put(KcClass, new HashMap<>());
                int procSize = sn.proc.get(station).get(sn.jobclasses.get(c)).size();
                for (int q = 0; q < procSize; q++) {
                  newProc
                      .get(station)
                      .get(KcClass)
                      .put(q, sn.proc.get(station).get(sn.jobclasses.get(c)).get(q).clone());
                }
              }

              // Routing/switching probabilities among basic classes
              int skipRows = 0;
              for (int row = 0; row < M * K; row++) {
                if ((row + skipRows) % Kc == Kc - 1) {
                  skipRows++;
                }
                int skipCols = 0;
                for (int col = 0; col < M * K; col++) {
                  if ((col + skipCols) % Kc == Kc - 1) {
                    skipCols++;
                  }
                  newRT.set(row + skipRows, col + skipCols, sn.rt.get(row, col));
                }
              }

              // Copy routing table from the original to the transient classes (forward)
              double tmpSum = 0;
              for (int row = c; row < M * K; row += K) {
                for (int col = idxClassesInChain.get(0); col < M * K; col += K) {
                  tmpSum += sn.rt.get(row, col);
                }
              }
              if (tmpSum > 0) {
                int rowIter = 0;
                for (int row = c; row < M * K; row += K) {
                  int colIter = 0;
                  for (int col = idxClassesInChain.get(0); col < M * K; col += K) {
                    newRT.set(K + rowIter, K + colIter, sn.rt.get(row, col));
                    colIter += Kc;
                  }
                  rowIter += Kc;
                }
              }

              // Phases of transient classes
              phases_c = sn.phases.clone();
              phases_c.expandMatrix(M, Kc, sn.phases.getNumElements() + M);
              for (int row = 0; row < M; row++) {
                phases_c.set(row, K, sn.phases.get(row, c));
              }

              // Routing matrix from a transient class that completes is diverted back into the
              // original classes
              for (Integer l : idxClassesInChain) { // For each completing class
                for (int j = 0; j < M; j++) {
                  // Return fluid to original class
                  newRT.set(
                      i * Kc + (int) idxTranCl.get(0, c),
                      j * Kc + l,
                      sn.rt.get(i * K + c, j * K + l));
                  // Delete corresponding transition among transient classes
                  newRT.set(
                      i * Kc + (int) idxTranCl.get(0, c), j * Kc + (int) idxTranCl.get(0, l), 0);
                }
              }

              // Setup Initial Point and Empty Array for Solution
              int stateLength = (int) phases_c.elementSum();
              double[] initialState = new double[stateLength];
              double[] nextState = new double[stateLength];
              for (int idx = 0; idx < stateLength; idx++) {
                initialState[idx] = 0;
                nextState[idx] = 0;
              }
              double fluid_c = 0;

              for (int j = 0; j < M; j++) {
                for (int l = 0; l < K; l++) {
                  int idxNew_jl =
                      (int) phases_c.sumSubMatrix(0, j, 0, phases_c.getNumCols())
                          + (int) phases_c.sumSubMatrix(j, j + 1, 0, l);
                  int idxNew_jt =
                      (int) phases_c.sumSubMatrix(0, j, 0, phases_c.getNumCols())
                          + (int) phases_c.sumSubMatrix(j, j + 1, 0, (int) idxTranCl.get(0, l));
                  int idx_jl =
                      (int) sn.phases.sumSubMatrix(0, j, 0, sn.phases.getNumCols())
                          + (int) sn.phases.sumSubMatrix(j, j + 1, 0, l);
                  if (i == j && l == c) {
                    initialState[idxNew_jt] = // mass in phases all moved back into phase 1
                        options.init_sol.sumSubMatrix(
                            0, 1, idx_jl, idx_jl + (int) sn.phases.get(j, l));
                    fluid_c +=
                        options.init_sol.sumSubMatrix(
                            0, 1, idx_jl, idx_jl + (int) sn.phases.get(j, l));
                  } else { // Leave mass as it is
                    int idx = idxNew_jl;
                    for (int q = idx_jl; q < idx_jl + sn.phases.get(j, l); q++) {
                      initialState[idx] = options.init_sol.get(0, q);
                      idx++;
                    }
                  }
                }
              }

              // Determine max integration time
              double minNonZeroRate = POSITIVE_INFINITY;
              for (int row = 0; row < M; row++) {
                for (int col = 0; col < K; col++) {
                  double val = slowrate.get(row, col);
                  if (val > GlobalConstants.CoarseTol && val < minNonZeroRate) {
                    minNonZeroRate = val;
                  }
                }
              }
              // Solve ODE until T = 100 events with slowest exit rate
              double[] tRange = {0, abs(100 / minNonZeroRate)};

              // Indices of new classes at Station i
              LinkedList<Integer> idxN = new LinkedList<>();
              double end = phases_c.sumSubMatrix(i, i + 1, K, Kc);
              for (int idx = 0; idx < end; idx++) {
                idxN.add(
                    idx
                        + (int) phases_c.sumSubMatrix(0, i, 0, phases_c.getNumCols())
                        + (int) phases_c.sumSubMatrix(i, i + 1, 0, K));
              }

              // Set-up the ODEs for the new QN
              FirstOrderDifferentialEquations ode =
                  new ClosingAndStateDepMethodsODE(
                      sn, newMu, newPi, newProc, newRT, S, options, initialState.length);

              // ODE analysis
              Matrix tFull = new Matrix(0, 0);
              Matrix stateFull = new Matrix(0, 0);
              int iter = 1;
              boolean finished = false;
              int tref = 0;
              while (iter <= options.iter_max && !finished) {

                // Solve ODE - y_mean_iter is the transient solution in stage e
                FirstOrderIntegrator odeSolver;
                if (options.stiff) {
                  odeSolver = options.odesolvers.accurateStiffODESolver;
                } else {
                  odeSolver = options.odesolvers.accurateODESolver;
                }
                if (options.tol > 0.001) {
                  odeSolver = options.odesolvers.fastStiffODESolver;
                }

                odeSolver.clearStepHandlers();
                TransientDataHandler stepHandler = new TransientDataHandler(initialState.length);
                odeSolver.addStepHandler(stepHandler);

                try {
                  odeSolver.integrate(ode, tRange[0], initialState, tRange[1], nextState);
                } catch (RuntimeException e) {
                  throw new RuntimeException("ODE Solver Failed.");
                }

                Matrix tmpTIter = stepHandler.tVec;
                Matrix tmpStateIter = stepHandler.xVec;

                iter++;

                if (tFull.isEmpty()) {
                  tFull = tmpTIter.clone();
                  stateFull = tmpStateIter.clone();
                } else {
                  int startRow = tFull.getNumRows();
                  tFull.expandMatrix(
                      tFull.getNumRows() + tmpTIter.getNumRows(),
                      1,
                      tFull.getNumElements() + tmpTIter.getNumElements());
                  int endRow = tFull.getNumRows();
                  for (int row = startRow; row < endRow; row++) {
                    tFull.set(row, 0, tmpTIter.get(row - startRow, 0) + tref);
                  }
                  stateFull = Matrix.concatRows(stateFull, tmpStateIter, null);
                }

                tmpSum = 0;
                int tmpTMax = tmpStateIter.getNumRows();
                for (Integer idx : idxN) {
                  tmpSum += tmpStateIter.get(tmpTMax - 1, idx);
                }
                if (tmpSum < 0.000000001) {
                  finished = true;
                }
                tref += tmpTIter.get(tmpTMax - 1, 0);
                for (int idx = 0; idx < initialState.length; idx++) {
                  initialState[idx] = tmpStateIter.get(tmpTMax - 1, idx);
                }
              }

              // Retrieve response time CDF for class k
              int fullTMax = tFull.getNumRows();
              tmpRT[i][c][0] = tFull;
              tmpRT[i][c][1] = tFull.clone();
              if (fluid_c > 0) {
                for (int row = 0; row < fullTMax; row++) {
                  tmpSum = 0;
                  for (Integer idx : idxN) {
                    tmpSum += stateFull.get(row, idx);
                  }
                  tmpRT[i][c][1].set(row, 0, 1 - tmpSum / fluid_c);
                }
              } else {
                tmpRT[i][c][1].ones();
              }

              if (iter > options.iter_max) {
                System.err.format(
                    "Maximum number of iterations reached when computing the response time distribution. "
                        + "Response time distributions may be inaccurate. Increase option.iter_max (currently at %d).",
                    options.iter_max);
              }
            }
          }
        }
      }
    }

    if (KcClass != null) {
      sn.jobclasses.remove(KcClass);
    }

    Matrix[][] RTret = new Matrix[M][K];
    for (int i = 0; i < M; i++) {
      for (int c = 0; c < K; c++) {
        RTret[i][c] = Matrix.concatColumns(tmpRT[i][c][1], tmpRT[i][c][0], null);
      }
    }
    return RTret;
  }

  public static void getFeatureSet() {
    // TODO: implementation - note return type should likely not be void
    throw new RuntimeException("getFeatureSet() has not yet been implemented in JLINE.");
  }

  public static boolean supports(Model model) {
    // TODO: implementation
    throw new RuntimeException("supports() has not yet been implemented in JLINE.");
  }

  public static SolverOptions defaultOptions() {
    return new SolverOptions(SolverType.FLUID);
  }
}
