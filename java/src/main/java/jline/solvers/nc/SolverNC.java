package jline.solvers.nc;

import jline.api.PFQN;
import jline.api.SN;
import jline.util.Maths;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.GlobalConstants;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Node;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.state.State;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.util.Matrix;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static jline.lib.KPCToolbox.*;
import static jline.api.PFQN.pfqn_nc;
import static jline.api.SN.snDeaggregateChainResults;

import static jline.api.SN.snGetDemandsChain;
import static jline.lang.state.State.toMarginal;


public class SolverNC extends NetworkSolver {

  public SolverNC(Network model, SolverOptions options) {
    super(model, "SolverNC", options);
    this.sn = model.getStruct(false);
    this.result = new SolverNCResult();
  }

  public SolverNC(Network model) {
    super(model, "SolverNC");
    this.sn = model.getStruct(false);
    this.result = new SolverNCResult();
  }

  public NetworkStruct getStruct() {
    if (this.sn == null)
      this.sn = this.model.getStruct(false);
    return this.sn;
  }

  public void setStruct(NetworkStruct sn) {
    this.sn = sn;
  }

  public static SolverOptions defaultOptions() {
    return new SolverOptions(SolverType.NC);
  }

  @Override
  public void runAnalyzer() throws IllegalAccessException {
    if (this.model == null)
      throw new RuntimeException("Model is not provided");
    if (this.sn == null)
      this.sn = this.model.getStruct(false);
    if (this.options == null)
      this.options = new SolverOptions(SolverType.NC);

    //NCRunner runner = new NCRunner(this);
    //this.result = runner.run();
  }


  public Double getProb(Node node, Matrix state) {
    if (GlobalConstants.DummyMode) {
      return Double.NaN;
    }

    Matrix state_new = state.clone();

    if (state_new == null || state_new.isEmpty()) {
      state_new = sn.state.get((StatefulNode) sn.nodes.get((int) sn.nodeToStateful.get(node.getNodeIdx())));
    }

    long startTimeMillis = System.currentTimeMillis();
    NetworkStruct sn = getStruct();
    int ist = (int) sn.nodeToStation.get(node.getNodeIdx());
    sn.state.put((StatefulNode) sn.nodes.get(ist), state_new);
    resetRandomGeneratorSeed(options.seed);

    Matrix Pnir;

    SolverNCResult ncResult = (SolverNCResult) this.result;
    if (ncResult != null && ncResult.prob != null
            && ncResult.prob.logNormConstAggr != null && Double.isFinite(ncResult.prob.logNormConstAggr)) {
      Pnir = solver_nc_marg(sn, this.options, ncResult.prob.logNormConstAggr).lPr;
    } else {
      SolverNCMargReturn ret = solver_nc_marg(sn, this.options, null);
      Pnir = ret.lPr;
      ((SolverNCResult) this.result).prob.logNormConstAggr = ret.G;
    }
    ((SolverNCResult) this.result).solver = this.name;
    ((SolverNCResult) this.result).prob.marginal = Pnir;
    long endTimeMillis = System.currentTimeMillis();
    double runtime = (endTimeMillis-startTimeMillis) / 1000.0;
    this.result.runtime = runtime;

    return Pnir.get(ist);
  }

  public static SolverNCMargReturn solver_nc_marg(NetworkStruct sn, SolverOptions options, Double lG) {
    int M = sn.nstations;
    int K = sn.nclasses;
    Map<StatefulNode, Matrix> state = sn.state;
    Matrix S = sn.nservers;
    Matrix V = new Matrix(sn.nstateful, K);
    for (int i = 0; i < sn.visits.size(); i++) {
      V = V.add(1, sn.visits.get(i));
    }
    //V = cellsum(sn.visits);
    Matrix rates = sn.rates;
    Matrix ST = rates.clone();
    for (int i = 0; i < ST.numRows; i++) {
      for (int j = 0; j < ST.numCols; j++) {
        ST.set(i, j, 1.0/ST.get(i, j));
      }
    }
    ST.removeNaN();

    SN.snGetDemandsChainReturn ret = snGetDemandsChain(sn);
    Matrix Lchain = ret.Lchain;
    Matrix STchain = ret.STchain;
    Matrix Nchain = ret.Nchain;

    long startTimeMillis = System.currentTimeMillis();
    M = STchain.numRows;
    K = STchain.numCols;

    Matrix mu = new Matrix(0, (int) Nchain.elementSum());
    for (int i = 0; i < M; i++) {
      Matrix tmp = new Matrix(1, (int) Nchain.elementSum());
      if (Double.isInfinite(S.get(i))) {
        for (int j = 0; j < tmp.length(); j++) {
          tmp.set(j, j+1);
        }
      } else {
        for (int j = 0; j < tmp.length(); j++) {
          tmp.set(j, Math.min(j+1, S.get(i)));
        }
      }
      mu = Matrix.concatRows(mu, tmp, null);
    }


    if (lG == null) {
      Matrix Z_tmp = Nchain.clone();
      Z_tmp.fill(0.0);
      lG = PFQN.pfqn_ncld(Lchain, Nchain, Z_tmp, mu, options).lG;
    }

    double G = Math.exp(lG);
    Matrix lPr = new Matrix(sn.nstations, 1);
    lPr.fill(0.0);

    for (int ist = 0; ist < sn.nstations; ist++) {
      int ind = (int) sn.stationToNode.get(ist);
      int isf = (int) sn.stationToStateful.get(ist);
      State.StateMarginalStatistics ret1 = toMarginal(sn, ind, state.get((StatefulNode) sn.nodes.get(isf)),
              null, null, null, null, null);
      Matrix nirvec = ret1.nir;
      Matrix sivec = ret1.sir;
      List<Matrix> kirvec = ret1.kir;

      if (nirvec.elementMin() < -1e-6) {
        lPr.set(M-1, Double.NaN);
      } else {
        //nivec_chain = nirvec * sn.chains'; Assume it to be matrix
        Matrix nivec_chain = nirvec.mult(sn.chains.transpose());

        Matrix Lchain_tmp = new Matrix(0, Lchain.numCols);
        Matrix mu_tmp = new Matrix(0, mu.numCols);
        for (int i = 0; i < sn.nstations; i++) {
          if (i != ist) {
            Matrix Lchain_row_i = Matrix.extractRows(Lchain, i, i+1, null);
            Matrix mu_row_i = Matrix.extractRows(mu, i, i+1, null);
            Lchain_tmp = Matrix.concatRows(Lchain_tmp, Lchain_row_i, null);
            mu_tmp = Matrix.concatRows(mu_tmp, mu_row_i, null);
          }
        }
        Matrix Nchain_tmp = Nchain.clone();
        for (int i = 0; i < Nchain_tmp.length(); i++) {
          Nchain_tmp.set(i, Nchain_tmp.get(i)-nivec_chain.get(i));
        }
        Matrix Zchain_tmp = Nchain.clone();
        Zchain_tmp.fill(0.0);

        double lG_minus_i = PFQN.pfqn_ncld(Lchain_tmp, Nchain_tmp, Zchain_tmp, mu_tmp, options).lG;
        double lF_i = 0.0;

        switch (sn.sched.get(sn.stations.get(ist))) {
          case FCFS:
            for (int r = 0; r < K; r++) {
              Map<Integer, Matrix> PHr = sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r));
              if (!PHr.isEmpty()) {
                Matrix kir = new Matrix(1, kirvec.size());
                for (int i = 0; i < kir.length(); i++) {
                  kir.set(i, kirvec.get(i).get(0, r));
                }
                if (kir.length() > 1) {
                  throw new RuntimeException("solver_nc_marg: Cannot return state probability " +
                          "because the product-form solution requires exponential service times at FCFS nodes.");
                }

                if (Math.abs(ST.get(ist, r)-Matrix.extractRows(ST,ist,ist+1,null).elementMax())>1e-6) {
                  throw new RuntimeException("solver_nc_marg: Cannot return state probability " +
                          "because the product-form solution requires identical service times at FCFS nodes.");
                }
              }
            }
            boolean allZero = true;
            for (int i = 0; allZero && i < sivec.numRows; i++) {
              for (int j = 0; allZero && j < sivec.numCols; j++) {
                if (Math.abs(sivec.get(i, j)) > 1e-6) {
                  allZero = false;
                }
              }
            }
            if (!allZero) {
              Matrix tmp = Matrix.extractRows(nirvec, 0, 1, null);
              for (int i = 0; i < tmp.length(); i++) {
                tmp.set(i, tmp.get(i)*Math.log(V.get(ist,K-1)));
              }
              double sum_kirvec = 0.0;
              for (int i = 0; i < kirvec.size(); i++) {
                sum_kirvec += kirvec.get(i).elementSum();
              }
              Matrix mu_row_ist = new Matrix(1, (int) sum_kirvec);
              Matrix.extract(mu, ist, ist+1, 0, (int) sum_kirvec, mu_row_ist, 0, 0);
              for (int i = 0; i < mu_row_ist.length(); i++) {
                mu_row_ist.set(i, Math.log(mu_row_ist.get(i)));
              }
              lF_i += (tmp.elementSum() - mu_row_ist.elementSum());
            } else {
              lF_i = 0.0;
            }
            break;
            /* Need pfqn_ncldld
          case SIRO: {
            break;
          }
             */
          case PS:
            for (int r = 0; r < K; r++) {
              Map<Integer, Matrix> PHr = sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r));
              if (!PHr.isEmpty()) {
                Matrix kir = new Matrix(1, kirvec.size());
                for (int i = 0; i < kir.length(); i++) {
                  kir.set(i, kirvec.get(i).get(0, r));
                }

                Matrix PHr_tmp = PHr.get(0);
                for (int i = 0; i < PHr_tmp.numRows; i++) {
                  for (int j = 0; j < PHr_tmp.numCols; j++) {
                    PHr_tmp.set(i, j, -1.0*PHr_tmp.get(i, j));
                  }
                }
                PHr_tmp = PHr_tmp.inv();
                Matrix Ar = map_pie(PHr.get(0), PHr.get(1)).mult(PHr_tmp);

                Matrix kir_tmp = Ar.clone();
                for (int i = 0; i < kir_tmp.numRows; i++) {
                  for (int j = 0; j < kir_tmp.numCols; j++) {
                    kir_tmp.set(i, j, kir.get(i, j) * Math.log(V.get(ist, r)*kir_tmp.get(i, j)));
                  }
                }

                lF_i += (kir_tmp.elementSum() - Matrix.factln(kir).elementSum());
              }
            }

            double sum_kirvec = 0.0;
            for (int i = 0; i < kirvec.size(); i++) {
              sum_kirvec += kirvec.get(i).elementSum();
            }
            Matrix mu_row_ist = new Matrix(1, (int) sum_kirvec);
            Matrix.extract(mu, ist, ist+1, 0, (int) sum_kirvec, mu_row_ist, 0, 0);
            for (int i = 0; i < mu_row_ist.numCols; i++) {
              mu_row_ist.set(i, Math.log(mu_row_ist.get(i)));
            }
            lF_i += (Maths.factln(sum_kirvec) - mu_row_ist.elementSum());
            break;
          case INF:
            for (int r = 0; r < K; r++) {
              Map<Integer, Matrix> PHr = sn.proc.get(sn.stations.get(ist)).get(sn.jobclasses.get(r));
              if (!PHr.isEmpty()) {
                Matrix kir = new Matrix(1, kirvec.size());
                for (int i = 0; i < kir.length(); i++) {
                  kir.set(i, kirvec.get(i).get(0, r));
                }
                Matrix PHr_tmp = PHr.get(0);
                for (int i = 0; i < PHr_tmp.numRows; i++) {
                  for (int j = 0; j < PHr_tmp.numCols; j++) {
                    PHr_tmp.set(i, j, -1.0*PHr_tmp.get(i, j));
                  }
                }
                PHr_tmp = PHr_tmp.inv();
                Matrix Ar = map_pie(PHr.get(0), PHr.get(1)).mult(PHr_tmp);

                Matrix kir_tmp = Ar.clone();
                for (int i = 0; i < kir_tmp.numRows; i++) {
                  for (int j = 0; j < kir_tmp.numCols; j++) {
                    kir_tmp.set(i, j, kir.get(i, j) * Math.log(V.get(ist, r)*kir_tmp.get(i, j)));
                  }
                }

                lF_i += (kir_tmp.elementSum() - Matrix.factln(kir).elementSum());
              }
            }
            break;
        }
        lPr.set(ist, lF_i + lG_minus_i - lG);
      }
    }
    long endTimeMillis = System.currentTimeMillis();
    double runtime = (endTimeMillis-startTimeMillis) / 1000.0;
    lPr.removeNaN();
    return new SolverNCMargReturn(lPr, G, runtime);
  }

  public static SolverNCJointReturn solver_nc_joint(NetworkStruct sn, SolverOptions options) {
    Map<StatefulNode, Matrix> state = sn.state;
    Matrix S = sn.nservers;
    Matrix rates = sn.rates;
    Matrix ST = rates.clone();
    for (int i = 0; i < ST.numRows; i++) {
      for (int j = 0; j < ST.numCols; j++) {
        ST.set(i, j, 1.0/ST.get(i, j));
      }
    }
    ST.removeNaN();
    SN.snGetDemandsChainReturn ret = snGetDemandsChain(sn);
    Matrix STchain = ret.STchain;
    Matrix Vchain = ret.Vchain;
    Matrix alpha = ret.alpha;
    Matrix Nchain = ret.Nchain;
    long startTimeMillis = System.currentTimeMillis();
    int M = STchain.numRows;
    int K = ST.numCols;
    Matrix Lchain = new Matrix(0, K);

    Matrix mu_chain = new Matrix(0, (int) Nchain.elementSum());
    for (int i = 0; i < M; i++) {
      Matrix tmp = new Matrix(1, (int) Nchain.elementSum());
      if (Double.isInfinite(S.get(i))) {
        for (int j = 0; j < tmp.length(); j++) {
          tmp.set(j, j+1);
        }
      } else {
        for (int j = 0; j < tmp.length(); j++) {
          tmp.set(j, Math.min(j+1, S.get(i)));
        }
      }
      mu_chain = Matrix.concatRows(mu_chain, tmp, null);

      Matrix tmp1 = Matrix.extractRows(STchain, i, i+1, null)
              .elementMult(Matrix.extractRows(Vchain, i, i+1, null), null);
      Lchain = Matrix.concatRows(Lchain, tmp1, null);
    }

    Matrix Z_tmp = Nchain.clone();
    Z_tmp.fill(0.0);
    double lG = PFQN.pfqn_ncld(Lchain, Nchain, Z_tmp, mu_chain, options).lG;
    double lPr = 0.0;

    for (int i = 0; i < M; i++) {
      int isf = (int) sn.stationToStateful.get(i);
      State.StateMarginalStatistics ret1 = toMarginal(sn, i, state.get((StatefulNode) sn.nodes.get(isf)),
              null, null, null, null, null);
      Matrix nivec = ret1.nir;
      Matrix nivec_chain = nivec.mult(sn.chains.transpose());
      Matrix Lchain_i = Matrix.extractRows(Lchain, i, i+1, null);
      Matrix ST_i = Matrix.extractRows(ST, i, i+1, null);
      Matrix alpha_i = Matrix.extractRows(alpha, i, i+1, null);
      Matrix STchain_i = Matrix.extractRows(STchain, i, i+1, null);
      Matrix mu_chain_i = Matrix.extractRows(mu_chain, i, i+1, null);
      Matrix Zvec_chain = nivec_chain.clone();
      Zvec_chain.fill(0.0);
      Matrix Zvec = nivec.clone();
      Zvec.fill(0.0);

      double lF_i = PFQN.pfqn_ncld(Lchain_i, nivec_chain, Zvec_chain, mu_chain_i, options).lG;
      double lg0_i = PFQN.pfqn_ncld(ST_i.elementMult(alpha_i, null), nivec, Zvec, mu_chain_i, options).lG;
      double lG0_i = PFQN.pfqn_ncld(STchain_i, nivec_chain, Zvec_chain, mu_chain_i, options).lG;
      lPr = lPr + lF_i + (lg0_i - lG0_i);
    }
    double Pr = Math.exp(lPr - lG);
    long endTimeMillis = System.currentTimeMillis();
    double runtime = (endTimeMillis-startTimeMillis) / 1000.0;
    double G = Math.exp(lG);
    return new SolverNCJointReturn(Pr, G, lG, runtime);
  }

  public static SolverNCReturn solver_nc(NetworkStruct sn, SolverOptions options) {
    int M = sn.nstations;
    Matrix nservers = sn.nservers;
    Matrix NK = sn.njobs.transpose();
    Map<Station, SchedStrategy> sched = sn.sched;
    int C = sn.nchains;
    Matrix SCV = sn.scv;
    int K = sn.nclasses;

    Matrix V = new Matrix(sn.nstateful, K);
    for (int i = 0; i < sn.visits.size(); i++) {
      V = V.add(1, sn.visits.get(i));
    }
    //V = cellsum(sn.visits);

    Matrix rates = sn.rates;
    Matrix ST = rates.clone();
    for (int i = 0; i < ST.numRows; i++) {
      for (int j = 0; j < ST.numCols; j++) {
        ST.set(i, j, 1.0/ST.get(i, j));
      }
    }
    ST.removeNaN();
    Matrix ST0 = ST.clone();
    Matrix Nchain = new Matrix(1, C);
    Nchain.fill(0.0);

    for (int c = 0; c < C; c++) {
      Matrix inchain = sn.inchain.get(c);
      double sum = 0;
      for (int col = 0; col < sn.inchain.get(c).getNumCols(); col++) {
        sum += NK.get((int) inchain.get(0, col));
      }
      Nchain.set(c, sum);
    }

    List<Integer> openChains = new ArrayList<>();
    List<Integer> closedChains = new ArrayList<>();

    for (int c = 0; c < C; c++) {
      if (Double.isInfinite(Nchain.get(c))) {
        openChains.add(c);
      } else {
        closedChains.add(c);
      }
    }

    Matrix gamma = new Matrix(1, M);
    gamma.fill(0.0);
    Matrix eta_1 = new Matrix(1, M);
    eta_1.fill(0.0);
    Matrix eta = new Matrix(1, M);
    eta.fill(1.0);
    C = sn.nchains;

    if (!sched.containsValue(SchedStrategy.FCFS)) {
      options.iter_max = 1;
    }

    int it = 0;
    Matrix tmp_eta = new Matrix(1, M);
    for (int i = 0; i < M; i++) {
      tmp_eta.set(i, Math.abs(1-eta.get(i)/eta_1.get(i)));
    }

    Matrix lambda = null;
    Matrix Lchain = null;
    Matrix STchain = null;
    Matrix Vchain = null;
    Matrix alpha = null;
    Matrix Q = null;
    Matrix U = null;
    Matrix R = null;
    Matrix T = null;
    Matrix X = null;
    Matrix STeff = null;
    Double lG = null;
    String method = null;

    while (tmp_eta.elementMax() > options.iter_tol && it < options.iter_max) {
      it++;
      eta_1 = eta;

      if (it == 1) {
        lambda = new Matrix(1, C);
        lambda.fill(0.0);
        SN.snGetDemandsChainReturn ret = snGetDemandsChain(sn);
        Lchain = ret.Lchain;
        STchain = ret.STchain;
        Vchain = ret.Vchain;
        alpha = ret.alpha;
        Nchain = ret.Nchain;

        for (int c = 0; c < C; c++) {
          Matrix inchain = sn.inchain.get(c);
          boolean isOpenChain = false;
          for (int col = 0; col < inchain.getNumCols(); col++) {
            if (Double.isInfinite(NK.get((int) inchain.get(0, col)))) {
              isOpenChain = true;
              break;
            }
          }
          for (int i = 0; i < M; i++) {
            if (isOpenChain && Math.abs(i-sn.refstat.get((int) inchain.get(0))) < 1e-6) {
              lambda.set(c, 1 / STchain.get(i, c));
            }
          }
        }
      } else {
        for (int c = 0; c < C; c++) {
          Matrix inchain = sn.inchain.get(c);
          for (int i = 0; i < M; i++) {
            Matrix ST_tmp = new Matrix(1, inchain.getNumCols());
            Matrix alpha_tmp = new Matrix(1, inchain.getNumCols());
            for (int col = 0; col < inchain.getNumCols(); col++) {
              ST_tmp.set(col, ST.get(i, (int) inchain.get(0, col)));
              alpha_tmp.set(col, alpha.get(i, (int) inchain.get(0, col)));
            }
            STchain.set(i, c, ST_tmp.mult(alpha_tmp.transpose()).get(0));
            Lchain.set(i, c, Vchain.get(i, c) * STchain.get(i, c));
          }
        }
      }

      STchain.removeInfinity();
      Lchain.removeInfinity();
      Matrix Lms = new Matrix(M, C);
      Matrix Z = new Matrix(M, C);
      Matrix Zms = new Matrix(M, C);
      Lms.fill(0.0);
      Z.fill(0.0);
      Zms.fill(0.0);

      List<Integer> infServers = new ArrayList<>();
      for (int i = 0; i < M; i++) {
        if (Double.isInfinite(nservers.get(i))) {
          infServers.add(i);
          for (int j = 0; j < C; j++) {
            Z.set(i, j, Lchain.get(i, j));
          }
        } else {
          if (options.method != null && options.method.equalsIgnoreCase("exact") && nservers.get(i) > 1) {
            System.out.println("SolverNC does not support exact multiserver yet. Switching to approximate method.");
          }
          for (int j = 0; j < C; j++) {
            Lms.set(i, j, Lchain.get(i, j) / nservers.get(i));
            Zms.set(i, j, Lchain.get(i, j) * (nservers.get(i)-1) / nservers.get(i));
          }
        }
      }

      Matrix Z_new = new Matrix(1, C);
      for (int i = 0; i < C; i++) {
        Z_new.set(i, Z.sumCols(i) + Zms.sumCols(i));
      }
      Matrix Z_tmp_append_0 = new Matrix(1, Z_new.numCols+1);
      Z_tmp_append_0.set(Z_new.numCols, 0.0);
      Matrix.extract(Z_new, 0, 1, 0, Z_new.length(), Z_tmp_append_0, 0, 0);

      PFQN.pfqnNcXQReturn ret = pfqn_nc(lambda,Lms,Nchain,Z_new, options);
      lG = ret.lG;
      Matrix Xchain = ret.X;
      Matrix Qchain = ret.Q;
      method = ret.method;

      if (Zms.sumCols().elementMin() > GlobalConstants.FineTol) {
        Xchain = new Matrix(0, 0);
        Qchain = new Matrix(0, 0);
      }

      if (Xchain.isEmpty()) {
        Xchain = lambda.clone();
        Qchain = new Matrix(M, C);
        Qchain.fill(0.0);
        for (int r: closedChains) {
          Matrix Nchain_tmp = Nchain.clone();
          Nchain_tmp.set(r, Nchain_tmp.get(r) - 1);
          Matrix Nchain_tmp_append_1 = new Matrix(1, Nchain.numCols+1);
          Nchain_tmp_append_1.set(Nchain.numCols, 1.0);
          Matrix.extract(Nchain_tmp, 0, 1, 0, Nchain_tmp.length(), Nchain_tmp_append_1, 0, 0);
          Xchain.set(r, Math.exp(pfqn_nc(lambda,Lms,Nchain_tmp,Z_new, options).lG - lG));
          for (int i = 0; i < M; i++) {
            if (Lchain.get(i, r) > 1e-6) {
              if (Double.isInfinite(nservers.get(i))) {
                Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r));
              } else {
                Matrix lambda_tmp = new Matrix(1, lambda.length()+1);
                lambda_tmp.fill(0.0);
                Matrix.extract(lambda, 0, 1, 0, lambda.length(), lambda_tmp, 0, 0);

                Matrix L_tmp = new Matrix(0, Lms.numCols+1);
                for (int row = 0; row < Lms.numRows; row++) {
                  if (row != i) {
                    Matrix L_tmp_row = new Matrix(1, Lms.numCols+1);
                    L_tmp_row.set(Lms.numCols, 0.0);
                    Matrix.extract(Lms, row, row+1, 0, Lms.numCols, L_tmp_row, 0, 0);
                    L_tmp = Matrix.concatRows(L_tmp, L_tmp_row, null);
                  }
                }
                Matrix L_tmp_row = new Matrix(1, Lms.numCols+1);
                L_tmp_row.set(Lms.numCols, 1.0);
                Matrix.extract(Lms, i, i+1, 0, Lms.numCols, L_tmp_row, 0, 0);
                L_tmp = Matrix.concatRows(L_tmp, L_tmp_row, null);

                double res = pfqn_nc(lambda_tmp,L_tmp,Nchain_tmp_append_1,Z_tmp_append_0, options).lG;
                Qchain.set(i, r, Zms.get(i, r)*Xchain.get(r)+Lms.get(i,r)*Math.exp(res - lG));
              }
            }
          }
          Qchain.removeNaN();
        }

        for (int r: openChains) {
          for (int i = 0; i < M; i++) {
            Matrix lambda_open = new Matrix(1, openChains.size());
            Matrix Lchain_i_open = new Matrix(1, openChains.size());
            Matrix Qchain_i_closed = new Matrix(1, closedChains.size());
            for (int j = 0; j < openChains.size(); j++) {
              lambda_open.set(j, lambda.get(openChains.get(j)));
              Lchain_i_open.set(j, Lchain.get(i, openChains.get(j)));
            }
            for (int j = 0; j < closedChains.size(); j++) {
              Qchain_i_closed.set(j, Qchain.get(i, closedChains.get(j)));
            }
            Qchain.set(i, r, lambda.get(r)*Lchain.get(i, r)/
                    (1-lambda_open.mult(Lchain_i_open.transpose()).get(0)/nservers.get(i))
                    *(1+Qchain_i_closed.elementSum()));
          }
        }
      } else {
        for (int r = 0; r < C; r++) {
          for (int i = 0; i < M; i++) {
            if (Lchain.get(i, r) > 1e-6) {
              if (Double.isInfinite(nservers.get(i))) {
                Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r));
              }
            }
          }
        }
      }

      boolean allNaN = true;
      for (int i = 0; i < Xchain.length(); i++) {
        if (!Double.isNaN(Xchain.get(i))) {
          allNaN = false;
          break;
        }
      }

      if (allNaN) {
        System.out.println("Normalizing constant computations produced a floating-point range exception. Model is likely too large.");
      }

      //Z = sum(Z(1:M,:),1);

      Matrix Rchain = new Matrix(Qchain.numRows, Qchain.numCols);
      for (int i = 0; i < Rchain.numRows; i++) {
        for (int j = 0; j < Rchain.numCols; j++) {
          Rchain.set(i, j, Qchain.get(i, j)/Xchain.get(j)/Vchain.get(i, j));
        }
      }

      for (int i = 0; i < infServers.size(); i++) {
        int row = infServers.get(i);
        for (int j = 0; j < Rchain.numCols; j++) {
          Rchain.set(row, j, Lchain.get(row, j) / Vchain.get(row, j));
        }
      }

      Matrix Tchain = new Matrix(Vchain.numRows, Vchain.numCols);
      for (int i = 0; i < Tchain.numRows; i++) {
        for (int j = 0; j < Tchain.numCols; j++) {
          Tchain.set(i, j, Xchain.get(j) * Vchain.get(i, j));
        }
      }

      SN.snDeaggregateChainResultsReturn ret1 = snDeaggregateChainResults(sn, Lchain, ST, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain);
      Q = ret1.Q;
      U = ret1.U;
      R = ret1.R;
      T = ret1.T;
      X = ret1.X;
      STeff = ST.clone();
      //TODO: Depend on npfqn_nonexp_approx
      //npfqn_nonexp_approx(options.config.highVar,sn,ST0,V,SCV,T,U,gamma,nservers);

      for (int i = 0; i < M; i++) {
        tmp_eta.set(i, Math.abs(1-eta.get(i)/eta_1.get(i)));
      }
    }
    Q.abs();
    R.abs();
    X.abs();
    U.abs();
    X.removeInfinity();
    U.removeInfinity();
    Q.removeInfinity();
    R.removeInfinity();
    return new SolverNCReturn(Q, U, R, T, C, X, lG, STeff, it, method);
  }

  public static class SolverNCMargReturn {
    public Matrix lPr;
    public double G;
    public double runtime;
    public SolverNCMargReturn(Matrix lPr, double G, double runtime) {
      this.lPr = lPr;
      this.G = G;
      this.runtime = runtime;
    }
  }

  public static class SolverNCJointReturn {
    public double Pr;
    public double G;
    public double lG;
    public double runtime;
    public SolverNCJointReturn(double Pr, double G, double lG, double runtime) {
      this.Pr = Pr;
      this.G = G;
      this.lG = lG;
      this.runtime = runtime;
    }
  }

  public static class SolverNCReturn {
    public Matrix Q;
    public Matrix U;
    public Matrix R;
    public Matrix T;
    public int C;
    public Matrix X;
    public double lG;
    public Matrix STeff;
    public int it;
    public String method;

    public SolverNCReturn(Matrix Q, Matrix U, Matrix R, Matrix T, int C,
                          Matrix X, double lG, Matrix STeff, int it, String method) {
      this.Q = Q;
      this.U = U;
      this.R = R;
      this.T = T;
      this.C = C;
      this.X = X;
      this.lG = lG;
      this.STeff = STeff;
      this.it = it;
      this.method = method;
    }
  }
}
