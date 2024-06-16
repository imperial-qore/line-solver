package jline.solvers.nc;

import jline.api.NPFQN;
import jline.api.PFQN;
import jline.api.SN;
import jline.lang.FeatureSet;
import jline.solvers.mam.SolverMAM;
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

import java.util.*;

import static jline.api.PFQN.*;
import static jline.api.SN.*;
import static jline.lib.KPCToolbox.*;

import static jline.lang.state.State.toMarginal;


public class SolverNC extends NetworkSolver {

  public SolverNC(Network model, SolverOptions options) {
    super(model, "SolverNC", options);
    this.sn = model.getStruct(false);
    this.result = new SolverNCResult();
  }

  public SolverNC(Network model) {
    super(model, "SolverNC", new NCOptions());
    this.sn = model.getStruct(false);
    this.result = new SolverNCResult();
  }

  public SolverNC(Network model, String method) {
    super(model, "SolverNC", new NCOptions().method(method));
    this.sn = model.getStruct(false);
    this.result = new SolverNCResult();
  }

  public SolverNC(Network model, Object... varargin) {
    super(model, "SolverNC", new NCOptions());
    this.options = SolverNC.parseOptions(this.options, varargin);
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
      this.options = new NCOptions();

    NCRunner runner = new NCRunner(this);
    this.result = runner.run();

    int M = sn.nstations;
    int R = sn.nclasses;
    if(!this.result.TN.isEmpty()){
      this.result.AN = new Matrix(M,R,M*R);
      for(int i=0;i<M;i++){
        for(int j=0;j<M;j++){
          for(int k=0;k<R;k++){
            for(int r=0;r<R;r++){
              this.result.AN.set(i,k,this.result.AN.get(i,k)+this.result.TN.get(j,r)*sn.rt.get(j*R+r,i*R+k));
            }
          }
        }
      }
    } else {
      this.result.AN = new Matrix(0,0,0);
    }
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
    for (int i = 0; i < ST.getNumRows(); i++) {
      for (int j = 0; j < ST.getNumCols(); j++) {
        ST.set(i, j, 1.0/ST.get(i, j));
      }
    }
    ST.removeNaN();

    SN.snGetDemandsChainReturn ret = snGetDemandsChain(sn);
    Matrix Lchain = ret.Lchain;
    Matrix STchain = ret.STchain;
    Matrix Nchain = ret.Nchain;

    long startTimeMillis = System.currentTimeMillis();
    M = STchain.getNumRows();
    K = STchain.getNumCols();

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

        Matrix Lchain_tmp = new Matrix(0, Lchain.getNumCols());
        Matrix mu_tmp = new Matrix(0, mu.getNumCols());
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
            for (int i = 0; allZero && i < sivec.getNumRows(); i++) {
              for (int j = 0; allZero && j < sivec.getNumCols(); j++) {
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
                for (int i = 0; i < PHr_tmp.getNumRows(); i++) {
                  for (int j = 0; j < PHr_tmp.getNumCols(); j++) {
                    PHr_tmp.set(i, j, -1.0*PHr_tmp.get(i, j));
                  }
                }
                PHr_tmp = PHr_tmp.inv();
                Matrix Ar = map_pie(PHr.get(0), PHr.get(1)).mult(PHr_tmp);

                Matrix kir_tmp = Ar.clone();
                for (int i = 0; i < kir_tmp.getNumRows(); i++) {
                  for (int j = 0; j < kir_tmp.getNumCols(); j++) {
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
            for (int i = 0; i < mu_row_ist.getNumCols(); i++) {
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
                for (int i = 0; i < PHr_tmp.getNumRows(); i++) {
                  for (int j = 0; j < PHr_tmp.getNumCols(); j++) {
                    PHr_tmp.set(i, j, -1.0*PHr_tmp.get(i, j));
                  }
                }
                PHr_tmp = PHr_tmp.inv();
                Matrix Ar = map_pie(PHr.get(0), PHr.get(1)).mult(PHr_tmp);

                Matrix kir_tmp = Ar.clone();
                for (int i = 0; i < kir_tmp.getNumRows(); i++) {
                  for (int j = 0; j < kir_tmp.getNumCols(); j++) {
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
    for (int i = 0; i < ST.getNumRows(); i++) {
      for (int j = 0; j < ST.getNumCols(); j++) {
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
    int M = STchain.getNumRows();
    int K = ST.getNumCols();
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
    for (int i = 0; i < ST.getNumRows(); i++) {
      for (int j = 0; j < ST.getNumCols(); j++) {
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
      Matrix Z_tmp_append_0 = new Matrix(1, Z_new.getNumCols()+1);
      Z_tmp_append_0.set(Z_new.getNumCols(), 0.0);
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
          Matrix Nchain_tmp_append_1 = new Matrix(1, Nchain.getNumCols()+1);
          Nchain_tmp_append_1.set(Nchain.getNumCols(), 1.0);
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

                Matrix L_tmp = new Matrix(0, Lms.getNumCols()+1);
                for (int row = 0; row < Lms.getNumRows(); row++) {
                  if (row != i) {
                    Matrix L_tmp_row = new Matrix(1, Lms.getNumCols()+1);
                    L_tmp_row.set(Lms.getNumCols(), 0.0);
                    Matrix.extract(Lms, row, row+1, 0, Lms.getNumCols(), L_tmp_row, 0, 0);
                    L_tmp = Matrix.concatRows(L_tmp, L_tmp_row, null);
                  }
                }
                Matrix L_tmp_row = new Matrix(1, Lms.getNumCols()+1);
                L_tmp_row.set(Lms.getNumCols(), 1.0);
                Matrix.extract(Lms, i, i+1, 0, Lms.getNumCols(), L_tmp_row, 0, 0);
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

      Matrix Rchain = new Matrix(Qchain.getNumRows(), Qchain.getNumCols());
      for (int i = 0; i < Rchain.getNumRows(); i++) {
        for (int j = 0; j < Rchain.getNumCols(); j++) {
          Rchain.set(i, j, Qchain.get(i, j)/Xchain.get(j)/Vchain.get(i, j));
        }
      }

      for (int i = 0; i < infServers.size(); i++) {
        int row = infServers.get(i);
        for (int j = 0; j < Rchain.getNumCols(); j++) {
          Rchain.set(row, j, Lchain.get(row, j) / Vchain.get(row, j));
        }
      }

      Matrix Tchain = new Matrix(Vchain.getNumRows(), Vchain.getNumCols());
      for (int i = 0; i < Tchain.getNumRows(); i++) {
        for (int j = 0; j < Tchain.getNumCols(); j++) {
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
      NPFQN.npfqnNonexpApproxReturn NPFQNret = NPFQN.npfqn_nonexp_approx(options.config.highvar == null ? "interp": options.config.highvar,sn,ST0,V,SCV,T,U,gamma,nservers);
      ST = NPFQNret.ST;
      gamma = NPFQNret.gamma;
      eta = NPFQNret.eta;
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

  public static SolverNCLDReturn solver_ncld(NetworkStruct sn, SolverOptions options) {
    int M = sn.nstations;
    int K = sn.nclasses;
    Matrix nservers = sn.nservers;
    Matrix nserversFinite = nservers.clone();
    nserversFinite.removeInfinity();
    double minFiniteServer = Double.MAX_VALUE;
    for (int i=0; i<nservers.getNumElements(); i++) {
      if (Double.isFinite(nservers.get(i)) && nservers.get(i) < minFiniteServer) {
        minFiniteServer = nservers.get(i);
      }
    }
    if (minFiniteServer > 1) {
      if (sn.lldscaling.isEmpty() && M == 2 && Double.isFinite(sn.njobs.elementMaxAbs())) {
        double Nt = sn.njobs.elementSum();
        sn.lldscaling = sn.lldscaling.concatCols(new Matrix(M, (int) Nt));
        for (int i=0; i<M; i++) {
          for (int j=0; j<Nt; j++) {
            sn.lldscaling.set(i, j, Math.min(j + 1, sn.nservers.get(i)));
          }
        }
      } else {
        throw new RuntimeException("The load-dependent solver does not support multi-server stations yet. Specify multi-server stations via limited load-dependence.");
      }
    }
    if (sn.cdscaling != null && !sn.cdscaling.isEmpty() && options.method.equalsIgnoreCase("exact")) {
      throw new RuntimeException("Exact class-dependent solver not yet available in NC.");
    }
    Matrix NK = sn.njobs.transpose();
    if (Double.isInfinite(NK.elementMax())) {
      throw new RuntimeException("The load-dependent solver does not support open classes yet.");
    }
    int C = sn.nchains;
    Matrix SCV = sn.scv;
    Matrix gamma = new Matrix(M, 1);
    gamma.zero();
    Matrix V = new Matrix(sn.nstateful, K);
    for (int i = 0; i < sn.visits.size(); i++) {
      V = V.add(1, sn.visits.get(i));
    }
    Matrix ST = Matrix.ones(sn.rates.getNumRows(), sn.rates.getNumCols()).elementDiv(sn.rates);
    for (int i=0; i<ST.getNumRows(); i++) {
      for (int j=0; j<ST.getNumCols(); j++) {
        if (Double.isNaN(ST.get(i, j))) {
          ST.set(i, j, 0);
        }
      }
    }
    Matrix ST0 = ST.clone();
    Matrix lldscaling = sn.lldscaling;
    Matrix NKfinite = NK.clone();
    NKfinite.removeInfinity();
    double Nt = NKfinite.elementSum();
    if (lldscaling.isEmpty()) {
      lldscaling = Matrix.ones(M, (int) Math.ceil(Nt));
    }
    SN.snGetDemandsChainReturn demandsChainReturn = snGetDemandsChain(sn);
    Matrix Vchain = demandsChainReturn.Vchain;
    Matrix alpha = demandsChainReturn.alpha;
    Matrix eta_1 = new Matrix(1, M);
    eta_1.zero();
    Matrix eta = Matrix.ones(1, M);
    if (!sn.sched.containsValue(SchedStrategy.FCFS)) {
      options.iter_max = 1;
    }
    int iter = 0;
    long Tstart = System.currentTimeMillis();
    SN.snDeaggregateChainResultsReturn snDeaggragatedChains = null;
    Matrix lambda = null;
    Matrix Lchain = null;
    Matrix STchain = null;
    Matrix Q = null;
    Matrix U = null;
    Matrix R = null;
    Matrix T = null;
    Matrix X = null;
    Double lG = null;
    String method = null;
    while (Matrix.ones(eta.getNumRows(), eta.getNumCols()).sub(1, eta.elementDiv(eta_1)).elementMaxAbs() > options.iter_tol && iter < options.iter_max) {
      iter +=1;
      eta_1 = eta;
      M = sn.nstations;
      K = sn.nclasses;
      C = sn.nchains;
      Lchain = new Matrix(M, C);
      Lchain.zero();
      STchain = new Matrix(M, C);
      STchain.zero();
      Matrix SCVchain = new Matrix(M, C);
      SCVchain.zero();
      Matrix Nchain = new Matrix(1, C);
      Nchain.zero();
      Matrix refstatchain = new Matrix(C, 1);
      refstatchain.zero();
      for (int c=0; c<C; c++) {
        Matrix inchain = sn.inchain.get(c);
        boolean isOpenChain = false;
        for (int i=0; i<inchain.getNumElements(); i++) {
          if (Double.isInfinite(sn.njobs.get((int) inchain.get(i)))) {
            isOpenChain = true;
          }
        }
        for (int i=0; i<M; i++) {
          Matrix STinchain = new Matrix(1, inchain.getNumElements());
          Matrix alphainchain = new Matrix(1, inchain.getNumElements());
          Matrix SCVinchain = new Matrix(1, inchain.getNumElements());
          for (int j=0; j<inchain.getNumElements(); j++) {
            STinchain.set(j, ST.get(i, (int) inchain.get(j)));
            alphainchain.set(j, alpha.get(i, (int) inchain.get(j)));
            SCVinchain.set(j, SCV.get(i, (int) inchain.get(j)));
          }
          Lchain.set(i, c, Vchain.get(i, c) * STinchain.mult(alphainchain.transpose()).toDouble());
          STchain.set(i, c, STinchain.mult(alphainchain.transpose()).toDouble());
          if (isOpenChain && i == sn.refstat.get((int) inchain.get(0))) {
            Matrix STinchainFinite = STinchain.clone();
            STinchainFinite.removeInfinity();
            STchain.set(i, c, STinchainFinite.elementSum());
          }
          else {
            STchain.set(i, c, STinchain.mult(alphainchain.transpose()).toDouble());
          }
          SCVchain.set(i, c, SCVinchain.mult(alphainchain.transpose()).toDouble());
        }
        Matrix NKinchain = new Matrix(inchain.getNumElements());
        for (int i=0; i<inchain.getNumElements(); i++) {
          NKinchain.set(i, NK.get((int) inchain.get(i)));
        }
        Nchain.set(c, NKinchain.elementSum());
        refstatchain.set(c, sn.refstat.get((int) inchain.get(0)));
        if ((sn.refstat.get((int) inchain.get(0)) - refstatchain.get(c)) != 0) {
          throw new RuntimeException(String.format("Classes in chain %d have different reference station.", c));
        }
      }
      for (int i=0; i<STchain.getNumRows(); i++) {
        for (int j=0; j< STchain.getNumCols(); j++) {
          if (!Double.isFinite(STchain.get(i, j))) {
            STchain.set(i, j, 0);
          }
          if (!Double.isFinite(Lchain.get(i, j))) {
            Lchain.set(i, j, 0);
          }
        }
      }
      Tstart = System.currentTimeMillis();
      Matrix Nchainfinite = Nchain.clone();
      Nchainfinite.removeInfinity();
      Nt = Nchainfinite.elementSum();
      Matrix L = new Matrix(M, C);
      L.zero();
      Matrix mu = new Matrix(M, (int) Math.ceil(Nt));
      List<Integer> infServers = new ArrayList<>();
      Matrix Z = L.clone();
      for (int i=0; i<M; i++) {
        if (Double.isInfinite(nservers.get(i))) {
          infServers.add(i);
          for (int j=0; j<C; j++) {
            L.set(i, j, Lchain.get(i, j));
            Z.set(i, j, Lchain.get(i, j));
          }
          for (int j=0; j<Math.ceil(Nt); j++) {
            mu.set(i, j, j+1);
          }
        } else {
          if (options.method.equalsIgnoreCase("exact") && nservers.get(i) > 1) {
            System.out.println("Warning: SolverNC does not support exact multiserver yet. Switching to approximate method.");
          }
          for (int j=0; j<C; j++) {
            L.set(i, j, Lchain.get(i, j));
          }
          for (int j=0; j<Math.ceil(Nt); j++) {
            mu.set(i, j, lldscaling.get(i, j));
          }
        }
      }
      Matrix Qchain = new Matrix(M, C);
      Qchain.zero();
      Matrix Nchain0 = Nchain.clone();
      Nchain0.zero();
      PFQN.pfqnNcReturn ret = pfqn_ncld(L, Nchain, Nchain0, mu, options);
      lG = ret.lG;
      method = ret.method;
      Matrix Xchain = new Matrix(1, 0);
      if (Xchain.isEmpty()) {
        Matrix lGr = new Matrix(1, C);
        Matrix lGhat_fnci = new Matrix(1, C);
        Matrix lGhatir = new Matrix(1, C);
        Matrix lGr_i = new Matrix(1, C);
        Matrix ldDemand = new Matrix(M, C);
        for (int r=0; r<C; r++) {
          Matrix Nchain_r = Matrix.oner(Nchain, Collections.singletonList(r));
          lGr.set(r, pfqn_ncld(L, Nchain_r, Nchain0, mu, options).lG);
          Xchain = Xchain.concatCols(new Matrix(Math.exp(lGr.get(r) - lG)));
          for (int i=0; i<M; i++) {
            Qchain.set(i, r, 0);
          }
          Matrix CQchain_r = new Matrix(M, 1);
          CQchain_r.zero();
          if (M == 2 && Double.isInfinite(sn.nservers.elementMaxAbs())) {
            int firstDelay = -1;
            for (int i=0; i<sn.nservers.getNumElements(); i++) {
              if (Double.isInfinite(sn.nservers.get(i))) {
                firstDelay = i;
                break;
              }
            }
            Qchain.set(firstDelay, r, Lchain.get(firstDelay, r) * Xchain.get(r));
            for (int i=0; i<M; i++) {
              if (i != firstDelay) {
                Qchain.set(i, r, Nchain.get(r) - Lchain.get(firstDelay, r) * Xchain.get(r));
              }
            }
          } else {
            for (int i=0; i<M; i++) {
              Matrix Lms_i = Matrix.extractRows(L, i, i+1, null);
              Matrix mu_i = Matrix.extractRows(mu, i, i+1, null);
              Matrix muhati = pfqn_mushift(mu, i);
              pfqnFncReturn fncret = pfqn_fnc(Matrix.extractRows(muhati, i, i+1, null));
              Matrix muhati_f = fncret.mu;
              double c = fncret.c.toDouble();
              if (Lchain.get(i, r) > 0) {
                if (Double.isInfinite(nservers.get(i))) {
                  Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r));
                }
                else {
                  if (i == (M - 1) && nserversFinite.elementSum() == 1) {
                    double Lchainsum = 0;
                    double Qchainsum = 0;
                    for (int j=0; j<nservers.getNumElements(); j++) {
                      if (Double.isInfinite(nservers.get(j))) {
                        Lchainsum += Lchain.get(j, r);
                      } else {
                        Qchainsum += Qchain.get(j, r);
                      }
                    }
                    Qchain.set(i, r, Math.max(0, Nchain.get(r) - Lchainsum * Xchain.get(r) - Qchainsum));
                  } else {
                    lGhat_fnci.set(r, pfqn_ncld(Matrix.concatRows(L, Matrix.extractRows(L, i, i+1, null), null), Nchain_r, Nchain0, Matrix.concatRows(muhati, muhati_f, null), options).lG);
                    lGhatir.set(r, pfqn_ncld(L, Nchain_r, Nchain0, muhati, options).lG);
                    lGr_i.set(r, pfqn_ncld(Lms_i, Nchain_r, Nchain0, muhati, options).lG);
                    double dlGa = lGhat_fnci.get(r) - lGhatir.get(r);
                    double dlG_i = lGr_i.get(r) - lGhatir.get(r);
                    CQchain_r.set(i, (Math.exp(dlGa) - 1) + c * (Math.exp(dlG_i) - 1));
                    ldDemand.set(i, r, Math.log(L.get(i, r)) + lGhatir.get(r) - Math.log(mu.get(i, 1)) - lGr.get(r));
                    Qchain.set(i, r, Math.exp(ldDemand.get(i, r)) * Xchain.get(r) * (1 + CQchain_r.get(i)));
                  }
                }
              }
            }
          }
        }
      } else {
        for (int r=0; r<C; r++) {
          for (int i=1; i<M; i++) {
            if (Lchain.get(i, r) > 0) {
              if (Double.isInfinite(nservers.get(i))) {
                Qchain.set(i, r, Lchain.get(i, r) * Xchain.get(r));
              }
            }
          }
        }
      }
      if (Arrays.stream(Xchain.toArray1D()).allMatch(Double::isNaN)) {
        System.out.println("Warning: Normalizing constant computations produced a floating-point range exception. Model is likely too large.");
      }
      Z = Z.sumCols();
      Matrix Rchain = Qchain.element_divide(Xchain.repmat(M, 1)).element_divide(Vchain);
      for (int i : infServers) {
        for (int j=0; j<Rchain.getNumCols(); j++) {
          Rchain.set(i, j, Lchain.get(i, j) / Vchain.get(i, j));
        }
      }
      Matrix Tchain = Xchain.repmat(M, 1).elementMult(Vchain, null);
      Matrix Uchain = Tchain.elementMult(Lchain, null);
      Matrix Cchain = Nchain.elementDiv(Xchain).sub(1, Z);
      snDeaggragatedChains = snDeaggregateChainResults(sn, Lchain, ST, STchain, Vchain, alpha, null, null, Rchain, Tchain, null, Xchain);
      Q = snDeaggragatedChains.Q;
      U = snDeaggragatedChains.U;
      R = snDeaggragatedChains.R;
      T = snDeaggragatedChains.T;
      C = (int) snDeaggragatedChains.C.get(0);
      X = snDeaggragatedChains.X;
      NPFQN.npfqnNonexpApproxReturn NPFQNret = NPFQN.npfqn_nonexp_approx(options.config.highvar == null ? "default": options.config.highvar,sn,ST0,V,SCV,T,U,gamma,nservers);
      ST = NPFQNret.ST;
      gamma = NPFQNret.gamma;
      eta = NPFQNret.eta.transpose();
    }
    SN.snGetProductFormChainParamsReturn snProductForm = snGetProductFormParams(sn);
    double runtime = (System.currentTimeMillis() - Tstart) / 1000.0;
    Q = snDeaggragatedChains.Q;
    Q.abs();
    R = snDeaggragatedChains.R;
    R.abs();
    U = snDeaggragatedChains.U;
    U.abs();
    List<Integer> openClasses = new ArrayList<>();
    for (int i=0; i<M; i++) {
      for (int j=0; j<NK.getNumElements(); j++) {
          if (Double.isInfinite(NK.get(j))) {
            openClasses.add(j);
          }
        }
      if (sn.nservers.get(i) > 1 && sn.nservers.get(i) < Double.POSITIVE_INFINITY) {
        for (int r=0; r<K; r++) {
          Matrix c = Matrix.extractColumn(sn.chains, r, null).find();
          if ((!openClasses.contains(r) && snDeaggragatedChains.X.get(r) > 0) || (openClasses.contains(r) && snProductForm.lambda.get(r) > 0)) {
              U.set(i, r, snDeaggragatedChains.X.get(r) * ST.get(i, r) / sn.nservers.get(i));
              for (int j = 0; j < c.getNumElements(); j++) {
                U.set(i, r, U.get(i, r) * sn.visits.get(j).get(i, r) / sn.visits.get(j).get((int) sn.refstat.get(r), r));
              }
            }
          }
        } else {

        for (int j=0; j<U.getNumCols(); j++) {
          U.set(i, j, U.get(i, j) / Matrix.extractRows(U, i, i+1, null).elementMax());
        }
        if (Matrix.extractRows(U, i, i+1, null).elementSum() > 1) {
          Matrix Uinotnan = Matrix.extractRows(U, i, i+1, null);
          Uinotnan.removeNaN();
          for (int j=0; j<U.getNumCols(); j++) {
            U.set(i, j, U.get(i, j) / U.elementSum());
          }
        }
      }
        // Continue from line 213
    }
    X = snDeaggragatedChains.X;
    X.apply(Double.POSITIVE_INFINITY, 0, "equal");
    X.apply(Double.NaN, 0, "equal");
    U.apply(Double.POSITIVE_INFINITY, 0, "equal");
    U.apply(Double.NaN, 0, "equal");
    Q.apply(Double.POSITIVE_INFINITY, 0, "equal");
    Q.apply(Double.NaN, 0, "equal");
    R.apply(Double.POSITIVE_INFINITY, 0, "equal");
    R.apply(Double.NaN, 0, "equal");
    return new SolverNCLDReturn(Q, U, R, snDeaggragatedChains.T, C, X, lG, runtime, iter, method);
  }

  /**
   * Returns the feature set supported by the NC solver
   * @return - the feature set supported by the NC solver
   */
  public static FeatureSet getFeatureSet(){
    FeatureSet s = new FeatureSet();
    // TODO: update with the features supported by JLINE. These are the features supported by LINE.
    String[] features = {"Sink", "Source",
        "ClassSwitch", "Delay", "Queue",
        "APH", "Coxian", "Erlang", "Det", "Exp", "HyperExp",
        "StatelessClassSwitcher", "InfiniteServer",
        "SharedServer", "Buffer", "Dispatcher",
        "Server", "Sink", "RandomSource", "ServiceTunnel",
        "SchedStrategy_INF", "SchedStrategy_PS", "SchedStrategy_SIRO",
        "RoutingStrategy_PROB", "RoutingStrategy_RAND",
        "SchedStrategy_FCFS", "ClosedClass", "ClosedClass",
        "Cache", "CacheClassSwitcher", "OpenClass"};
    s.setTrue(features);
    return s;
  }

  /**
   * Checks whether the given model is supported by the NC solver
   * @param model - the network model
   * @return - true if the model is supported, false otherwise
   */
  public static boolean supports(Network model){
    FeatureSet featUsed = model.getUsedLangFeatures();
    FeatureSet featSupported = SolverNC.getFeatureSet();
    return FeatureSet.supports(featSupported, featUsed);
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

  public static class SolverNCLDReturn {
    public Matrix Q;
    public Matrix U;
    public Matrix R;
    public Matrix T;
    public int C;
    public Matrix X;
    public double lG;
    public double runtime;
    public int it;
    public String method;

    public SolverNCLDReturn(Matrix Q, Matrix U, Matrix R, Matrix T, int C,
                            Matrix X, double lG, double runtime, int it, String method) {
      this.Q = Q;
      this.U = U;
      this.R = R;
      this.T = T;
      this.C = C;
      this.X = X;
      this.lG = lG;
      this.runtime = runtime;
      this.it = it;
      this.method = method;
    }
  }

}
