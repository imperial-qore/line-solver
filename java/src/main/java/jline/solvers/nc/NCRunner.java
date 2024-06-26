package jline.solvers.nc;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodes.Cache;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.analyzers.*;
import jline.util.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class NCRunner {

  protected NetworkStruct sn;
  protected SolverOptions options;
  protected SolverNCResult res;
  private final SolverNC solver;

  public NCRunner(SolverNC solver) {
    this.solver = solver;
    this.sn = solver.sn;
    this.options = solver.options;
    this.res = null;
  }

  /**
   * runAnalyzer() method from LINE.
   * This method is used to run the solver on the given network.
   * It returns the performance measures corresponding to the given network.
   * The method is called from the SolverNC class.
   *
   * @return - the performance measures corresponding to the given network
   */
  public SolverResult run() {
    long T0 = System.currentTimeMillis();
    int iter = 0;
    if (this.solver.enableChecks && !SolverNC.supports(this.solver.model)) {
      // TODO: not implemented
      throw new RuntimeException("This model contains features not supported by the solver.");
    }
    SolverNCResult ret = new SolverNCResult();
    NCAnalyzer analyzer;
    boolean wasDefault = false;
    // Method Selection and Preprocessing
    switch (options.method) {
      case "default":
        if (sn.nstations == 2 && !sn.nodetypes.contains(NodeType.Cache) && sn.nodetypes.contains(NodeType.Delay)) {
          boolean hasFiniteStations = false;
          for (int i = 0; i < sn.nservers.getNumRows(); i++) {
            for (int j = 0; j < sn.nservers.getNumCols(); j++) {
              if (Double.isFinite(sn.nservers.get(i, j)) && sn.nservers.get(i, j) > 1) {
                hasFiniteStations = true;
                break;
              }
            }
          }
          if (hasFiniteStations) {
            options.method = "comomld";
            wasDefault = true;
          }
        }
        break;
      case "exact":
        if (!solver.model.hasProductFormSolution()) {
          throw new RuntimeException("The exact method requires the model to have a product-form solution. This model does not have one. You can use Network.hasProductFormSolution() to check before running the solver.");
        } else if (sn.lldscaling.isEmpty()) {
          double Nt = sn.njobs.elementSum();
          if (Double.isFinite(Nt)) {
            sn.lldscaling = Matrix.ones(sn.nstations, (int) Nt);
            for (int i = 0; i < sn.nstations; i++) {
              if (sn.nservers.get(i) > 1 && Double.isFinite(sn.nservers.get(i))) {
                for (int j = 0; j < Nt; j++) {
                  sn.lldscaling.set(i, j, Math.min(j + 1, sn.nservers.get(i)));
                }
                sn.nservers.set(i, 1);
              }
            }
          }
        }
        break;
    }
    solver.resetRandomGeneratorSeed(options.seed);

    // Special Handling for Non-Reentrant Nodes: the system model meets certain criteria 
    // (e.g., no closed jobs and contains only source, cache, and sink nodes)
    List<NodeType> nonReentrant = new ArrayList<>(Arrays.asList(NodeType.Source, NodeType.Cache, NodeType.Sink));
    if (sn.nclosedjobs == 0 && sn.nodetypes.size() == 3 && sn.nodetypes.containsAll(nonReentrant)) {
      for (int ind = 0; ind < sn.nnodes; ind++) {
        // Cache Node Handling: For each cache node, initialize hit and miss probabilities and analyze the network.
        if (sn.nodetypes.get(ind) == NodeType.Cache) {
          Cache cacheNode = (Cache) solver.model.getNodes().get(ind);
          Matrix prob = new Matrix(cacheNode.getHitClass());
          for (int i = 0; i < prob.getNumRows(); i++) {
            for (int j = 0; j < prob.getNumCols(); j++) {
              if (prob.get(i, j) > 0) {
                prob.set(i, j, 0.5);
              }
            }
          }
          cacheNode.setResultHitProb(prob);
          Matrix missProb = new Matrix(prob.getNumRows(), prob.getNumCols());
          for (int i = 0; i < prob.getNumRows(); i++) {
            for (int j = 0; j < prob.getNumCols(); j++) {
              missProb.set(i, j, 1 - prob.get(i, j));
            }
          }
          cacheNode.setResultMissProb(missProb);
        }
      }
      solver.model.refreshChains(true);
      analyzer = new SolverNCAnalyzer();
      analyzer.analyze(this.sn, this.options, ret);
      for (int ind = 0; ind < this.sn.nnodes; ind++) {
        if (this.sn.nodetypes.get(ind) == NodeType.Cache) {
          Cache cacheNode = (Cache) solver.model.getNodes().get(ind);
          Matrix hitClass = cacheNode.getHitClass();
          Matrix missClass = cacheNode.getMissClass();
          Matrix hitProb = new Matrix(1, hitClass.length());
          for (int k = 0; k < hitClass.length(); k++) {
            int chain_k = 0;
            for (; chain_k < this.sn.chains.getNumRows(); chain_k++) {
              if (this.sn.chains.get(chain_k, k) > 0)
                break;
            }
            Matrix inchain = new Matrix(1, this.sn.chains.getNumCols());
            for (int i = 0; i < inchain.getNumCols(); i++) {
              inchain.set(0, i, this.sn.chains.get(chain_k, i) > 0 ? 1 : 0);
            }
            int h = (int) hitClass.get(k);
            int m = (int) missClass.get(k);
            if (h > -1 && m > -1) {
              double sumXN = 0;
              for (int i = 0; i < inchain.getNumCols(); i++) {
                if (inchain.get(i) > 0 && !Double.isNaN(ret.XN.get(i))) {
                  sumXN += ret.XN.get(i);
                }
              }
              hitProb.set(k, ret.XN.get(h) / sumXN);
            }
          }
          Matrix missProb = new Matrix(1, hitClass.length());
          for (int i = 0; i < hitClass.length(); i++) {
            missProb.set(i, 1 - hitProb.get(i));
          }
          cacheNode.setResultHitProb(hitProb);
          cacheNode.setResultMissProb(missProb);
        }
      }
      this.solver.model.refreshChains(true);
    } else {
      if (sn.nodetypes.contains(NodeType.Cache)) {
        analyzer = new SolverNCCacheQNAnalyzer();
        analyzer.analyze(this.sn, this.options, ret);
        for (int ind = 0; ind < sn.nnodes; ind++) {
          if (sn.nodetypes.get(ind) == NodeType.Cache) {
            Cache cacheNode = (Cache) solver.model.getNodes().get(ind);
            cacheNode.setResultHitProb(Matrix.extractRows(ret.hitProb, ind, ind + 1, null));
            cacheNode.setResultMissProb(Matrix.extractRows(ret.missProb, ind, ind + 1, null));
          }
        }
        solver.model.refreshChains(true);
      } else {
        if ((!(sn.lldscaling == null) && !sn.lldscaling.isEmpty()) || (!(sn.cdscaling == null) && !sn.cdscaling.isEmpty())) {
          analyzer = new SolverNCLDAnalyzer();
          analyzer.analyze(this.sn, this.options, ret);
        } else {
          switch (options.method) {
            case "exact":
              if (!solver.model.hasOpenClasses()) {
                analyzer = new SolverNCLDAnalyzer();
                analyzer.analyze(this.sn, this.options, ret);
              } else {
                analyzer = new SolverNCAnalyzer();
                analyzer.analyze(this.sn, this.options, ret);
              }
              break;
            case "rd": case "nrp": case "nr.probit": case "nrl": case "nr.logit": case "comomld":
              analyzer = new SolverNCLDAnalyzer();
              analyzer.analyze(this.sn, this.options, ret);
            default:
              analyzer = new SolverNCAnalyzer();
              analyzer.analyze(this.sn, this.options, ret);
          }
        }

      }
    }
    return ret;
  }
}
