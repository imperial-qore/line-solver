package jline.solvers.nc;

import jline.solvers.SolverResult;
import jline.util.Matrix;

public class SolverNCResult extends SolverResult {
  public String solver;
  public Prob prob;
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

  /*
   * The following fields are only used by the SolverNCCacheQNAnalyzer
   */
  public Matrix hitProb;
  public Matrix missProb;


  class Prob {
    public Double logNormConstAggr;
    public Matrix marginal;
  }
}
