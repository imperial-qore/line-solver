package jline.solvers.mva;

import jline.solvers.SolverResult;
import jline.util.Matrix;

public class SolverMVAResult extends SolverResult{
	  public double logNormConstAggr;
	  /*
	  * The following two fields are used only by the runAnalyzer() method of the MVA solver.
	  * */
	  public Matrix AN;
	  public Matrix WN;
	  public int iter;
	/*
	* The following fields are only used by the SolverMVACacheQNAnalyzer
	 */
	public Matrix hitProb;
	public Matrix missProb;
}
