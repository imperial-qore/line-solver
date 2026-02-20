/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

public class NCResult extends SolverResult {
    public String solver;
    public Prob prob;
    public double lG;
    public Matrix STeff;
    public int it;
    public String method;
    public Matrix pij;  // Item probabilities for cache analyzer
    public int iter;    // Iteration count

    /*
     * The following fields are only used by the SolverNCCacheQNAnalyzer
     */
    public Matrix hitProb;
    public Matrix missProb;

    public NCResult() {
        this.prob = new Prob();
    }

    class Prob {
        public Double logNormConstAggr;
        public Matrix marginal;
        public Double joint;
        public Matrix itemProb;
    }
}
