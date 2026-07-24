/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

public class NCResult extends SolverResult {
    // solver/method/iter are inherited from SolverResult, not redeclared: see _kb/11-conventions-and-gotchas.md
    public Prob prob;
    public double lG;
    public Matrix STeff;
    public int it;
    public Matrix pij;  // Item probabilities for cache analyzer

    // hitProb/missProb are used only by SolverNCCacheQNAnalyzer
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
