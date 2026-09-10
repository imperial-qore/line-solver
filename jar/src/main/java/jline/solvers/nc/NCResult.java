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
    /** Product-form certificate of a Petri net, set by Solver_nc_spn_analyzer. */
    public jline.api.spn.Spn_pf.SpnPfResult spnpf;

    /** Per-item occupancy of each cache node, (nitems x lists+1) keyed by node index; column 0 = miss. */
    public java.util.Map<Integer, Matrix> cacheItemProb = new java.util.HashMap<Integer, Matrix>();

    public NCResult() {
        this.prob = new Prob();
    }

    /**
     * Log normalizing constant of the analyzed model, or NaN when none was
     * computed. The analyzers park it in the package-private {@code prob} block
     * and leave the flat {@link #lG} at its zero default, so a caller outside
     * this package that reads the field alone gets 0 for every solve.
     *
     * @return log g(N), or NaN
     */
    public double logNormConstAggr() {
        if (this.prob != null && this.prob.logNormConstAggr != null) {
            return this.prob.logNormConstAggr.doubleValue();
        }
        return this.lG == 0 ? Double.NaN : this.lG;
    }

    class Prob {
        public Double logNormConstAggr;
        public Matrix marginal;
        public Double joint;
        public Matrix itemProb;
    }
}
