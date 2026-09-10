/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

import java.util.HashMap;
import java.util.Map;

/**
 * Result container for Mean Value Analysis (MVA) solver computations.
 * 
 * <p>MVAResult extends the base SolverResult to include MVA-specific
 * performance metrics and intermediate computation results. This includes
 * normalizing constants, iteration counts, and cache-specific probabilities.</p>
 * 
 * @see SolverMVA
 * @see SolverResult
 * @since 1.0
 */
public class MVAResult extends SolverResult {
    
    /** Logarithm of the aggregate normalizing constant */
    public double logNormConstAggr;
    
    // The iteration count lives on SolverResult.iter. It used to be redeclared
    // here as well, which SHADOWED the base field: setAvgResults writes the base
    // one through a SolverResult reference while the analyzer wrote this one, so
    // a caller that cast the result to MVAResult -- the natural thing to do --
    // read 0 on every method. Do not reintroduce it.
    /**
     * Whether the fixed point met its tolerance, or {@code null} when the handler
     * does not report one. This is NOT {@code iter < iter_max}: in Solver_amvald
     * {@code totiter} aggregates the NESTED inner sweeps, so it routinely reaches
     * the cap on a solve whose outer residual is exactly 0. Only the residual
     * decides convergence there.
     *
     * <p>Left null by the single-loop handlers (the pfqn_linearizer family,
     * pfqn_bs, pfqn_aql), for which the count IS a sound signal because it is not
     * an aggregate: those exit on {@code norm(Q-Qlast)<tol || iter > maxiter}, so
     * reaching maxiter does mean the tolerance was never met. Callers should
     * prefer this flag and fall back to the count only when it is null.
     */
    public Boolean converged = null;
    
    /** Cache hit probabilities [items x classes] (used by cache analyzers) */
    public Matrix hitProb;
    
    /** Cache miss probabilities [items x classes] (used by cache analyzers) */
    public Matrix missProb;

    /** Per-cache per-item occupancy [nitems x (lists+1)], keyed by node index
     *  (used by the cache+queueing analyzer to populate getAvgItemTable) */
    public Map<Integer, Matrix> cacheItemProb = new HashMap<Integer, Matrix>();
}
