/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;

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
    
    /** Number of iterations performed by the MVA algorithm (used by runAnalyzer method) */
    public int iter;
    
    /** Cache hit probabilities [items x classes] (used by cache analyzers) */
    public Matrix hitProb;
    
    /** Cache miss probabilities [items x classes] (used by cache analyzers) */
    public Matrix missProb;
}
