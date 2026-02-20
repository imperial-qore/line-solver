/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ssa;

import jline.lang.constant.SolverType;
import jline.solvers.SolverOptions;

/**
 * Configuration options for Stochastic Simulation Algorithm (SSA) solver.
 *
 * <p>SSAOptions extends the base solver options to provide SSA-specific
 * configuration parameters. SSA is a discrete-event simulation method that
 * generates sample paths of the queueing network state evolution over time.</p>
 * 
 * <p>Key SSA characteristics:
 * <ul>
 *   <li>Stochastic simulation approach</li>
 *   <li>Handles general service and arrival processes</li>
 *   <li>Supports both transient and steady-state analysis</li>
 *   <li>Can analyze complex networks beyond product-form</li>
 *   <li>Provides statistical estimates with confidence intervals</li>
 * </ul>
 * </p>
 * 
 * @see SolverSSA
 * @see SolverOptions
 * @since 1.0
 */
public class SSAOptions extends SolverOptions {
    
    /**
     * Constructs a new SSAOptions instance with default SSA solver configuration.
     */
    public SSAOptions() {
        super(SolverType.SSA);
    }
}
