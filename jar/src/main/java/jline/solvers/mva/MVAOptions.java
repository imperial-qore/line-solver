/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import jline.lang.constant.SolverType;
import jline.solvers.SolverOptions;

/**
 * Configuration options for Mean Value Analysis (MVA) solver.
 * 
 * <p>MVAOptions extends the base solver options to provide MVA-specific
 * configuration parameters. MVA is an exact analytical method for computing
 * steady-state performance measures of closed queueing networks with
 * product-form solutions.</p>
 * 
 * <p>Key MVA characteristics:
 * <ul>
 *   <li>Exact results for product-form networks</li>
 *   <li>Efficient for moderate network sizes</li>
 *   <li>Supports multi-class closed networks</li>
 *   <li>Handles load-independent and load-dependent servers</li>
 * </ul>
 * </p>
 * 
 * @see SolverMVA
 * @see SolverOptions
 * @since 1.0
 */
public class MVAOptions extends SolverOptions {
    
    /**
     * Constructs a new MVAOptions instance with default MVA solver configuration.
     */
    public MVAOptions() {
        super(SolverType.MVA);
    }
}
