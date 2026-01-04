/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.env;

import jline.lang.Environment;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;

/**
 * ENV is an alias for SolverENV (Ensemble environment solver).
 */
public class ENV extends SolverENV {

    public ENV(Environment renv, NetworkSolver[] solvers) {
        super(renv, solvers);
    }

    public ENV(Environment renv, NetworkSolver[] solvers, SolverOptions options) {
        super(renv, solvers, options);
    }
}