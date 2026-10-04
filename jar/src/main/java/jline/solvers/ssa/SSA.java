/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ssa;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * SSA is an alias for SolverSSA (Stochastic State-space Analysis solver).
 */
public class SSA extends SolverSSA {

    public SSA(Network model, String method) {
        super(model, method);
    }

    public SSA(Network model, SolverOptions options) {
        super(model, options);
    }

    public SSA(Network model, Object... varargin) {
        super(model, varargin);
    }

    public SSA(Network model) {
        super(model);
    }
}
