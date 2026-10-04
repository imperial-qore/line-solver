/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ag;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * AG is an alias for SolverAG (RCAT agent-decomposition solver).
 */
public class AG extends SolverAG {

    public AG(Network model, String method) {
        super(model, method);
    }

    public AG(Network model, SolverOptions options) {
        super(model, options);
    }

    public AG(Network model, Object... varargin) {
        super(model, varargin);
    }

    public AG(Network model) {
        super(model);
    }
}
