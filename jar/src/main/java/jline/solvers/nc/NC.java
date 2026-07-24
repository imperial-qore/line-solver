/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * NC is an alias for SolverNC (Normalizing Constant solver).
 */
public class NC extends SolverNC {

    public NC(Network model, String method) {
        super(model, method);
    }

    public NC(Network model, SolverOptions options) {
        super(model, options);
    }

    public NC(Network model, Object... varargin) {
        super(model, varargin);
    }

    public NC(Network model) {
        super(model);
    }
}
