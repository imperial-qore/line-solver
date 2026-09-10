/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.auto;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * AUTO is an alias for SolverAUTO (Automatic solver selection).
 */
public class AUTO extends SolverAUTO {

    public AUTO(Network model, SolverOptions options) {
        super(model, options);
    }

    public AUTO(Network model, Object... varargin) {
        super(model, varargin);
    }

    public AUTO(Network model) {
        super(model);
    }
}
