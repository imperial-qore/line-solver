/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ba;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * BA is an alias for SolverBA (Bound Analysis solver).
 */
public class BA extends SolverBA {

    public BA(Network model, String method) {
        super(model, method);
    }

    public BA(Network model, SolverOptions options) {
        super(model, options);
    }

    public BA(Network model, Object... varargin) {
        super(model, varargin);
    }

    public BA(Network model) {
        super(model);
    }
}
