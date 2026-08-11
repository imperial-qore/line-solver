/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * LDES is an alias for SolverLDES (Discrete Event Simulation solver).
 */
public class LDES extends SolverLDES {

    public LDES(Network model, String method) {
        super(model, method);
    }

    public LDES(Network model, SolverOptions options) {
        super(model, options);
    }

    public LDES(Network model, Object... varargin) {
        super(model, varargin);
    }

    public LDES(Network model) {
        super(model);
    }
}
