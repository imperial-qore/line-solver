/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.des;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * DES is an alias for SolverDES (Discrete Event Simulation solver).
 */
public class DES extends SolverDES {

    public DES(Network model, String method) {
        super(model, method);
    }

    public DES(Network model, SolverOptions options) {
        super(model, options);
    }

    public DES(Network model, Object... varargin) {
        super(model, varargin);
    }

    public DES(Network model) {
        super(model);
    }
}
