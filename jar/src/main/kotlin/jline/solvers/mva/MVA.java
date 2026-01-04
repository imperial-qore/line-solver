/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * MVA is an alias for SolverMVA (Mean Value Analysis solver).
 */
public class MVA extends SolverMVA {

    public MVA(Network model, String method) {
        super(model, method);
    }

    public MVA(Network model, SolverOptions options) {
        super(model, options);
    }

    public MVA(Network model, Object... varargin) {
        super(model, varargin);
    }

    public MVA(Network model) {
        super(model);
    }
}
