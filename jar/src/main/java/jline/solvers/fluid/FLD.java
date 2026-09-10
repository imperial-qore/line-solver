/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * FLD is an alias for SolverFluid (Fluid/Mean-Field Approximation solver).
 */
public class FLD extends SolverFluid {

    public FLD(Network model, String method) {
        super(model, method);
    }

    public FLD(Network model, SolverOptions options) {
        super(model, options);
    }

    public FLD(Network model, Object... varargin) {
        super(model, varargin);
    }

    public FLD(Network model) {
        super(model);
    }
}
