/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * Fluid is an alias for SolverFluid (Fluid/Mean-Field Approximation solver).
 */
public class Fluid extends SolverFluid {

    public Fluid(Network model, String method) {
        super(model, method);
    }

    public Fluid(Network model, SolverOptions options) {
        super(model, options);
    }

    public Fluid(Network model, Object... varargin) {
        super(model, varargin);
    }

    public Fluid(Network model) {
        super(model);
    }
}
