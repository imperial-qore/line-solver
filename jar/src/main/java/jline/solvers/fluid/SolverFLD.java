/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * SolverFLD is an alias for SolverFluid (Fluid/Mean-Field Approximation solver).
 * It mirrors the MATLAB class name, where SolverFLD is the primary class and
 * SolverFluid the alias.
 */
public class SolverFLD extends SolverFluid {

    public SolverFLD(Network model, String method) {
        super(model, method);
    }

    public SolverFLD(Network model, SolverOptions options) {
        super(model, options);
    }

    public SolverFLD(Network model, Object... varargin) {
        super(model, varargin);
    }

    public SolverFLD(Network model) {
        super(model);
    }
}
