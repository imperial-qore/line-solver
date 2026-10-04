/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ctmc;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * CTMC is an alias for SolverCTMC (Continuous Time Markov Chain solver).
 */
public class CTMC extends SolverCTMC {

    public CTMC(Network model, String method) {
        super(model, method);
    }

    public CTMC(Network model, SolverOptions options) {
        super(model, options);
    }

    public CTMC(Network model, Object... varargin) {
        super(model, varargin);
    }

    public CTMC(Network model) {
        super(model);
    }
}
