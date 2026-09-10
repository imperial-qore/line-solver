/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ln;

import jline.lang.constant.SolverType;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.SolverOptions;

/**
 * LN is an alias for SolverLN (Layered Network solver).
 */
public class LN extends SolverLN {

    public LN(LayeredNetwork model, SolverOptions options) {
        super(model, options);
    }

    public LN(LayeredNetwork model) {
        super(model);
    }

    public LN(LayeredNetwork model, SolverType solverType) {
        super(model, solverType);
    }

    public LN(LayeredNetwork model, SolverType solverType, SolverOptions options) {
        super(model, solverType, options);
    }

    public LN(LayeredNetwork model, SolverFactory solverFactory) {
        super(model, solverFactory);
    }

    public LN(LayeredNetwork model, SolverFactory solverFactory, SolverOptions options) {
        super(model, solverFactory, options);
    }
}
