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

    /**
     * A factory whose layer solver is NAMED, which a lambda cannot be asked for.
     *
     * <p>Routed call groups (synchCallRoundRobin, synchCallJSQ) are refused under
     * a layer solver without state-dependent routing, and that is a property of
     * the factory rather than of the model, so nothing the ensemble carries can
     * answer it. MATLAB reads the name off the function handle's text; a Java
     * lambda has no such text, so the type is passed here instead.
     */
    public LN(LayeredNetwork model, SolverFactory solverFactory, SolverOptions options,
              SolverType layerSolverType) {
        super(model, solverFactory, options, layerSolverType);
    }
}
