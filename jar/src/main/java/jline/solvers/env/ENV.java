/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.env;

import jline.lang.Environment;
import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.ln.SolverFactory;

/**
 * ENV is an alias for SolverENV (Ensemble environment solver).
 */
public class ENV extends SolverENV {

    public ENV(Environment renv, NetworkSolver[] solvers) {
        super(renv, solvers);
    }

    public ENV(Environment renv, NetworkSolver[] solvers, SolverOptions options) {
        super(renv, solvers, options);
    }

    /**
     * Creates an ENV solver from a factory that produces an inner solver for
     * each stage sub-model.
     *
     * @param renv    the random-environment model
     * @param factory a factory mapping a stage model to its inner solver
     */
    public ENV(Environment renv, SolverFactory factory) {
        super(renv, buildSolvers(renv, factory));
    }

    /**
     * Creates an ENV solver from a per-stage solver factory, with options.
     *
     * @param renv    the random-environment model
     * @param factory a factory mapping a stage model to its inner solver
     * @param options the solver options
     */
    public ENV(Environment renv, SolverFactory factory, SolverOptions options) {
        super(renv, buildSolvers(renv, factory), options);
    }

    private static NetworkSolver[] buildSolvers(Environment renv, SolverFactory factory) {
        Network[] models = renv.getEnsemble().toArray(new Network[0]);
        NetworkSolver[] solvers = new NetworkSolver[models.length];
        for (int e = 0; e < models.length; e++) {
            solvers[e] = factory.at(models[e]);
        }
        return solvers;
    }
}