/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.uq;

import jline.lang.Network;
import jline.solvers.SolverOptions;

/**
 * UQ is an alias for {@link SolverUQ} (Bayesian-style parameter
 * uncertainty analysis). The wrapper expands a model containing a {@link
 * jline.lang.processes.Prior} distribution into one network per alternative,
 * runs an inner solver on each, and aggregates the results into prior-weighted
 * means and per-alternative empirical posterior distributions.
 */
public class UQ extends SolverUQ {

    public UQ(Network model, SolverFactory solverFactory) {
        super(model, solverFactory);
    }

    public UQ(Network model, SolverFactory solverFactory, SolverOptions options) {
        super(model, solverFactory, options);
    }
}
