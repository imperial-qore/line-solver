/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.lqns;

import jline.lang.layered.LayeredNetwork;
import jline.solvers.SolverOptions;

/**
 * LQNS is an alias for SolverLQNS (Layered Queueing Network Solver).
 */
public class LQNS extends SolverLQNS {

    public LQNS(LayeredNetwork model, SolverOptions options) {
        super(model, options);
    }

    public LQNS(LayeredNetwork model) {
        super(model);
    }
}
