/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.wrappers.lqns;

import jline.lang.constant.SolverType;
import jline.solvers.SolverOptions;

public class LQNSOptions extends SolverOptions {
    public LQNSOptions() {
        super(SolverType.LQNS);
    }
}
