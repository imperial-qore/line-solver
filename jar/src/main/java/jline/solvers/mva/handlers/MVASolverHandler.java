/*
 * Copyright (c) 2012-2025, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.mva.handlers;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mva.MVAResult;

public interface MVASolverHandler {
    MVAResult solve(NetworkStruct sn, SolverOptions options);
}
