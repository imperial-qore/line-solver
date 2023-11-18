package jline.solvers.mva.handlers;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;

public interface MVASolverHandler {
    SolverMVAResult solve(NetworkStruct sn, SolverOptions options);
}
