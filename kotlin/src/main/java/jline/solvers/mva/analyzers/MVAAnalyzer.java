package jline.solvers.mva.analyzers;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVAResult;

public interface MVAAnalyzer {
    void analyze(NetworkStruct sn, SolverOptions options, SolverMVAResult res);
}
