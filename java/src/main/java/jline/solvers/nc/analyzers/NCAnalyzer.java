package jline.solvers.nc.analyzers;

import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.solvers.nc.SolverNCResult;

public interface NCAnalyzer {
  void analyze(NetworkStruct sn, SolverOptions options, SolverNCResult res);

}
