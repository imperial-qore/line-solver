package jline.solvers.ln;

import jline.lang.Network;
import jline.solvers.NetworkSolver;

public interface SolverFactory {
    NetworkSolver at(Network model);     // return solver associated to model
}
