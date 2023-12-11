package jline.solvers.ctmc;

import jline.solvers.ssa.state.SSAStateMatrix;
import jline.util.Matrix;

import java.util.ArrayList;

public class SolverCTMCResult {
    protected ArrayList<SSAStateMatrix> stateSpace;
    protected ArrayList<EventData>  eventSpace;
    protected Matrix piVector;
    protected Matrix infGen;
    public double runTime;
}
