package jline.solvers.jmt;

import jline.lang.Metric;
import jline.solvers.SolverResult;

import java.util.ArrayList;
import java.util.List;

public class SolverJMTResult extends SolverResult {

    // The following field is used only by the JMT solver.
    public List<Metric> metrics = new ArrayList<>();
    public double logNormConstAggr;
}
