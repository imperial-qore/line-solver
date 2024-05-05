package jline.solvers.ssa;

import jline.examples.ClosedModel;
import jline.examples.OpenModel;
import jline.lang.Network;
import jline.solvers.SolverOptions;

public class SolverSSAOpenExamplesBenchmark {

    public static void main(String[] args) {
        SolverOptions options1 = new SolverOptions();
        options1.keep = true;
        options1.cutoff = 10;
        timeSolver(OpenModel.ex1_line_v(), "Open Example 1, FCFS", options1);

        SolverOptions options3 = new SolverOptions();
        options3.keep = true;
        options3.cutoff = 7;
        timeSolver(OpenModel.ex3_line(), "Open Example 3, PS", options3);

        SolverOptions options6 = new SolverOptions();
        options6.cutoff = 1;
        timeSolver(OpenModel.ex6(), "Open Example 6, FCFS and PS", options6);
    }

    private static void timeSolver(Network network, String modelName, SolverOptions options) {
        SolverSSA solver = new SolverSSA(network, options);
        solver.getAvgTable();
        System.out.println("Model: " + modelName);
        System.out.println("Runtime: " + solver.result.runtime);
    }
}
