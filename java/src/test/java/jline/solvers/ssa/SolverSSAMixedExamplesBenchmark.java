package jline.solvers.ssa;

import jline.examples.MixedModel;
import jline.examples.OpenModel;
import jline.lang.Network;
import jline.solvers.SolverOptions;

public class SolverSSAMixedExamplesBenchmark {

    public static void main(String[] args) {
        SolverOptions options1 = new SolverOptions();
        options1.keep = true;
        options1.cutoff = 3;
        timeSolver(MixedModel.ex1_line(), "Mixed Example 1, PS", options1);

        SolverOptions options2 = new SolverOptions();
        options2.keep = false;
        options2.cutoff = 3;
        timeSolver(MixedModel.ex2(), "Mixed Example 2, PS", options2);

        SolverOptions options3 = new SolverOptions();
        options3.keep = false;
        options3.cutoff = 3;
        timeSolver(MixedModel.ex3(), "Mixed Example 3, FCFS", options3);

        SolverOptions options5 = new SolverOptions();
        options5.keep = false;
        options5.cutoff = 3;
        options5.samples = 20000;
        timeSolver(MixedModel.ex5(), "Mixed Example 5, PS",options5);

    }





    private static void timeSolver(Network network, String modelName, SolverOptions options) {
        SolverSSA solver = new SolverSSA(network, options);
        solver.getAvgTable();
        System.out.println("Model: " + modelName);
        System.out.println("Runtime: " + solver.result.runtime);
    }


}
