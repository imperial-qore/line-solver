package jline.solvers.ssa;

import jline.examples.ClosedModel;
import jline.lang.Network;
import jline.lang.constant.VerboseLevel;
import jline.solvers.SolverOptions;
import jline.util.Maths;

public class SolverSSAClosedExamplesBenchmark {


    public static void main(String[] args) {
        timeSolver(ClosedModel.example_closedModel_1(), "Closed Example 1, FCFS");
        timeSolver(ClosedModel.example_closedModel_2(), "Closed Example 2, PS");
        timeSolver(ClosedModel.example_closedModel_3(), "Closed Example 3, PS");
        timeSolver(ClosedModel.example_closedModel_4(), "Closed Example 4, FCFS");
        timeSolver(ClosedModel.example_closedModel_7ps(), "Closed Example 7, PS");
        timeSolver(ClosedModel.example_closedModel_7fcfs(), "Closed Example 7, FCFS");
        timeSolver(ClosedModel.example_closedModel_8(), "Closed Example 8, FCFS");

        SolverOptions solverOptions = new SolverOptions();
        solverOptions.samples = 5000;
        timeSolver(ClosedModel.example_closedModel_9(), "Closed Example 9, FCFS", solverOptions);

    }

    private static void timeSolver(Network network, String modelName) {
        SolverOptions options = new SolverOptions();
        options.method= "para";
        options.verbose = VerboseLevel.DEBUG;
        options.seed = 23000;
        Maths.setRandomNumbersMatlab(true);
        SolverSSA solver = new SolverSSA(network, options);
        solver.getAvgTable().print();
        System.out.println("Model: " + modelName);
        System.out.println("Runtime: " + solver.result.runtime);
    }

    private static void timeSolver(Network network, String modelName, SolverOptions solverOptions) {
        SolverSSA solver = new SolverSSA(network, solverOptions);
        solver.getAvgTable();
        System.out.println("Model: " + modelName);
        System.out.println("Runtime: " + solver.result.runtime);
    }


}
