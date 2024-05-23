package jline.solvers.ssa;

import jline.examples.ClosedModel;
import jline.lang.Network;
import jline.lang.state.ThreadLocalRandom;
import jline.solvers.NetworkAvgTable;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.util.Maths;
import org.apache.commons.math3.random.MersenneTwister;

public class SolverSSAClosedExamplesBenchmark {


    public static void main(String[] args) {
        timeSolver(ClosedModel.ex1(), "Closed Example 1, FCFS");
        timeSolver(ClosedModel.ex2_line(), "Closed Example 2, PS");
        timeSolver(ClosedModel.ex3_line(), "Closed Example 3, PS");
        timeSolver(ClosedModel.ex4_line(), "Closed Example 4, FCFS");
        timeSolver(ClosedModel.ex7_line_ps(), "Closed Example 7, PS");
        timeSolver(ClosedModel.ex7_line_fcfs(), "Closed Example 7, FCFS");
        timeSolver(ClosedModel.ex8_line(), "Closed Example 8, FCFS");

        SolverOptions solverOptions = new SolverOptions();
        solverOptions.samples = 5000;
        timeSolver(ClosedModel.ex9_line(), "Closed Example 9, FCFS", solverOptions);

    }

    private static void timeSolver(Network network, String modelName) {
        SolverSSA solver = new SolverSSA(network);
        solver.getAvgTable();
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
