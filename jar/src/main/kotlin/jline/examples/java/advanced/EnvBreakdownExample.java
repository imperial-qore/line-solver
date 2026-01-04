/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.examples.java.advanced;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.env.ENV;
import jline.solvers.fluid.FLD;
import jline.util.matrix.Matrix;
import jline.VerboseLevel;

/**
 * Example demonstrating the node breakdown/repair API for random environments.
 *
 * This example shows how to easily model a server that can break down and be
 * repaired using the new Environment convenience methods.
 */
public class EnvBreakdownExample {

    /**
     * Holder class for returning both model and queue from createBaseModel.
     */
    private static class ModelWithQueue {
        final Network model;
        final Queue queue;

        ModelWithQueue(Network model, Queue queue) {
            this.model = model;
            this.queue = queue;
        }
    }

    /**
     * Create a base queueing network model with a single server.
     *
     * @return ModelWithQueue containing the network model and queue node
     */
    private static ModelWithQueue createBaseModel() {
        Network model = new Network("ServerWithFailures");

        // Define nodes
        Source source = new Source(model, "Arrivals");
        Queue queue = new Queue(model, "Server", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Departures");

        // Define job class
        OpenClass jobclass = new OpenClass(model, "Jobs");

        // Set service and arrival rates (UP state)
        source.setArrival(jobclass, new Exp(0.8));  // Arrival rate
        queue.setService(jobclass, new Exp(2.0));   // Service rate when UP
        queue.setNumberOfServers(1);

        // Set routing
        RoutingMatrix P = model.initRoutingMatrix();
        P.addConnection(source, queue, jobclass);
        P.addConnection(queue, sink, jobclass);
        model.link(P);

        return new ModelWithQueue(model, queue);
    }

    /**
     * Example 1: Using addNodeFailureRepair convenience method with node object (recommended).
     */
    public static void example1_basicFailureRepair() {
        System.out.println("=".repeat(70));
        System.out.println("Example 1: Using addNodeFailureRepair with node object");
        System.out.println("=".repeat(70));

        ModelWithQueue modelWithQueue = createBaseModel();

        // Create environment with 2 stages (UP and DOWN)
        Environment env = new Environment("ServerEnv1", 2);

        // Add failure and repair for the server node in one call
        // Parameters: baseModel, nodeOrName, breakdownDist, repairDist, downServiceDist
        // Note: Can pass either node object (queue) or node name ("Server")
        env.addNodeFailureRepair(modelWithQueue.model, modelWithQueue.queue, new Exp(0.1), new Exp(1.0), new Exp(0.5));

        // Initialize and print stage table
        env.init();
        System.out.println("\nStage table for env1:");
        env.printStageTable();
        System.out.println();
    }

    /**
     * Example 2: Using separate breakdown and repair calls with node object.
     */
    public static void example2_separateCalls() {
        System.out.println("=".repeat(70));
        System.out.println("Example 2: Using separate breakdown and repair calls with node object");
        System.out.println("=".repeat(70));

        ModelWithQueue modelWithQueue = createBaseModel();

        Environment env = new Environment("ServerEnv2", 2);

        // Add breakdown (creates UP and DOWN stages) - using node object
        env.addNodeBreakdown(modelWithQueue.model, modelWithQueue.queue, new Exp(0.1), new Exp(0.5));

        // Add repair transition - using node object
        env.addNodeRepair(modelWithQueue.queue, new Exp(1.0));

        // Initialize and print stage table
        env.init();
        System.out.println("\nStage table for env2:");
        env.printStageTable();
        System.out.println();
    }

    /**
     * Example 3: With custom reset policies using node object.
     */
    public static void example3_customResetPolicies() {
        System.out.println("=".repeat(70));
        System.out.println("Example 3: With custom reset policies using node object");
        System.out.println("=".repeat(70));

        ModelWithQueue modelWithQueue = createBaseModel();

        Environment env = new Environment("ServerEnv3", 2);

        // Reset policy: clear all queues on breakdown
        Environment.ResetQueueLengthsFunction resetBreakdown = input -> {
            Matrix result = new Matrix(input.getNumRows(), input.getNumCols());
            result.zero();  // Clear all queues when server breaks down
            return result;
        };

        // Reset policy: keep all jobs on repair
        Environment.ResetQueueLengthsFunction resetRepair = input -> input;  // Keep jobs when repaired

        // Using node object instead of string name
        env.addNodeFailureRepair(modelWithQueue.model, modelWithQueue.queue, new Exp(0.1), new Exp(1.0), new Exp(0.5),
                               resetBreakdown, resetRepair);

        env.init();
        System.out.println("\nStage table for env3:");
        env.printStageTable();
        System.out.println();
    }

    /**
     * Example 4: Modifying reset policies after creation using node object.
     */
    public static void example4_modifyResetPolicies() {
        System.out.println("=".repeat(70));
        System.out.println("Example 4: Modifying reset policies after creation using node object");
        System.out.println("=".repeat(70));

        ModelWithQueue modelWithQueue = createBaseModel();

        Environment env = new Environment("ServerEnv4", 2);

        // Create environment with default reset policies - using node object
        env.addNodeFailureRepair(modelWithQueue.model, modelWithQueue.queue, new Exp(0.1), new Exp(1.0), new Exp(0.5));

        // Update breakdown reset policy to clear queues - using node object
        env.setBreakdownResetPolicy(modelWithQueue.queue, input -> {
            Matrix result = new Matrix(input.getNumRows(), input.getNumCols());
            result.zero();
            return result;
        });

        // Update repair reset policy (keep jobs) - using node object
        env.setRepairResetPolicy(modelWithQueue.queue, input -> input);

        env.init();
        System.out.println("\nStage table for env4:");
        env.printStageTable();
        System.out.println();
    }

    /**
     * Example 5: Solve the environment model and display results.
     */
    public static void example5_solveEnvironment() throws Exception {
        System.out.println("=".repeat(70));
        System.out.println("Example 5: Solving environment model with ENV solver");
        System.out.println("=".repeat(70));

        ModelWithQueue modelWithQueue = createBaseModel();

        Environment env = new Environment("ServerEnv5", 2);
        // Using node object
        env.addNodeFailureRepair(modelWithQueue.model, modelWithQueue.queue, new Exp(0.1), new Exp(1.0), new Exp(0.5));
        env.init();

        // Create solver options
        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 0.01;
        options.iter_max = 100;
        options.timespan[0] = 0;
        options.timespan[1] = Double.POSITIVE_INFINITY;
        options.verbose = VerboseLevel.STD;

        // Create fluid solver options for each stage
        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.setODEMaxStep(0.25);

        // Create solvers for each stage
        int numStages = env.getEnsemble().size();
        NetworkSolver[] solvers = new NetworkSolver[numStages];
        for (int e = 0; e < numStages; e++) {
            solvers[e] = new FLD(env.getModel(e));
            solvers[e].options = fluidOptions;
        }

        // Create and run ENV solver
        ENV solver = new ENV(env, solvers, options);

        System.out.println("\nAverage Performance Metrics:");
        solver.getAvgTable().print();

        System.out.println("\nInterpretation:");
        System.out.println("- The system alternates between UP (operational) and DOWN (failed) states");
        System.out.println("- UP state: Server processes jobs at rate 2.0");
        System.out.println("- DOWN state: Server processes jobs at reduced rate 0.5");
        System.out.println("- Breakdown occurs at rate 0.1 (mean time to failure = 10 time units)");
        System.out.println("- Repair occurs at rate 1.0 (mean time to repair = 1 time unit)");
        System.out.println("- Results show averaged performance across both states");
        System.out.println();
    }

    /**
     * Main method to run all examples.
     */
    public static void main(String[] args) {
        try {
            example1_basicFailureRepair();
            example2_separateCalls();
            example3_customResetPolicies();
            example4_modifyResetPolicies();
            example5_solveEnvironment();

            System.out.println("=".repeat(70));
            System.out.println("All examples completed successfully!");
            System.out.println("=".repeat(70));
        } catch (Exception e) {
            System.err.println("Error running examples:");
            e.printStackTrace();
        }
    }
}
