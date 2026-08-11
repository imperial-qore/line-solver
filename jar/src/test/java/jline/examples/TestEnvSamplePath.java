package jline.examples;

import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.VerboseLevel;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.env.SolverENV;
import jline.solvers.fluid.SolverFluid;

/**
 * Test for Environment with contrasting utilizations using Fluid solver.
 *
 * This example demonstrates a two-stage random environment where:
 * - Fast stage: service rate = 10.0 (low utilization ~45%)
 * - Slow stage: service rate = 0.8 (high utilization ~100%, saturated)
 *
 * Corresponds to matlab/test/test_env_sample_path.m
 */
public class TestEnvSamplePath {

    public static void main(String[] args) {
        System.out.println("Creating two-stage environment model with contrasting utilizations...\n");

        // Model parameters (same as MATLAB)
        int N = 5;  // Number of jobs
        double thinkRate = 1.0;  // Think time service rate
        double fastRate = 10.0;  // Fast stage service rate
        double slowRate = 0.8;   // Slow stage service rate

        // Create environment with 2 stages
        int E = 2;
        Env envModel = new Env("ServerModes", E);

        // Create stage networks with contrasting service rates
        Network fastModel = createStageModel("FastModel", N, thinkRate, fastRate);
        Network slowModel = createStageModel("SlowModel", N, thinkRate, slowRate);

        envModel.addStage(0, "Fast", "operational", fastModel);
        envModel.addStage(1, "Slow", "degraded", slowModel);

        // Add transitions between stages (same as MATLAB)
        envModel.addTransition(0, 1, new Exp(0.5));  // Fast -> Slow
        envModel.addTransition(1, 0, new Exp(1.0));  // Slow -> Fast

        System.out.println("Environment model created with stages: Fast, Slow");
        System.out.println("Fast stage: service rate = 10.0 (low utilization)");
        System.out.println("Slow stage: service rate = 0.8 (high utilization)");
        System.out.println();

        // Analyze each stage using Fluid solver
        System.out.println("--- Individual Stage Analysis (Fluid Solver) ---");
        System.out.println("Showing steady-state metrics for each stage:\n");

        String[] stageNames = {"Fast", "Slow"};
        double[] serviceRates = {fastRate, slowRate};

        for (int e = 0; e < E; e++) {
            Network stageModel = envModel.getModel(e);

            SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
            fluidOptions.timespan[1] = 100;
            fluidOptions.verbose = VerboseLevel.SILENT;
            fluidOptions.method = "closing";

            SolverFluid solver = new SolverFluid(stageModel);
            solver.options = fluidOptions;
            solver.runAnalyzer();
            SolverResult result = solver.result;

            System.out.printf("Stage %d (%s, service rate = %.1f):%n", e, stageNames[e], serviceRates[e]);
            System.out.printf("  %-12s %-10s %-10s %-10s%n", "Station", "QLen", "Util", "Tput");
            System.out.printf("  %-12s %-10.4f %-10.4f %-10.4f%n",
                "ThinkTime", result.QN.get(0, 0), result.UN.get(0, 0), result.TN.get(0, 0));
            System.out.printf("  %-12s %-10.4f %-10.4f %-10.4f%n",
                "Server", result.QN.get(1, 0), result.UN.get(1, 0), result.TN.get(1, 0));
            System.out.println();
        }

        // Show utilization summary
        System.out.println("--- Utilization Summary ---");
        System.out.println("Fast stage (service rate 10.0): Server utilization ~45%");
        System.out.println("Slow stage (service rate 0.8): Server utilization ~100% (saturated)");
        System.out.println("\nTest completed successfully.");
    }

    /**
     * Create a closed queueing network with delay and server.
     */
    private static Network createStageModel(String name, int N, double thinkRate, double serviceRate) {
        Network model = new Network(name);

        Delay delay = new Delay(model, "ThinkTime");
        Queue queue = new Queue(model, "Server", SchedStrategy.FCFS);

        ClosedClass jobClass = new ClosedClass(model, "Jobs", N, delay, 0);

        delay.setService(jobClass, new Exp(thinkRate));
        queue.setService(jobClass, new Exp(serviceRate));

        // Serial routing: delay -> queue -> delay
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(delay, queue);
        P.set(queue, delay);
        model.link(P);

        return model;
    }
}
