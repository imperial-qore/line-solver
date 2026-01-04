package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.Environment;
import jline.lang.Network;
import jline.lang.constant.SolverType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.solvers.NetworkSolver;
import jline.solvers.env.ENV;
import jline.solvers.fluid.FLD;
import jline.solvers.jmt.JMT;
import jline.solvers.mva.MVA;
import jline.solvers.SolverOptions;
import jline.VerboseLevel;
import java.util.Scanner;

/**
 * Examples demonstrating queueing networks in random environments.
 * 
 * This class provides Java implementations corresponding to the Kotlin notebooks
 * in jline.examples.kotlin.advanced.randomEnv package.
 */
public class RandomEnvExamples {

    private static final Scanner scanner = new Scanner(System.in);

    private static void pauseForUser() {
        // Skip pause if running in non-interactive mode (e.g., Maven exec)
        if (System.console() == null) {
            System.out.println("\n[Running in non-interactive mode, continuing...]");
            return;
        }
        System.out.println("\nPress Enter to continue to next example...");
        try {
            scanner.nextLine();
        } catch (Exception e) {
            // Ignore scanner errors in case of pipe or redirection
        }
    }

    /**
     * Demonstrates a random environment model with two stages (renv_twostages_repairmen.ipynb).
     * 
     * This example models a system that alternates between two environmental states,
     * such as normal operation and degraded mode. The repairmen model captures how
     * the system transitions between states and how performance varies in each state.
     * 
     * Features:
     * - Two-stage random environment
     * - State-dependent service rates
     * - Environmental state transitions
     * - Analysis of availability and performance trade-offs
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void renv_twostages_repairmen() throws Exception {
        Environment envModel = RandomEnvironmentModel.renv_twostages_repairmen();
        int E = envModel.getEnsemble().size();
        
        // Create solver options
        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 0.01;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.STD;
        
        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan[1] = 1000;
        fluidOptions.stiff = false;
        fluidOptions.setODEMaxStep(0.25);
        
        // Create solvers for each stage
        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new FLD(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }
        
        // Create environment solver
        ENV solver = new ENV(envModel, solvers, options);
        
        try {
            solver.getAvgTable().print();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates a random environment model with three stages (renv_threestages_repairmen.ipynb).
     * 
     * This example extends the two-stage model to include three environmental states,
     * allowing for more complex failure and recovery patterns. This could model systems
     * with multiple failure modes or degradation levels.
     * 
     * Features:
     * - Three-stage random environment
     * - Multiple degradation levels
     * - Complex state transition patterns
     * - Performance analysis across environmental states
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void renv_threestages_repairmen() throws Exception {
        Environment envModel = RandomEnvironmentModel.renv_threestages_repairmen();
        int E = envModel.getEnsemble().size();
        
        // Create solver options
        SolverOptions envOptions = new SolverOptions(SolverType.ENV);
        envOptions.iter_tol = 0.05;
        envOptions.timespan[0] = 0;
        
        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.stiff = false;
        fluidOptions.setODEMaxStep(0.25);
        fluidOptions.verbose = VerboseLevel.SILENT;
        
        // Create solvers for each stage
        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new FLD(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }
        
        // Create environment solver
        ENV solver = new ENV(envModel, solvers, envOptions);
        
        try {
            solver.getAvgTable().print();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates a random environment model with four stages (renv_fourstages_repairmen.ipynb).
     * 
     * This example shows a more complex random environment with four states, suitable
     * for modeling systems with multiple components that can fail independently or
     * systems with graduated performance levels based on environmental conditions.
     * 
     * Features:
     * - Four-stage random environment
     * - Rich state space for complex systems
     * - Analysis of multi-level degradation
     * - Optimization of repair strategies
     * 
     * @throws Exception if the solver encounters an error
     */
    public static void renv_fourstages_repairmen() throws Exception {
        Environment envModel = RandomEnvironmentModel.renv_fourstages_repairmen();
        int E = envModel.getEnsemble().size();
        
        // Create solver options
        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.iter_tol = 0.05;
        options.timespan[0] = 0;
        options.verbose = VerboseLevel.STD;
        
        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.stiff = false;
        fluidOptions.setODEMaxStep(0.25);
        fluidOptions.verbose = VerboseLevel.SILENT;
        
        // Create solvers for each stage
        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new FLD(envModel.getModel(e));
            solvers[e].options = fluidOptions;
        }
        
        // Create environment solver
        ENV solver = new ENV(envModel, solvers, options);
        
        try {
            solver.getAvgTable().print();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        pauseForUser();
    }

    /**
     * Demonstrates a basic random environment model (example_random_env_basics).
     *
     * This is the fundamental example for random environments, showing a simple
     * queueing system with a server that switches between two modes: Fast and Slow.
     *
     * Features:
     * - Simple closed network with delay and server
     * - Two-stage environment (Fast mode: rate 4.0, Slow mode: rate 1.0)
     * - Exponential transitions (Fast->Slow at rate 0.5, Slow->Fast at rate 1.0)
     * - Analysis using Fluid solver
     * - Environment-averaged and stage-wise performance metrics
     *
     * @return the configured environment model
     * @throws Exception if the solver encounters an error
     */
    public static Environment example_random_env_basics() throws Exception {
        // Block 1: Create base network model
        Network baseModel = new Network("BaseModel");
        Delay delay = new Delay(baseModel, "ThinkTime");
        Queue queue = new Queue(baseModel, "Fast/Slow Server", SchedStrategy.FCFS);

        // Closed class with 5 jobs
        int N = 5;
        ClosedClass jobclass = new ClosedClass(baseModel, "Jobs", N, delay, 0);
        delay.setService(jobclass, new Exp(1.0));  // Think time = 1.0
        queue.setService(jobclass, new Exp(2.0));  // Placeholder service rate

        // Connect nodes in a cycle
        baseModel.link(Network.serialRouting(delay, queue));

        // Block 2: Create the random environment
        int E = 2;
        Environment env = new Environment("ServerModes", E);

        // Stage 0: Fast mode (service rate = 4.0)
        Network fastModel = baseModel.copy();
        Queue fastQueue = (Queue) fastModel.getNodeByName("Fast/Slow Server");
        fastQueue.setService(fastModel.getClasses().get(0), new Exp(4.0));
        env.addStage(0, "Fast", "operational", fastModel);

        // Stage 1: Slow mode (service rate = 1.0)
        Network slowModel = baseModel.copy();
        Queue slowQueue = (Queue) slowModel.getNodeByName("Fast/Slow Server");
        slowQueue.setService(slowModel.getClasses().get(0), new Exp(1.0));
        env.addStage(1, "Slow", "degraded", slowModel);

        // Define transitions between stages
        env.addTransition(0, 1, new Exp(0.5));
        env.addTransition(1, 0, new Exp(1.0));

        // Block 3: Inspect the environment structure
        System.out.println("Environment stages:");
        env.printStageTable();

        // Block 4: Solve using SolverENV
        SolverOptions envOptions = new SolverOptions(SolverType.ENV);
        envOptions.iter_tol = 0.01;
        envOptions.iter_max = 50;
        envOptions.verbose = VerboseLevel.SILENT;

        SolverOptions fldOptions = new SolverOptions(SolverType.FLUID);
        fldOptions.timespan[1] = 100;
        fldOptions.verbose = VerboseLevel.SILENT;

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new FLD(env.getModel(e));
            solvers[e].options = fldOptions;
        }

        ENV envSolver = new ENV(env, solvers, envOptions);

        // Display average results weighted by environment probabilities
        System.out.println("\n--- Environment-Averaged Results ---");
        envSolver.getAvgTable().print();

        // Block 5: Compare with individual stage analysis
        System.out.println("\n--- Individual Stage Analysis (MVA) ---");
        for (int e = 0; e < E; e++) {
            System.out.println("\nStage " + e + ":");
            new MVA(env.getModel(e)).getAvgTable().print();
        }

        return env;
    }

    /**
     * Main method to run all random environment examples.
     *
     * @param args command line arguments (not used)
     */
    public static void main(String[] args) {
        System.out.println("\n=== Running example: renv_twostages_repairmen ===");
        try {
            renv_twostages_repairmen();
        } catch (Exception e) {
            System.err.println("renv_twostages_repairmen failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: renv_threestages_repairmen ===");
        try {
            renv_threestages_repairmen();
        } catch (Exception e) {
            System.err.println("renv_threestages_repairmen failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        System.out.println("\n=== Running example: renv_fourstages_repairmen ===");
        try {
            renv_fourstages_repairmen();
        } catch (Exception e) {
            System.err.println("renv_fourstages_repairmen failed: " + e.getMessage());
            e.printStackTrace();
        }
        
        scanner.close();
    }
}