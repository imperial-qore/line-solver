package jline.examples.java.advanced;

import jline.api.mc.Ctmc_solve;
import jline.io.MAPQN2RENV;
import jline.lang.ClosedClass;
import jline.lang.Environment;
import jline.lang.Model;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SolverType;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.env.ENV;
import jline.solvers.env.SolverENV;
import jline.solvers.ctmc.CTMC;
import jline.solvers.ctmc.ResultCTMC;
import jline.solvers.ctmc.handlers.Ctmc_avg_from_pi;
import jline.solvers.ctmc.handlers.Solver_ctmc;
import jline.solvers.fluid.FLD;
import jline.solvers.ln.LNTranAvgResult;
import jline.solvers.ln.SolverFactory;
import jline.solvers.ln.SolverLN;
import jline.solvers.mam.MAM;
import jline.solvers.wrappers.jmt.JMT;
import jline.solvers.mva.MVA;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.VerboseLevel;
import java.util.Scanner;

/**
 * Examples demonstrating queueing networks in random environments.
 * 
 * This class provides Java implementations corresponding to the example notebooks
 * in jline.examples.java.advanced.randomEnv package.
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
     * Solver options that will print their table.
     *
     * <p>A solver built with {@code VerboseLevel.SILENT} writes that level into
     * {@code GlobalConstants}, and a default-constructed {@code SolverOptions}
     * reads it back, so an unqualified {@code getAvgTable().print()} after a
     * silent ENV solve emits NOTHING. The per-stage tables below are part of
     * what the reference shows, so they ask for {@code STD} by name.
     *
     * @param type the solver the options are for
     * @return options at {@code VerboseLevel.STD}
     */
    private static SolverOptions printing(SolverType type) {
        SolverOptions o = new SolverOptions(type);
        o.verbose = VerboseLevel.STD;
        return o;
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
        
        // THE ITERATION CONTROLS AND THE STAGE HORIZON ARE PART OF THE GOLDEN.
        // An environment couples TRANSIENT stage solves, so a row given a
        // different horizon or integrator answers a different question: the
        // reference states timespan [0, 1e3] on the fluid solver and leaves the
        // integrator at its default, and overriding `stiff` or the ODE max step
        // walks a different trajectory.
        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.timespan = new double[]{0, Double.POSITIVE_INFINITY};
        options.iter_max = 100;
        options.iter_tol = 0.01;
        options.method = "default";
        options.verbose = VerboseLevel.STD;
        
        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan = new double[]{0, 1000};
        fluidOptions.verbose = VerboseLevel.SILENT;
        
        // Create solvers for each stage
        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new FLD(envModel.getModel(e));
            solvers[e].options = fluidOptions.copy();
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
        
        // CTMC STAGES, as the reference pins: each stage is solved exactly and
        // the horizon is the whole line, so the stage solve is a steady state.
        SolverOptions envOptions = new SolverOptions(SolverType.ENV);
        envOptions.timespan = new double[]{0, Double.POSITIVE_INFINITY};
        envOptions.iter_max = 100;
        envOptions.iter_tol = 0.05;
        envOptions.method = "default";
        
        // The stage horizon stays at the CTMC default, which is open at both
        // ends: SolverENV reads timespan[1] to choose between a steady-state and
        // a transient stage solve, and an open horizon is the steady state the
        // reference asks each stage for.
        SolverOptions ctmcOptions = new SolverOptions(SolverType.CTMC);
        ctmcOptions.stiff = false;
        ctmcOptions.verbose = VerboseLevel.SILENT;
        
        // Create solvers for each stage
        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new CTMC(envModel.getModel(e));
            solvers[e].options = ctmcOptions.copy();
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
        
        // THIS FIXED POINT IS TOLERANCE-DEPENDENT: run to convergence it lands
        // on Queue1 Tput 0.97136 where the reference's iter_tol = 0.05 stops it
        // at 0.9716, and the golden holds the latter. So the controls are stated
        // rather than left at the engine default.
        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.timespan = new double[]{0, Double.POSITIVE_INFINITY};
        options.iter_max = 100;
        options.iter_tol = 0.05;
        options.method = "default";
        options.verbose = VerboseLevel.STD;
        
        SolverOptions fluidOptions = new SolverOptions(SolverType.FLUID);
        fluidOptions.timespan = new double[]{0, Double.POSITIVE_INFINITY};
        fluidOptions.verbose = VerboseLevel.SILENT;
        
        // Create solvers for each stage
        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new FLD(envModel.getModel(e));
            solvers[e].options = fluidOptions.copy();
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
     * Demonstrates a basic random environment model (renv_basic).
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
    public static Environment renv_basic() throws Exception {
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
            solvers[e].options = fldOptions.copy();
        }

        ENV envSolver = new ENV(env, solvers, envOptions);

        // Display average results weighted by environment probabilities
        System.out.println("\n--- Environment-Averaged Results ---");
        envSolver.getAvgTable().print();

        // Block 5: Compare with individual stage analysis
        System.out.println("\n--- Individual Stage Analysis (MVA) ---");
        for (int e = 0; e < E; e++) {
            System.out.printf("%nStage: %s (prob = %.4f)%n",
                              env.getStageName(e), env.probEnv.get(0, e));
            new MVA(env.getModel(e), printing(SolverType.MVA)).getAvgTable().print();
        }

        return env;
    }


    /**
     * `renv_genqn.m`: the stage generator of the repairmen environments, solved
     * on its own by MVA.
     *
     * <p>The reference is a FUNCTION rather than a script, and it is what every
     * environment above hands to `addStage`. Publishing it under its own name
     * makes the two-station closed network the environments are built from
     * visible, and solving it once by MVA states the steady state each stage
     * would reach if the environment never switched.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void renv_genqn() throws Exception {
        Network model = RandomEnvironmentModel.renv_genqn(1.0, 0.5, 4);
        new MVA(model).getAvgTable().print();
        pauseForUser();
    }

    /**
     * `renv_map_fallback.m`: the random-environment fallback a solver takes
     * when it cannot consume a non-renewal service process.
     *
     * <p>A solver with no MAP/MMPP support does not reject the model: the
     * NetworkSolver base class intercepts it, replaces every modulating chain
     * by a set of random-environment stages in which the process is exponential
     * with its phase-conditional intensity (`MAPQN2RENV.map2renv`), and solves
     * the stages with the same solver through SolverENV. The interception is
     * METHOD-AWARE, so a method that handles the process natively (MVA 'rqna',
     * MAM, CTMC, FLD, SSA, JMT, LDES) runs unchanged.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void renv_map_fallback() throws Exception {
        Network model = RandomEnvironmentModel.mmppClosed("mmppClosed");

        // The environment image the solver builds for itself, shown here.
        SolverOptions imgOptions = new SolverOptions(SolverType.MVA);
        MAPQN2RENV.RenvImage image = MAPQN2RENV.map2renvImage(model, imgOptions);
        System.out.printf("Random-environment image: %d stages, phase orders %s, MMPP image: %b%n",
                image.nstages, java.util.Arrays.toString(image.orders), image.isMMPP);
        image.env.printStageTable();

        // MVA reaches the answer THROUGH the fallback; CTMC takes the MMPP
        // natively and is exact, so it is the reference the fallback is read
        // against.
        System.out.println("\n--- MVA (through the random-environment fallback) ---");
        new MVA(model).getAvgTable().print();
        System.out.println("\n--- CTMC (native MMPP, exact) ---");
        new CTMC(model).getAvgTable().print();

        // The two environment LIMITS, forced. 'dec' is the quasi-stationary
        // limit of a SLOW environment (each phase reaches its own steady state
        // before the next switch); 'avg' is the rate-averaged limit of a FAST
        // one (the phases blur into their mean rate).
        for (String limit : new String[]{"dec", "avg"}) {
            SolverOptions options = new SolverOptions(SolverType.MVA);
            options.config.map_env_method = limit;
            System.out.println("\n--- MVA, map_env_method = '" + limit + "' ---");
            new MVA(model, options).getAvgTable().print();
        }

        // Opting out restores the feature rejection, which is the point: the
        // fallback is a convenience, not a claim that MVA understands a MAP.
        SolverOptions off = new SolverOptions(SolverType.MVA);
        off.config.map_env = "off";
        try {
            new MVA(model, off).getAvgTable();
            System.out.println("map_env='off': the model was accepted, which it should not be");
        } catch (Exception ex) {
            System.out.println("\nmap_env='off': " + ex.getMessage());
        }
        pauseForUser();
    }

    /**
     * `example_mapqn2renv.m`: the same conversion reached through MAPQN2RENV,
     * the named entry point, with the stage networks then solved one at a time.
     *
     * <p>`mapqn2renv` is a pure delegate to `map2renv`; it is kept because the
     * reference names it, and because the per-stage solve below is the part
     * that shows what the image is FOR. Each stage is an ordinary closed
     * network that any solver can take: MVA cannot take the MMPP, but it can
     * take these.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void example_mapqn2renv() throws Exception {
        System.out.println("=== MMPP2 Closed QN to Random Environment Transformation ===\n");
        Network model = RandomEnvironmentModel.mmppClosed("MMPP_ClosedQN");
        System.out.println("Original Network:");
        System.out.println("  Topology: Delay -> MMPP2 Queue -> Delay (cyclic)");
        System.out.println("  Population: N = 5 jobs");
        System.out.println("  Think time: Exp(1.0)");
        System.out.println("  Service: MMPP2 with phases 0,1\n");

        Environment envModel = MAPQN2RENV.mapqn2renv(model);
        System.out.println("Environment stages:");
        envModel.printStageTable();

        int E = envModel.getNumberOfStages();

        // The fluid coupling, over a finite stage horizon: the mean-field exit
        // metrics are a Riemann-Stieltjes sum of the stage trajectory against
        // the holding-time CDF, so the horizon is part of the answer.
        SolverOptions envOptions = new SolverOptions(SolverType.ENV);
        envOptions.iter_max = 100;
        envOptions.iter_tol = 0.01;
        envOptions.method = "default";
        envOptions.verbose = VerboseLevel.SILENT;

        SolverOptions fldOptions = new SolverOptions(SolverType.FLUID);
        fldOptions.timespan = new double[]{0, 100};
        fldOptions.verbose = VerboseLevel.SILENT;

        // The reference opens with the ORIGINAL MMPP2 model solved by JMT, the
        // simulation the two environment solves are read against.
        try {
            System.out.println("\nSolving original MMPP2 model using JMT...");
            System.out.println("Original MMPP2 Queue results:");
            new JMT(model).getAvgTable().print();
        } catch (Exception ex) {
            System.out.println("JMT solver error: " + ex.getMessage());
        }

        System.out.println("\nSolving environment model using ENV (with FLD)...");
        try {
            NetworkSolver[] fldSolvers = new NetworkSolver[E];
            for (int e = 0; e < E; e++) {
                fldSolvers[e] = new FLD(envModel.getModel(e));
                fldSolvers[e].options = fldOptions.copy();
            }
            new ENV(envModel, fldSolvers, envOptions.copy()).getAvgTable().print();
        } catch (Exception ex) {
            System.out.println("ENV error: " + ex.getMessage() + "\n");
            System.out.println("Note: ENV requires transient analysis support.");
            System.out.println("The environment model was created successfully.\n");
        }

        // The CTMC coupling reads the SAME stage models the fluid one just
        // wrote a purposely FRACTIONAL initial state onto: a fluid initial
        // condition is a mean, not a state. A discrete stage solver can only
        // start from a state its enumeration contains, so this arm is guarded
        // and its refusal reported, exactly as the reference guards it.
        System.out.println("\nSolving environment model using ENV (with CTMC, cutoff=100)...");
        try {
            SolverOptions ctmcOptions = new SolverOptions(SolverType.CTMC);
            ctmcOptions.method = "exact";
            ctmcOptions.timespan = new double[]{0, 100};
            ctmcOptions.cutoff = Matrix.singleton(100);
            ctmcOptions.verbose = VerboseLevel.SILENT;
            NetworkSolver[] ctmcSolvers = new NetworkSolver[E];
            for (int e = 0; e < E; e++) {
                ctmcSolvers[e] = new CTMC(envModel.getModel(e), ctmcOptions.copy());
            }
            new ENV(envModel, ctmcSolvers, envOptions.copy()).getAvgTable().print();
        } catch (Exception ex) {
            System.out.println("ENV/CTMC error: " + ex.getMessage() + "\n");
        }

        System.out.println("\nSolving stage networks individually (steady-state with MVA)...");
        try {
            for (int e = 0; e < E; e++) {
                System.out.println("\n  Stage " + envModel.getStageName(e) + ":");
                new MVA(envModel.getModel(e), printing(SolverType.MVA)).getAvgTable().print();
            }
        } catch (Exception ex) {
            System.out.println("Stage solver error: " + ex.getMessage());
        }
        System.out.println("\n=== Transformation Complete ===");
        pauseForUser();
    }

    /** The three metrics the exact joint chain is read for. */
    private static final class JointAvg {
        Matrix QN;
        Matrix UN;
        Matrix TN;
    }

    /**
     * Day-averaged metrics from the exact joint (hour x network-state) CTMC.
     *
     * <p>`SolverENV.getGenerator` returns the FLATTENED generator of the
     * environment: each stage's own generator on its diagonal block and the
     * environment's rate on the identity of the off-diagonal one, since a stage
     * switch carries the network state across unchanged. Solving that chain is
     * the ground truth both couplings are read against.
     *
     * <p>The per-stage metrics come from `ctmc_avg_from_pi` on the CONDITIONAL
     * law of each block. One consequence is worth naming: that function reports
     * `max` of the arrival- and departure-based utilization estimates, and the
     * two coincide only in equilibrium. A stage conditional law is not one, so
     * the Util column sits slightly above `T*E[S]/c` while QLen and Tput are
     * exact.
     */
    private static JointAvg exactJointMetrics(Environment env, double horizon) {
        int E = env.getNumberOfStages();
        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new CTMC(env.getModel(e), ctmcStageOptions(horizon));
        }
        SolverENV solverX = new SolverENV(env, solvers);
        Matrix piJoint = Ctmc_solve.ctmc_solve(solverX.getGenerator().renvInfGen);

        Matrix QN = null;
        Matrix UN = null;
        Matrix TN = null;
        int off = 0;
        for (int e = 0; e < E; e++) {
            NetworkStruct sn = env.getModel(e).getStruct();
            ResultCTMC r = Solver_ctmc.solver_ctmc(sn, ctmcStageOptions(horizon));
            int ns = r.getQ().getNumRows();
            Matrix blk = new Matrix(1, ns);
            double prob = 0;
            for (int i = 0; i < ns; i++) {
                double v = piJoint.get(off + i);
                blk.set(0, i, v);
                prob += v;
            }
            off += ns;
            if (prob > 0) {
                for (int i = 0; i < ns; i++) {
                    blk.set(0, i, blk.get(0, i) / prob);
                }
            }
            Ctmc_avg_from_pi.Result a = Ctmc_avg_from_pi.ctmc_avg_from_pi(
                    r.getSn(), blk, r.getStateSpace(), r.getStateSpaceAggr(),
                    r.getArvRates(), r.getDepRates());
            if (QN == null) {
                QN = new Matrix(a.QN.getNumRows(), a.QN.getNumCols());
                UN = new Matrix(a.QN.getNumRows(), a.QN.getNumCols());
                TN = new Matrix(a.QN.getNumRows(), a.QN.getNumCols());
            }
            for (int i = 0; i < QN.getNumRows(); i++) {
                for (int k = 0; k < QN.getNumCols(); k++) {
                    QN.set(i, k, QN.get(i, k) + prob * a.QN.get(i, k));
                    UN.set(i, k, UN.get(i, k) + prob * a.UN.get(i, k));
                    TN.set(i, k, TN.get(i, k) + prob * a.TN.get(i, k));
                }
            }
        }
        JointAvg out = new JointAvg();
        out.QN = QN;
        out.UN = UN;
        out.TN = TN;
        return out;
    }

    /** The stage solver of the terminal example: an exact CTMC over a finite horizon. */
    private static SolverOptions ctmcStageOptions(double horizon) {
        SolverOptions o = new SolverOptions(SolverType.CTMC);
        o.method = "exact";
        o.timespan = new double[]{0, horizon};
        o.verbose = VerboseLevel.SILENT;
        return o;
    }

    /** One station's metric, totalled over the classes it carries. */
    private static double stationSum(Matrix A, int st) {
        double s = 0;
        for (int k = 0; k < A.getNumCols(); k++) {
            s += A.get(st, k);
        }
        return s;
    }

    /** One row of the analyzer comparison the reference prints. */
    private static void compareRow(String label, Matrix Q, Matrix U, Matrix T, int st) {
        System.out.printf("%-10s %12.5f %12.5f %12.5f%n", label, stationSum(Q, st),
                stationSum(U, st), stationSum(T, st));
    }

    /** The 24 hourly stages wired into the fixed daily cycle 1 -> 2 -> ... -> 24 -> 1. */
    private static Environment dailyCycle(String name, Network[] stage, double[] durationHr) {
        int E = stage.length;
        Environment env = new Environment(name, E);
        for (int h = 0; h < E; h++) {
            env.addStage(h, String.format("Hour%02d", h), "operational", stage[h]);
        }
        for (int h = 0; h < E; h++) {
            env.addTransition(h, (h + 1) % E, new Exp(1.0 / durationHr[h]));
        }
        env.init();
        return env;
    }

    /**
     * `renv_container_terminal.m`: a daily-cycle container terminal, solved by
     * the STATE-VECTOR analyzer and read against the exact joint chain.
     *
     * <p>The Rotterdam terminal of Dhingra et al.: container handling demand
     * varies over the 24 hours of a day, each hour is one environment stage
     * with its own demand intensity and random (exponential) duration, and the
     * environment visits the stages in a fixed daily cycle. The semi-open SOQN
     * of the original study is rendered as the equivalent finite-token CLOSED
     * network that an environment stage must be: N straddle carriers cycle
     * between a yard staging Delay and a multi-server quay-crane Queue, and the
     * hourly demand modulates the yard staging rate, so the cranes congest
     * during the peaks.
     *
     * <p>WHAT THE THREE ROWS MEASURE. `statevec` carries the whole joint
     * distribution across a stage switch, `meanfield` collapses it to marginal
     * mean queue lengths, and `exact` is the stationary law of the full (hour x
     * network-state) chain. The state-vector blend reproduces the exact answer
     * to nine digits here, which is the property the example exists to show;
     * the mean-field collapse does not.
     *
     * <p>The second section repeats it with the internal handling (quay cranes
     * -> stacking cranes) collapsed by Norton's theorem into one closed
     * load-dependent FES, so `sn.lldscaling` reaches the state-vector
     * analyzer's CTMC backend. The third analyses the OPEN counterpart, where
     * the hourly demand is an MMPP(24) into a multi-server queue with an
     * unbounded buffer: that is outside the closed-CTMC state-vector path and
     * is what MAM solves exactly as a QBD.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void renv_container_terminal() throws Exception {
        double[] dailyRateHr = {6, 30, 40, 62, 76, 79, 119, 164, 152, 130, 79, 70,
                                57, 57, 113, 130, 162, 202, 148, 118, 92, 62, 36, 8};
        double[] dailyDurationHr = {0.51, 0.84, 0.76, 0.98, 0.67, 2.33, 1.80, 0.62,
                                    0.36, 1.30, 1.02, 0.98, 0.86, 1.56, 1.30, 1.57,
                                    0.34, 0.22, 0.47, 1.33, 1.56, 0.93, 0.71, 1.11};
        int E = dailyRateHr.length;   // hourly environment stages
        int N = 6;                    // straddle carriers circulating (the closed tokens)
        int nCranes = 2;              // quay cranes, a multi-server FCFS station
        double craneRate = 16;        // moves per hour served by one crane

        Network[] stage = new Network[E];
        for (int h = 0; h < E; h++) {
            // The per-token yard completion rate of hour h: the closed-network
            // image of that hour's external demand intensity.
            stage[h] = RandomEnvironmentModel.terminalModel(dailyRateHr[h] / N, craneRate, nCranes, N);
        }
        Environment env = dailyCycle("RotterdamDailyCycle", stage, dailyDurationHr);
        System.out.printf("Rotterdam container terminal: %d hourly stages, %d carriers, %d cranes.%n",
                E, N, nCranes);

        // The stage solver is an exact CTMC over a 12-hour transient horizon,
        // which covers the hour-duration CDF past 99.9% for every stage, so the
        // quadrature is not truncating the blend.
        double T = 12;
        SolverOptions statevecOpt = envOptions("statevec");
        SolverOptions meanfieldOpt = envOptions("meanfield");

        SolverENV sv = envOverCtmcStages(env, T, statevecOpt);
        SolverENV mf = envOverCtmcStages(env, T, meanfieldOpt);
        sv.getEnsembleAvg();
        mf.getEnsembleAvg();
        Matrix Qsv = sv.result.QN;
        Matrix Usv = sv.result.UN;
        Matrix Tsv = sv.result.TN;
        Matrix Qmf = mf.result.QN;
        Matrix Umf = mf.result.UN;
        Matrix Tmf = mf.result.TN;
        JointAvg ex = exactJointMetrics(env, T);

        int crane = 1;   // the QuayCranes station
        System.out.println("\n=== Day-averaged quay-crane metrics (container terminal) ===");
        System.out.printf("%-10s %12s %12s %12s%n", "analyzer", "QLen", "Util", "Tput");
        compareRow("exact", ex.QN, ex.UN, ex.TN, crane);
        compareRow("statevec", Qsv, Usv, Tsv, crane);
        compareRow("meanfield", Qmf, Umf, Tmf, crane);
        System.out.printf("%nQLen error vs exact:  statevec = %.3e , meanfield = %.3e%n",
                Math.abs(stationSum(Qsv, crane) - stationSum(ex.QN, crane)),
                Math.abs(stationSum(Qmf, crane) - stationSum(ex.QN, crane)));
        System.out.printf("Util error vs exact:  statevec = %.3e , meanfield = %.3e%n",
                Math.abs(stationSum(Usv, crane) - stationSum(ex.UN, crane)),
                Math.abs(stationSum(Umf, crane) - stationSum(ex.UN, crane)));

        System.out.println("\nDay-averaged crane-queue table (state-vector analyzer):");
        sv.getAvgTable().print();

        System.out.println("\n=== Closed FES variant (Norton-aggregated terminal internals) ===");
        Matrix fesRate = RandomEnvironmentModel.fesRateCurve(16, 20, N);
        StringBuilder curve = new StringBuilder("Norton FES rate curve mu(n) = [");
        for (int i = 0; i < fesRate.length(); i++) {
            curve.append(i > 0 ? " " : "").append(String.format("%.5f", fesRate.get(i)));
        }
        System.out.println(curve.append(']').toString());

        Network[] stageFES = new Network[E];
        for (int h = 0; h < E; h++) {
            stageFES[h] = RandomEnvironmentModel.terminalFESModel(dailyRateHr[h] / N, fesRate);
        }
        Environment envFES = dailyCycle("RotterdamDailyCycleFES", stageFES, dailyDurationHr);

        SolverENV svf = envOverCtmcStages(envFES, T, statevecOpt);
        SolverENV mff = envOverCtmcStages(envFES, T, meanfieldOpt);
        svf.getEnsembleAvg();
        mff.getEnsembleAvg();
        Matrix Qsf = svf.result.QN;
        Matrix Usf = svf.result.UN;
        Matrix Tsf = svf.result.TN;
        Matrix Qff = mff.result.QN;
        Matrix Uff = mff.result.UN;
        Matrix Tff = mff.result.TN;
        JointAvg exf = exactJointMetrics(envFES, T);

        int fes = 1;   // the load-dependent FES station
        System.out.printf("%-10s %12s %12s %12s%n", "analyzer", "QLen", "Util", "Tput");
        compareRow("exact", exf.QN, exf.UN, exf.TN, fes);
        compareRow("statevec", Qsf, Usf, Tsf, fes);
        compareRow("meanfield", Qff, Uff, Tff, fes);
        System.out.printf("%nQLen error vs exact:  statevec = %.3e , meanfield = %.3e%n",
                Math.abs(stationSum(Qsf, fes) - stationSum(exf.QN, fes)),
                Math.abs(stationSum(Qff, fes) - stationSum(exf.QN, fes)));
        System.out.printf("Util error vs exact:  statevec = %.3e , meanfield = %.3e%n",
                Math.abs(stationSum(Usf, fes) - stationSum(exf.UN, fes)),
                Math.abs(stationSum(Uff, fes) - stationSum(exf.UN, fes)));

        // The open counterpart: the 24-hour demand as a Markov-modulated
        // Poisson stream rather than a finite carrier pool. The ring of hours
        // IS the modulating generator and the hourly demand its per-phase
        // arrival rate, so the model is an MMPP(24)/M/c with an unbounded
        // buffer -- outside the closed-CTMC state-vector path, and exactly the
        // QBD that MAM takes.
        System.out.println("\n=== Open system via MAM (MMPP(24)/M/c, infinite buffer) ===");
        Matrix D0 = new Matrix(E, E);
        Matrix D1 = new Matrix(E, E);
        for (int h = 0; h < E; h++) {
            D0.set(h, h, -1.0 / dailyDurationHr[h] - dailyRateHr[h]);
            D0.set(h, (h + 1) % E, 1.0 / dailyDurationHr[h]);
            D1.set(h, h, dailyRateHr[h]);
        }
        int cOpen = 8;   // an open stream is not throttled by a token pool
        Network openModel = new Network("OpenTerminal");
        Source ships = new Source(openModel, "Ships");
        Queue quayOpen = new Queue(openModel, "QuayCranes", SchedStrategy.FCFS);
        Sink done = new Sink(openModel, "Departures");
        OpenClass oc = new OpenClass(openModel, "Containers");
        ships.setArrival(oc, new MAP(D0, D1));
        quayOpen.setService(oc, new Exp(craneRate));
        quayOpen.setNumberOfServers(cOpen);
        openModel.link(Network.serialRouting(ships, quayOpen, done));

        double num = 0;
        double den = 0;
        for (int h = 0; h < E; h++) {
            num += dailyDurationHr[h] * dailyRateHr[h];
            den += dailyDurationHr[h];
        }
        double meanLambda = num / den;
        System.out.printf("MMPP mean arrival = %.2f/hr, %d cranes @ %.0f/hr, rho = %.3f%n",
                meanLambda, cOpen, craneRate, meanLambda / (cOpen * craneRate));
        new MAM(openModel).getAvgTable().print();
        pauseForUser();
    }

    /** The ENV options the terminal example drives both couplings with. */
    private static SolverOptions envOptions(String method) {
        SolverOptions o = new SolverOptions(SolverType.ENV);
        o.iter_max = 100;
        o.iter_tol = 1e-5;
        o.method = method;
        o.verbose = VerboseLevel.SILENT;
        return o;
    }

    /** An ENV solver over per-stage exact CTMC solvers with a finite horizon. */
    private static SolverENV envOverCtmcStages(Environment env, double horizon, SolverOptions envOpt) {
        int E = env.getNumberOfStages();
        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            solvers[e] = new CTMC(env.getModel(e), ctmcStageOptions(horizon));
        }
        return new SolverENV(env, solvers, envOpt.copy());
    }

    /** `PoissonSOQNSpec`: the geometry of the Rotterdam semi-open network. */
    private static final class SoqnParams {
        int numEntryServers = 6;
        int numStacks = 29;
        int numExitServers = 6;
        double entryServiceTime = 6.0;
        double travelToStackTime = 5.6;
        double stackServiceTime = 6.0;
        double travelToExitTime = 5.6;
        double exitServiceTime = 6.0;
    }

    /**
     * The flow-equivalent throughput curve of one SOQN subnetwork at l = 1..N
     * tokens, by exact MVA on the underlying closed network.
     *
     * <p>`upstream` is S1, entry gates -> travel -> 29 stacks -> back; the
     * downstream S2 is travel-to-exit against exit gates. `mu[0]` is left at
     * zero so the curve is indexed by occupancy directly, which is how the QBD
     * blocks read it.
     */
    private static double[] fesCurve(int N, SoqnParams p, boolean upstream) {
        double[] mu = new double[N + 1];
        for (int l = 1; l <= N; l++) {
            Network m = new Network(upstream ? "S1" : "S2");
            int gate;
            if (upstream) {
                Queue eg = new Queue(m, "EntryGates", SchedStrategy.FCFS);
                eg.setNumberOfServers(p.numEntryServers);
                Delay tv = new Delay(m, "TravelToStack");
                Queue[] st = new Queue[p.numStacks];
                for (int i = 0; i < p.numStacks; i++) {
                    st[i] = new Queue(m, "Stack" + (i + 1), SchedStrategy.FCFS);
                    st[i].setNumberOfServers(1);
                }
                ClosedClass cls = new ClosedClass(m, "Trucks", l, eg, 0);
                eg.setService(cls, new Exp(1.0 / p.entryServiceTime));
                tv.setService(cls, new Exp(1.0 / p.travelToStackTime));
                for (int i = 0; i < p.numStacks; i++) {
                    st[i].setService(cls, new Exp(1.0 / p.stackServiceTime));
                }
                RoutingMatrix P = m.initRoutingMatrix();
                P.set(cls, cls, eg, tv, 1.0);
                for (int i = 0; i < p.numStacks; i++) {
                    P.set(cls, cls, tv, st[i], 1.0 / p.numStacks);
                    P.set(cls, cls, st[i], eg, 1.0);
                }
                m.link(P);
                gate = 0;                      // EntryGates
            } else {
                Delay tv = new Delay(m, "TravelToExit");
                Queue xg = new Queue(m, "ExitGates", SchedStrategy.FCFS);
                xg.setNumberOfServers(p.numExitServers);
                ClosedClass cls = new ClosedClass(m, "Trucks", l, tv, 0);
                tv.setService(cls, new Exp(1.0 / p.travelToExitTime));
                xg.setService(cls, new Exp(1.0 / p.exitServiceTime));
                RoutingMatrix P = m.initRoutingMatrix();
                P.set(cls, cls, tv, xg, 1.0);
                P.set(cls, cls, xg, tv, 1.0);
                m.link(P);
                gate = 1;                      // ExitGates
            }
            SolverOptions mvaOpt = new SolverOptions(SolverType.MVA);
            mvaOpt.method = "exact";
            mvaOpt.verbose = VerboseLevel.SILENT;
            // The gate station's throughput IS the subnetwork's completion rate
            // at this occupancy, since every token passes it once per cycle.
            mu[l] = new MVA(m, mvaOpt).getAvgTput().get(gate, 0);
        }
        return mu;
    }

    /**
     * A banded LU factorization with NO PIVOTING, in column-major band storage.
     *
     * <p>The unpivoted form is what makes the blend affordable: every one of the
     * 24 resolvents is factored ONCE and then applied on every sweep of the
     * fixed point, so the iteration costs band solves rather than 24 fresh
     * factorizations per sweep. Skipping the pivot search is safe here and not
     * a shortcut: the matrix is `(sI - Q)^T` with `s > 0` and `Q` a generator,
     * so `sI - Q` is strictly row diagonally dominant and its transpose is
     * strictly column diagonally dominant, which is exactly the condition under
     * which Gaussian elimination without pivoting is stable.
     */
    private static final class Banded {
        private final int n;
        private final int kl;
        private final int ku;
        private final int ld;
        private final double[] a;

        Banded(int n, int kl, int ku) {
            this.n = n;
            this.kl = kl;
            this.ku = ku;
            this.ld = kl + ku + 1;
            this.a = new double[(kl + ku + 1) * n];
        }

        double at(int i, int j) {
            return a[j * ld + (i + ku - j)];
        }

        void put(int i, int j, double v) {
            a[j * ld + (i + ku - j)] = v;
        }

        void factor() {
            for (int j = 0; j + 1 < n; j++) {
                double d = at(j, j);
                int imax = Math.min(n - 1, j + kl);
                int kmax = Math.min(n - 1, j + ku);
                for (int i = j + 1; i <= imax; i++) {
                    double l = at(i, j) / d;
                    put(i, j, l);
                    if (l == 0.0) {
                        continue;
                    }
                    for (int k = j + 1; k <= kmax; k++) {
                        put(i, k, at(i, k) - l * at(j, k));
                    }
                }
            }
        }

        void solve(double[] b) {
            for (int j = 0; j + 1 < n; j++) {
                int imax = Math.min(n - 1, j + kl);
                for (int i = j + 1; i <= imax; i++) {
                    b[i] -= at(i, j) * b[j];
                }
            }
            for (int jj = n - 1; jj >= 0; jj--) {
                int kmax = Math.min(n - 1, jj + ku);
                for (int k = jj + 1; k <= kmax; k++) {
                    b[jj] -= at(jj, k) * b[k];
                }
                b[jj] /= at(jj, jj);
            }
        }
    }

    /**
     * `(sI - Q)^T` of the level-dependent QBD, banded and ready to factor.
     *
     * <p>The state is `(n, k)`: level `n` counts the jobs upstream of S2 (those
     * inside S1 plus the external backlog) and phase `k` the S2 occupancy, so
     * an S1 completion moves `(n,k) -> (n-1,k+1)`, an S2 completion `(n,k) ->
     * (n,k-1)` and an arrival `(n,k) -> (n+1,k)`. S1 serves at `mu1(min(n,
     * N-k))`, since a token can only be inside S1 if one of the pool is free to
     * hold it, and the top level `Mtr` is the truncation, where arrivals are
     * dropped.
     *
     * <p>Everything is written TRANSPOSED, because the resolvent is applied to a
     * ROW vector `pi` and the band solve wants a column system.
     */
    private static Banded soqnResolventMatrix(int N, double lam, double[] mu1, double[] mu2,
                                              int tailFactor, double s) {
        int M = N + 1;
        int Mtr = N + tailFactor * N;
        int dim = (Mtr + 1) * M;
        Banded B = new Banded(dim, M, M - 1);
        for (int n = 0; n <= Mtr; n++) {
            double lamEff = (n == Mtr) ? 0.0 : lam;
            for (int k = 0; k <= N; k++) {
                int row = n * M + k;
                double m1 = mu1[Math.min(n, N - k)];
                double m2 = mu2[k];
                B.put(row, row, s + lamEff + m1 + m2);
                if (k >= 1) {
                    B.put(n * M + (k - 1), row, -m2);
                }
                if (n < Mtr) {
                    B.put((n + 1) * M + k, row, -lam);
                }
                if (n >= 1 && k + 1 <= N) {
                    B.put((n - 1) * M + (k + 1), row, -m1);
                }
            }
        }
        return B;
    }

    /**
     * The cyclic resolvent blend over the 24 hourly environments.
     *
     * <p>One visit to stage h applies `s*pi*(sI-Q_h)^{-1}` with `s =
     * 1/E[duration]`, which is the LAW OF THE STATE AT A MEMORYLESS EXIT: the
     * resolvent is the time-average over an exponential horizon, so no
     * quadrature grid is needed and the exit vector is exact rather than
     * sampled. Entry vectors are chained around the cycle to an L1 fixed point,
     * and the day average weights each stage by its time fraction.
     *
     * @return the two Table-13 measures, mean external wait W and mean external
     *         queue length Qex, in that order
     */
    private static double[] blendSOQN(int N, double[] lambda, double[] durMin, double[] fracW,
                                      double[] mu1, double[] mu2, int tailFactor) {
        int K = lambda.length;
        int M = N + 1;
        int Mtr = N + tailFactor * N;
        int dim = (Mtr + 1) * M;
        Banded[] A = new Banded[K];    // factored ONCE; the sweeps only back-solve
        for (int h = 0; h < K; h++) {
            A[h] = soqnResolventMatrix(N, lambda[h], mu1, mu2, tailFactor, 1.0 / durMin[h]);
            A[h].factor();
        }
        double[][] piEnter = new double[K][dim];
        for (int h = 0; h < K; h++) {
            piEnter[h][0] = 1.0;       // start every stage empty
        }

        for (int it = 0; it < 200; it++) {
            double[][] prev = new double[K][];
            for (int h = 0; h < K; h++) {
                prev[h] = piEnter[h].clone();
            }
            for (int h = 0; h < K; h++) {
                double s = 1.0 / durMin[h];
                double[] rhs = new double[dim];
                for (int i = 0; i < dim; i++) {
                    rhs[i] = s * piEnter[h][i];
                }
                A[h].solve(rhs);
                piEnter[(h + 1) % K] = rhs;   // the exit law of h is the entry law of h+1
            }
            double l1 = 0;
            for (int h = 0; h < K; h++) {
                double d = 0;
                for (int i = 0; i < dim; i++) {
                    d += Math.abs(piEnter[h][i] - prev[h][i]);
                }
                l1 = Math.max(l1, d);
            }
            if (l1 < 1e-10) {
                break;
            }
        }

        double[] piAvg = new double[dim];
        for (int h = 0; h < K; h++) {
            double s = 1.0 / durMin[h];
            double[] rhs = new double[dim];
            for (int i = 0; i < dim; i++) {
                rhs[i] = s * piEnter[h][i];
            }
            A[h].solve(rhs);
            for (int i = 0; i < dim; i++) {
                piAvg[i] += fracW[h] * rhs[i];
            }
        }
        double tot = 0;
        for (int i = 0; i < dim; i++) {
            // The truncation leaves the blend a hair off a probability vector;
            // the clamp and the renormalization are the reference's own step.
            if (piAvg[i] < 0) {
                piAvg[i] = 0;
            }
            tot += piAvg[i];
        }
        for (int i = 0; i < dim; i++) {
            piAvg[i] /= tot;
        }

        double qlex = 0;
        double ql1 = 0;
        double ql2 = 0;
        double tput = 0;
        for (int n = 0; n <= Mtr; n++) {
            for (int k = 0; k <= N; k++) {
                double p = piAvg[n * M + k];
                if (p == 0) {
                    continue;
                }
                // A token pool of N holds n + k jobs at most; the excess waits
                // OUTSIDE the network, and that is the queue the study reports.
                qlex += Math.max(0, n + k - N) * p;
                ql1 += Math.min(n, N - k) * p;
                ql2 += k * p;
                if (k >= 1) {
                    tput += mu2[k] * p;
                }
            }
        }
        return new double[]{(qlex + ql1 + ql2) / tput, qlex};
    }

    /**
     * `renv_rotterdam_blending.m`: the 24-environment exponential blending
     * accuracy result of the Rotterdam container-terminal study, against its
     * published simulation.
     *
     * <p>The terminal is a SEMI-OPEN network: trucks arrive as a Poisson stream
     * at the hour's rate, a pool of N tokens admits them, and a truck that
     * finds no token waits outside. Inside, two flow-equivalent servers in
     * tandem stand for the upstream half (entry gates, travel, 29 stacks) and
     * the downstream half (travel to exit, exit gates), each calibrated by
     * exact MVA on the closed subnetwork it replaces. That gives a
     * level-dependent QBD whose blocks are the same ones the ENV state-vector
     * analyzer assembles internally, which is why this example builds them
     * directly rather than through SolverENV: the quantities it reports
     * (external wait and external queue) live OUTSIDE the token pool and so
     * outside any stage network's own metric table.
     *
     * <p>The comparison is against Table 13 of the study's discrete-event
     * simulation.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void renv_rotterdam_blending() throws Exception {
        double[] rateHr = {6, 30, 40, 62, 76, 79, 119, 164, 152, 130, 79, 70,
                           57, 57, 113, 130, 162, 202, 148, 118, 92, 62, 36, 8};
        double[] durHr = {0.51, 0.84, 0.76, 0.98, 0.67, 2.33, 1.80, 0.62,
                          0.36, 1.30, 1.02, 0.98, 0.86, 1.56, 1.30, 1.57,
                          0.34, 0.22, 0.47, 1.33, 1.56, 0.93, 0.71, 1.11};
        int K = rateHr.length;
        double[] lambda = new double[K];
        double[] durMin = new double[K];
        double[] fracW = new double[K];
        double totDur = 0;
        for (int h = 0; h < K; h++) {
            lambda[h] = rateHr[h] * 0.5 / 60.0;   // arrivals per MINUTE
            durMin[h] = durHr[h] * 60.0;
            totDur += durMin[h];
        }
        for (int h = 0; h < K; h++) {
            fracW[h] = durMin[h] / totDur;
        }

        // The study's simulated ground truth, N = 24..34 (Table 13).
        double[] simW = {154.0661, 118.3936, 99.962, 88.2379, 80.9272, 75.0887,
                         70.745, 67.4326, 64.9719, 62.7913, 61.2534};
        double[] simQex = {89.8508, 63.3970, 49.3650, 40.8783, 35.3546, 30.6543,
                           27.4111, 24.7218, 22.3622, 20.3893, 18.9252};
        int tailFactor = 15;                      // M_trunc = N*(1 + tailFactor)
        SoqnParams p = new SoqnParams();

        System.out.printf("Rotterdam SOQN exponential blending vs Dhingra simulation (tailFactor=%d)%n",
                tailFactor);
        System.out.printf("%4s | %10s %10s %7s | %10s %10s %7s%n", "N", "W_blend", "W_sim", "err%",
                "Qex_bl", "Qex_sim", "err%");
        int[] nlist = {24, 28, 34};
        for (int i = 0; i < nlist.length; i++) {
            int N = nlist[i];
            double[] mu1 = fesCurve(N, p, true);
            double[] mu2 = fesCurve(N, p, false);
            double[] r = blendSOQN(N, lambda, durMin, fracW, mu1, mu2, tailFactor);
            int j = N - 24;
            System.out.printf("%4d | %10.4f %10.4f %6.2f | %10.4f %10.4f %6.2f%n", N, r[0], simW[j],
                    100 * Math.abs(r[0] - simW[j]) / simW[j], r[1], simQex[j],
                    100 * Math.abs(r[1] - simQex[j]) / simQex[j]);
        }
        pauseForUser();
    }

    /**
     * `renv_lqn_twostages.m`: the LQN stage, two structurally identical copies
     * of which differ only in the database activity's host demand.
     *
     * @param name   the model name
     * @param dbMean mean host demand of the DB activity
     * @return the configured layered network
     */
    private static LayeredNetwork buildLQN(String name, double dbMean) {
        LayeredNetwork model = new LayeredNetwork(name);
        Processor p1 = new Processor(model, "ClientProcessor", 1, SchedStrategy.PS);
        Processor p2 = new Processor(model, "DBProcessor", 1, SchedStrategy.PS);
        Task t1 = new Task(model, "ClientTask", 5, SchedStrategy.REF).on(p1);
        t1.setThinkTime(Exp.fitMean(5.0));
        Task t2 = new Task(model, "DBTask", Integer.MAX_VALUE, SchedStrategy.INF).on(p2);
        Entry e1 = new Entry(model, "ClientEntry").on(t1);
        Entry e2 = new Entry(model, "DBEntry").on(t2);
        Activity a1 = new Activity(model, "ClientActivity", Exp.fitMean(1.0)).on(t1);
        a1.boundTo(e1);
        a1.synchCall(e2, 2.5);
        Activity a2 = new Activity(model, "DBActivity", Exp.fitMean(dbMean)).on(t2);
        a2.boundTo(e2);
        a2.repliesTo(e2);
        return model;
    }

    /**
     * An LQN-in-ENV solver over SolverLN(., FLD).
     *
     * <p>The transient window is set on SolverLN, NOT on the layer factory: the
     * layered fixed point solves each layer in steady state, and SolverLN
     * applies the timespan only to the per-layer transient call SolverENV
     * reads.
     */
    private static SolverENV lqnEnvSolver(Environment env, double horizon, int iterMax, double iterTol) {
        int E = env.getNumberOfStages();
        SolverOptions fldOpt = new SolverOptions(SolverType.FLUID);
        fldOpt.verbose = VerboseLevel.SILENT;
        SolverFactory fldFactory = new SolverFactory() {
            public NetworkSolver at(Network m) {
                return new FLD(m, fldOpt.copy());
            }
        };
        SolverOptions lnOpt = new SolverOptions(SolverType.LN);
        lnOpt.timespan = new double[]{0, horizon};
        lnOpt.verbose = VerboseLevel.SILENT;

        Solver[] solvers = new Solver[E];
        Model[] stageModels = env.getStageModels();
        for (int e = 0; e < E; e++) {
            solvers[e] = new SolverLN((LayeredNetwork) stageModels[e], fldFactory, lnOpt.copy());
        }
        SolverOptions envOpt = new SolverOptions(SolverType.ENV);
        envOpt.iter_max = iterMax;
        envOpt.iter_tol = iterTol;
        envOpt.verbose = VerboseLevel.SILENT;
        return new SolverENV(env, solvers, envOpt);
    }

    /** The finite entries of a matrix, summed; a disabled cell carries NaN. */
    private static double finiteSum(Matrix A) {
        double s = 0;
        for (int i = 0; i < A.getNumRows(); i++) {
            for (int k = 0; k < A.getNumCols(); k++) {
                double v = A.get(i, k);
                if (!Double.isNaN(v) && !Double.isInfinite(v)) {
                    s += v;
                }
            }
        }
        return s;
    }

    /** Env-averaged total throughput of a two-stage UP/DOWN environment. */
    private static double env2StageTput(double a, double b, double horizon) {
        Environment env = new Environment("R", 2);
        env.addStage(0, "UP", "operational", buildLQN("UP", 0.8));
        env.addStage(1, "DOWN", "degraded", buildLQN("DOWN", 3.0));
        env.addTransition(0, 1, new Exp(a));
        env.addTransition(1, 0, new Exp(b));
        env.init();
        SolverENV s = lqnEnvSolver(env, horizon, 10, 0.03);
        s.getEnsembleAvg();
        return finiteSum(s.result.TN);
    }

    /** Env-averaged total throughput of a three-stage UP/MID/DOWN environment. */
    private static double env3StageTput(double horizon) {
        Environment env = new Environment("R3", 3);
        env.addStage(0, "UP", "operational", buildLQN("UP", 0.8));
        env.addStage(1, "MID", "degraded", buildLQN("MID", 1.6));
        env.addStage(2, "DOWN", "failed", buildLQN("DOWN", 3.0));
        env.addTransition(0, 1, new Exp(0.3));
        env.addTransition(1, 2, new Exp(0.3));
        env.addTransition(2, 1, new Exp(0.6));
        env.addTransition(1, 0, new Exp(0.6));
        env.init();
        SolverENV s = lqnEnvSolver(env, horizon, 10, 0.03);
        s.getEnsembleAvg();
        return finiteSum(s.result.TN);
    }

    /**
     * The steady transient aggregate of one LQN stage, in the same
     * block-diagonal layout SolverENV consumes.
     *
     * @return the aggregate queue lengths and throughputs, in that order
     */
    private static Matrix[] stageAggregate(LayeredNetwork model, double horizon) {
        SolverOptions fldOpt = new SolverOptions(SolverType.FLUID);
        fldOpt.verbose = VerboseLevel.SILENT;
        SolverFactory fldFactory = new SolverFactory() {
            public NetworkSolver at(Network m) {
                return new FLD(m, fldOpt.copy());
            }
        };
        SolverOptions lnOpt = new SolverOptions(SolverType.LN);
        lnOpt.timespan = new double[]{0, horizon};
        lnOpt.verbose = VerboseLevel.SILENT;
        LNTranAvgResult tr = new SolverLN(model, fldFactory, lnOpt).getTranAvg();
        int M = tr.QNt.length;
        int K = tr.QNt[0].length;
        Matrix Q = new Matrix(M, K);
        Matrix T = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                Q.set(i, k, tailValue(tr.QNt[i][k]));
                T.set(i, k, tailValue(tr.TNt[i][k]));
            }
        }
        return new Matrix[]{Q, T};
    }

    /** The last sample of a transient series, or zero where the cell is empty. */
    private static double tailValue(Matrix series) {
        if (series == null || series.length() == 0) {
            return 0;
        }
        return series.get(series.length() - 1);
    }

    /**
     * `renv_lqn_twostages.m`: a LayeredNetwork operating in a two-stage random
     * environment.
     *
     * <p>SolverENV runs over an LQN base model through the uniform model/solver
     * interface, with no branching in ENV itself. The environment alternates
     * between an UP stage (fast database) and a DOWN stage (slow database), and
     * the environment-averaged total throughput must lie between the two
     * single-stage LQN solutions: a slower database slows the whole system, so
     * the average is bracketed rather than free.
     *
     * <p>Three further properties are checked because each would fail
     * differently: population conservation (the aggregate is still the closed
     * stage population), monotonicity in P(UP) (the blend weights the stages by
     * their stationary probability), and the same bracket under a THREE-stage
     * environment, which exercises the E &gt; 2 coupling rather than the
     * two-stage special case.
     *
     * @throws Exception if the solver encounters an error
     */
    public static void renv_lqn_twostages() throws Exception {
        LayeredNetwork upModel = buildLQN("LQN_UP", 0.8);      // fast DB activity
        LayeredNetwork downModel = buildLQN("LQN_DOWN", 3.0);  // slow DB activity (degraded)

        Environment env = new Environment("DBReliability", 2);
        env.addStage(0, "UP", "operational", upModel);
        env.addStage(1, "DOWN", "degraded", downModel);
        env.addTransition(0, 1, new Exp(0.2));   // mean UP time = 5
        env.addTransition(1, 0, new Exp(1.0));   // mean DOWN time = 1
        env.init();

        double T = 50;
        SolverENV envSolver = lqnEnvSolver(env, T, 20, 0.02);
        envSolver.getEnsembleAvg();
        Matrix QN = envSolver.result.QN;
        Matrix TN = envSolver.result.TN;
        System.out.printf("ENV over LQN ran: aggregate size = %d x %d%n",
                QN.getNumRows(), QN.getNumCols());
        envSolver.getAvgTable().print();

        Matrix[] up = stageAggregate(upModel, T);
        Matrix[] down = stageAggregate(downModel, T);

        // (1) The solver runs and returns finite, population-conserving metrics.
        for (int i = 0; i < QN.getNumRows(); i++) {
            for (int k = 0; k < QN.getNumCols(); k++) {
                if (Double.isNaN(QN.get(i, k)) || Double.isInfinite(QN.get(i, k))) {
                    throw new RuntimeException("ENV aggregate Q contains non-finite values");
                }
            }
        }
        double sumQ = finiteSum(QN);
        if (Math.abs(sumQ - finiteSum(up[0])) >= 1e-2 || Math.abs(sumQ - finiteSum(down[0])) >= 1e-2) {
            throw new RuntimeException("ENV aggregate does not conserve the closed population of the stages");
        }

        // (2) A physically monotone scalar (total throughput) must lie between
        // the two single-stage solutions. Disabled station/class cells carry
        // NaN and contribute zero, matching the stage references.
        double xUp = finiteSum(up[1]);
        double xDown = finiteSum(down[1]);
        double xEnv = finiteSum(TN);
        double lo = Math.min(xUp, xDown);
        double hi = Math.max(xUp, xDown);
        double tol = 1e-2 * Math.max(1, hi);
        System.out.printf("Total throughput  UP=%.4f  DOWN=%.4f  ENV=%.4f%n", xUp, xDown, xEnv);
        System.out.printf("Aggregate Q (sum) UP=%.4f  DOWN=%.4f  ENV=%.4f%n",
                finiteSum(up[0]), finiteSum(down[0]), sumQ);
        if (xEnv < lo - tol || xEnv > hi + tol) {
            throw new RuntimeException("ENV-averaged throughput is not bracketed by the single-stage solutions");
        }

        // (3) Quantitative coupling: the env-averaged throughput increases with
        // the stationary probability of the fast UP stage, P(UP) = b/(a+b) for
        // switch rates a (UP->DOWN) and b (DOWN->UP).
        double xLow = env2StageTput(1.0, 0.2, T);    // P(UP) = 0.167 (mostly slow DOWN)
        double xHigh = env2StageTput(0.2, 1.0, T);   // P(UP) = 0.833 (mostly fast UP)
        System.out.printf("Monotonicity   P(UP)=0.167 -> %.4f   P(UP)=0.833 -> %.4f%n", xLow, xHigh);
        if (xHigh <= xLow + 1e-3) {
            throw new RuntimeException("ENV throughput is not monotone in P(UP)");
        }
        if (xLow < lo - tol || xHigh > hi + tol) {
            throw new RuntimeException("ENV throughputs escape the single-stage bracket");
        }

        // (4) A three-stage environment stays bracketed by the extreme
        // single-stage solutions, which exercises the E > 2 coupling.
        double x3 = env3StageTput(T);
        System.out.printf("Three-stage ENV throughput = %.4f (bracket [%.4f, %.4f])%n", x3, lo, hi);
        if (x3 < lo - tol || x3 > hi + tol) {
            throw new RuntimeException("Three-stage ENV-averaged throughput is not bracketed");
        }

        System.out.println("PASS: ENV-over-LQN meanfield ran; throughput bracketed, monotone in P(UP), 3-stage bracketed.");
        pauseForUser();
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
        
        for (String name : new String[]{"renv_basic", "renv_genqn", "renv_map_fallback",
                                        "example_mapqn2renv", "renv_lqn_twostages",
                                        "renv_container_terminal", "renv_rotterdam_blending"}) {
            System.out.println("\n=== Running example: " + name + " ===");
            try {
                RandomEnvExamples.class.getMethod(name).invoke(null);
            } catch (Exception e) {
                System.err.println(name + " failed: " + e.getMessage());
                e.printStackTrace();
            }
        }

        scanner.close();
    }
}