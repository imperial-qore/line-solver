/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.des;

import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SolverType;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkElement;
import jline.lang.layered.LayeredNetworkStruct;
import jline.solvers.AvgHandle;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import jline.streaming.Collector;
import jline.streaming.StreamingOptions;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.DistributionResult;
import jline.lang.nodes.StatefulNode;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static jline.api.sn.SnGetArvRFromTputKt.snGetArvRFromTput;
import static jline.io.InputOutputKt.line_debug;
import static jline.solvers.des.analyzers.Solver_des_analyzerKt.solver_des_analyzer;
import static jline.solvers.des.analyzers.Solver_des_analyzer_parallelKt.solver_des_analyzer_parallel;
import static jline.solvers.des.analyzers.Solver_des_ln_analyzerKt.solver_des_ln_analyzer;

/**
 * @brief Discrete Event Simulation (DES) solver using SSJ library.
 *
 * @details SolverDES implements a discrete-event simulation solver that uses the SSJ
 * (Stochastic Simulation in Java) library to analyze queueing networks and
 * stochastic Petri nets.
 *
 * @section des_features Supported Features
 *
 * @subsection des_nodes Node Types
 * - **Source**: Job arrival points for open classes
 * - **Sink**: Job departure points for open classes
 * - **Queue**: Service stations with various scheduling disciplines
 * - **Delay**: Infinite-server stations (think time)
 * - **Fork/Join**: Parallel processing with synchronization (with quorum support)
 * - **Router**: Routing decisions with multiple strategies
 * - **ClassSwitch**: Dynamic job class switching
 * - **Logger/LogTunnel**: Job passage logging for trace collection
 * - **Place/Transition**: Stochastic Petri net modeling
 *
 * @subsection des_sched Scheduling Strategies
 * - **FCFS**: First-Come-First-Served
 * - **LCFS/LCFSPR/LCFSPI**: Last-Come-First-Served (non-preemptive, preemptive-resume, preemptive-independent)
 * - **PS/DPS/GPS**: Processor Sharing variants (standard, discriminatory, generalized)
 * - **HOL**: Head-of-Line priority scheduling
 * - **SIRO**: Service In Random Order
 * - **SJF/LJF**: Shortest/Longest Job First
 * - **SEPT/LEPT**: Shortest/Longest Expected Processing Time
 * - **INF**: Infinite server (for Delay nodes)
 * - Priority variants: FCFSPRIO, LCFSPRIO, LCFSPRPRIO, LCFSPIPRIO, PSPRIO, DPSPRIO, GPSPRIO
 *
 * @subsection des_routing Routing Strategies
 * - **PROB**: Probabilistic routing based on routing matrix
 * - **RAND**: Uniform random selection among destinations
 * - **RROBIN**: Round-robin cycling through destinations
 * - **WRROBIN**: Weighted round-robin
 * - **KCHOICES**: Power-of-K-choices (select shortest queue among K random samples)
 *
 * @subsection des_dist Service Time Distributions
 * - **Exp**: Exponential distribution
 * - **Erlang**: Erlang distribution (sum of exponentials)
 * - **HyperExp**: Hyperexponential (mixture of exponentials)
 * - **PH/APH**: Phase-type distributions (general and acyclic)
 * - **Coxian**: Coxian distribution
 * - **Immediate**: Zero service time (instantaneous)
 * - **Disabled**: No service (jobs bypass station)
 * - **Replayer**: Trace-driven service times from file
 *
 * @subsection des_classes Job Class Types
 * - **OpenClass**: Jobs arrive from external source, depart to sink
 * - **ClosedClass**: Fixed population circulating in the network
 * - **SelfLoopingClass**: Jobs remain at reference station
 *
 * @subsection des_advanced Advanced Features
 * - **Load Dependence**: Service rates varying with queue population
 * - **Finite Capacity**: Buffer limits with blocking
 * - **Multiserver**: Multiple parallel servers at a station
 * - **Class Priorities**: Priority-based scheduling across classes
 *
 * @section des_petri Stochastic Petri Net Support
 * DES supports stochastic Petri nets with:
 * - **Place nodes**: Token storage with capacity constraints
 * - **Transition nodes**: Token firing with enabling/inhibiting conditions
 * - **Timed transitions**: Exponential or phase-type firing delays
 * - **Immediate transitions**: Zero-delay firing with priority and weight
 *
 * @section des_analysis Analysis Modes
 * - **Steady-state**: Long-run average performance metrics (default)
 * - **Transient**: Time-varying metrics over specified time horizon
 *
 * @section des_metrics Output Metrics
 * - Queue length (QN): Average number of jobs at each station
 * - Utilization (UN): Fraction of time servers are busy
 * - Response time (RN): Average time from arrival to departure
 * - Throughput (TN): Average job completion rate
 * - Arrival rate (AN): Average job arrival rate
 *
 * @section des_example Example Usage
 * @code{.java}
 * Network model = new Network("M/M/1");
 * Source source = new Source(model, "Source");
 * Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
 * Sink sink = new Sink(model, "Sink");
 *
 * OpenClass jobClass = new OpenClass(model, "Jobs");
 * source.setArrival(jobClass, Exp.fitMean(1.0));
 * queue.setService(jobClass, Exp.fitMean(0.5));
 *
 * model.link(model.serialRouting(source, queue, sink));
 *
 * SolverDES solver = new SolverDES(model, "seed", 12345, "samples", 100000);
 * SolverResult result = solver.getAvg();
 * @endcode
 *
 * @see DESResult Result container for DES metrics
 * @see DESOptions Configuration options for DES solver
 * @see NetworkSolver Base class for network solvers
 * @since 1.0
 * @author QORE Lab, Imperial College London
 */
public class SolverDES extends NetworkSolver {

    /** LayeredNetwork model (non-null if solving LQN) */
    private LayeredNetwork lnModel;

    /** LayeredNetwork structure (populated if solving LQN) */
    private LayeredNetworkStruct lsn;

    /** Streaming collector for pushing metrics via OTLP */
    private Collector stream;

    /** Thread pool for parallel replication execution */
    private ExecutorService threadPool;

    /** Number of parallel threads for replication execution */
    private int numThreads;

    /**
     * Constructs a SolverDES with the given model using default options.
     *
     * @param model the queueing network model to solve
     */
    public SolverDES(Network model) {
        this(model, new DESOptions());
        this.result = new DESResult();
    }

    /**
     * Constructs a SolverDES with the given LayeredNetwork model using default options.
     *
     * @param model the layered queueing network model to solve
     */
    public SolverDES(LayeredNetwork model) {
        this(model, new DESOptions());
    }

    /**
     * Constructs a SolverDES with the given LayeredNetwork model and variable arguments.
     *
     * @param model the layered queueing network model to solve
     * @param args variable arguments for solver options
     */
    public SolverDES(LayeredNetwork model, Object... args) {
        super(null, "SolverDES");
        this.lnModel = model;
        this.lsn = model.getStruct(true);
        this.setOptions(Solver.parseOptions(new DESOptions(), args));
        this.result = new LNDESResult();
    }

    /**
     * Constructs a SolverDES with the given LayeredNetwork model and options.
     *
     * @param model the layered queueing network model to solve
     * @param options solver configuration options
     */
    public SolverDES(LayeredNetwork model, SolverOptions options) {
        super(null, "SolverDES", options);
        this.lnModel = model;
        this.lsn = model.getStruct(true);
        this.result = new LNDESResult();
    }

    /**
     * Constructs a SolverDES with the given model and variable arguments.
     *
     * @param model the queueing network model to solve
     * @param args variable arguments for solver options
     */
    public SolverDES(Network model, Object... args) {
        super(model, "SolverDES");
        this.setOptions(Solver.parseOptions(new DESOptions(), args));
        this.result = new DESResult();
    }

    /**
     * Constructs a SolverDES with the given model and method.
     *
     * @param model the queueing network model to solve
     * @param method the solution method to use
     */
    public SolverDES(Network model, String method) {
        super(model, "SolverDES", new DESOptions().method(method));
        this.result = new DESResult();
    }

    /**
     * Constructs a SolverDES with the given model and options.
     *
     * @param model the queueing network model to solve
     * @param options solver configuration options
     */
    public SolverDES(Network model, SolverOptions options) {
        super(model, "SolverDES", options);
        this.result = new DESResult();
    }

    // =====================================================
    // Thread Pool Management for Parallel Replications
    // =====================================================

    /**
     * Returns the thread pool for parallel replication execution.
     *
     * @return the thread pool
     */
    public ExecutorService getThreadPool() {
        return threadPool;
    }

    /**
     * Returns the number of threads for parallel execution.
     *
     * @return the number of threads
     */
    public int getNumThreads() {
        return numThreads;
    }

    /**
     * Sets the parallelism level for replication execution.
     *
     * @param numThreads number of parallel threads
     */
    public void setNumThreads(int numThreads) {
        this.numThreads = numThreads;
    }

    /**
     * Returns the feature set supported by the DES solver.
     * Supports multiclass Jackson queueing networks with FCFS queues.
     *
     * @return the feature set supported by the DES solver
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Sink", "Source",
                "Queue", "Delay",
                "Fork", "Join", "Forker", "Joiner",  // Fork-Join node support
                "Place", "Transition",  // Stochastic Petri Net support
                "Linkage", "Enabling", "Timing", "Firing", "Storage",  // Petri net section support
                "Logger", "LogTunnel",  // Logger node support
                "Cache", "CacheClassSwitcher",  // Cache node support with hit/miss class switching
                "ReplacementStrategy_LRU", "ReplacementStrategy_FIFO", "ReplacementStrategy_RR", "ReplacementStrategy_SFIFO",  // Cache replacement strategies
                "Buffer",  // Finite buffer capacity support
                "Region",  // Finite capacity region support
                "Exp", "Erlang", "HyperExp", "PH", "APH", "Coxian", "Cox2", "MAP", "MMAP", "BMAP", "MMPP2", "ME", "RAP", "Immediate", "Disabled", "Replayer",
                "Det", "Uniform", "Gamma", "Pareto", "Weibull", "Lognormal",  // Additional continuous distributions
                "Server", "JobSink", "RandomSource",
                "SchedStrategy_FCFS", "SchedStrategy_INF",
                "SchedStrategy_HOL",  // Priority scheduling (FCFS with priorities)
                "SchedStrategy_PS",   // Processor Sharing
                "SchedStrategy_DPS",  // Discriminatory Processor Sharing
                "SchedStrategy_GPS",  // Generalized Processor Sharing
                "SchedStrategy_PSPRIO",   // PS with priorities
                "SchedStrategy_DPSPRIO",  // DPS with priorities
                "SchedStrategy_GPSPRIO",  // GPS with priorities
                "SchedStrategy_LCFS",     // Last Come First Served (non-preemptive)
                "SchedStrategy_LCFSPR",   // LCFS Preemptive Resume
                "SchedStrategy_LCFSPI",   // LCFS Preemptive Independent
                "SchedStrategy_LCFSPRIO",     // LCFS with priorities (non-preemptive)
                "SchedStrategy_LCFSPRPRIO",   // LCFSPR with priorities
                "SchedStrategy_LCFSPIPRIO",   // LCFSPI with priorities
                "SchedStrategy_SIRO",
                "SchedStrategy_SJF",  // Shortest Job First
                "SchedStrategy_LJF",  // Longest Job First
                "SchedStrategy_LEPT",
                "SchedStrategy_SEPT",
                "SchedStrategy_SRPT",  // Shortest Remaining Processing Time (preemptive)
                "SchedStrategy_SRPTPRIO",  // SRPT with priorities
                "SchedStrategy_PSJF",  // Preemptive Shortest Job First
                "SchedStrategy_FB",  // Feedback / Least Attained Service
                "SchedStrategy_LRPT",  // Longest Remaining Processing Time
                "SchedStrategy_EXT",
                "SchedStrategy_POLLING",  // Polling scheduling (GATED, EXHAUSTIVE, KLIMITED)
                "Router",  // Router node support
                "ClassSwitch", "StatelessClassSwitcher",  // Class switching node support
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "RoutingStrategy_RROBIN", "RoutingStrategy_WRROBIN",
                "RoutingStrategy_KCHOICES",  // Power of K Choices routing
                "OpenClass",
                "ClosedClass",
                "SelfLoopingClass",
                "LoadDependence"  // Load-dependent service rates
        });
        return featSupported;
    }

    /**
     * Returns the network structure for this solver.
     *
     * @return the network structure
     */
    public NetworkStruct getStruct() {
        return this.model.getStruct(true);
    }

    /**
     * Performs a transient analysis of the model using discrete event simulation.
     * The simulation is run for the time horizon specified in options.timespan.
     * If options.timespan is not specified, it defaults to [0, 30/min_rate].
     * The number of replications is determined by options.samples.
     */
    public void getTranAvg() {
        // Ensure sn is initialized
        this.sn = this.model.getStruct(true);

        // Default timespan if not set or infinite
        if (this.options.timespan == null || Double.isInfinite(this.options.timespan[1])) {
            double minRate = Double.MAX_VALUE;
            // Check service rates and arrival rates
            for (int i = 0; i < this.sn.rates.length(); i++) {
                double r = this.sn.rates.get(i);
                if (r > 0 && r < minRate) {
                    minRate = r;
                }
            }
            if (minRate == Double.MAX_VALUE) {
                minRate = 1.0;
            }
            this.options.timespan = new double[]{0.0, 30.0 / minRate};
        }

        try {
            runAnalyzer();
        } catch (Exception e) {
            throw new RuntimeException("Transient analysis failed.", e);
        }
    }

    /**
     * Computes and returns average station metrics at steady-state.
     * For LayeredNetwork models, this delegates to the LN-specific analyzer.
     *
     * @return solver result containing average metrics
     */
    @Override
    public SolverResult getAvg() {
        // For LayeredNetwork models, use dedicated LN analyzer
        if (this.lnModel != null) {
            try {
                runAnalyzer();
            } catch (Exception e) {
                throw new RuntimeException("DES LN analysis failed: " + e.getMessage(), e);
            }
            return this.result;
        }
        // For Network models, use parent implementation
        return super.getAvg();
    }

    /**
     * Returns the list of valid solution methods for DES.
     *
     * @return list of valid method names
     */
    public List<String> listValidMethods() {
        return listValidMethods(null);
    }

    /**
     * Returns the list of valid solution methods for DES.
     *
     * @param model the network model (unused, for interface compatibility)
     * @return list of valid method names
     */
    public List<String> listValidMethods(Network model) {
        return Arrays.asList("default", "parallel");
    }

    @Override
    public void runAnalyzer() throws IllegalAccessException, ParserConfigurationException, IOException {
        long T0 = System.nanoTime();
        if (this.options == null) {
            this.options = new SolverOptions(SolverType.DES);
        }

        // Check if we're solving a LayeredNetwork model
        if (this.lnModel != null) {
            runLNAnalyzer();
            return;
        }

        if (this.enableChecks && !supports(this.model)) {
            throw new RuntimeException("This model is not supported by the DES solver.");
        }
        this.resetRandomGeneratorSeed(options.seed);
        String method = options.method;

        // Check for parallel replication mode (only if options is DESOptions)
        DESOptions desOptions = null;
        boolean useParallel = false;
        if (this.options instanceof DESOptions) {
            desOptions = (DESOptions) this.options;
            useParallel = desOptions.replications > 1;
        }

        if (useParallel) {
            // Initialize thread pool for parallel replications
            this.numThreads = Math.min(desOptions.numThreads, desOptions.replications);
            this.threadPool = Executors.newFixedThreadPool(this.numThreads);
            line_debug(options.verbose, String.format("DES solver starting (parallel): replications=%d, threads=%d, samples=%d, seed=%d",
                    desOptions.replications, this.numThreads, options.samples, options.seed));
        } else {
            line_debug(options.verbose, String.format("DES solver starting: method=%s, samples=%d, seed=%d",
                    method, options.samples, options.seed));
        }
        line_debug(options.verbose, "Running DES simulation, calling solver_des_analyzer");

        DESResult result;
        try {
            if (useParallel) {
                result = solver_des_analyzer_parallel(this.sn, desOptions, this);
            } else {
                result = solver_des_analyzer(this.sn, this.options, this);
            }
        } catch (RuntimeException e) {
            e.printStackTrace();
            throw new RuntimeException("DES simulation failed.", e);
        } finally {
            // Shutdown thread pool if it was used
            if (this.threadPool != null) {
                this.threadPool.shutdown();
                try {
                    this.threadPool.awaitTermination(10, TimeUnit.SECONDS);
                } catch (InterruptedException ex) {
                    this.threadPool.shutdownNow();
                }
                this.threadPool = null;
            }
        }

        // Extract result matrices
        Matrix QN = result.QN;
        Matrix UN = result.UN;
        Matrix RN = result.RN;
        Matrix TN = result.TN;
        Matrix TardN = result.TardN;
        Matrix SysTardN = result.SysTardN;
        Matrix CN = result.CN;
        Matrix XN = result.XN;

        NetworkStruct sn = result.sn;

        double runtime = result.runtime;
        int M = sn.nstations;
        int R = sn.nclasses;
        // Use arrival rate from DES result (includes dropped jobs), fallback to computed if null
        Matrix AN = result.AN != null ? result.AN : snGetArvRFromTput(sn, TN, getAvgTputHandles());
        Matrix WN = new Matrix(0, 0);

        this.result.method = result.method;
        this.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN, runtime, this.result.method, options.samples);

        // Set tardiness results
        this.result.TardN = TardN;
        this.result.SysTardN = SysTardN;

        // Transfer confidence interval and convergence data if available
        if (result instanceof DESResult) {
            DESResult desResult = (DESResult) result;
            ((DESResult) this.result).QNCI = desResult.QNCI;
            ((DESResult) this.result).UNCI = desResult.UNCI;
            ((DESResult) this.result).RNCI = desResult.RNCI;
            ((DESResult) this.result).TNCI = desResult.TNCI;
            ((DESResult) this.result).ANCI = desResult.ANCI;
            ((DESResult) this.result).WNCI = desResult.WNCI;
            // Transfer convergence data
            ((DESResult) this.result).converged = desResult.converged;
            ((DESResult) this.result).stoppingReason = desResult.stoppingReason;
            ((DESResult) this.result).convergenceBatches = desResult.convergenceBatches;
            ((DESResult) this.result).QNRelPrec = desResult.QNRelPrec;
            ((DESResult) this.result).UNRelPrec = desResult.UNRelPrec;
            ((DESResult) this.result).RNRelPrec = desResult.RNRelPrec;
            ((DESResult) this.result).TNRelPrec = desResult.TNRelPrec;
            // Transfer impatience statistics
            ((DESResult) this.result).renegedCustomers = desResult.renegedCustomers;
            ((DESResult) this.result).avgRenegingWaitTime = desResult.avgRenegingWaitTime;
            ((DESResult) this.result).renegingRate = desResult.renegingRate;
            ((DESResult) this.result).balkedCustomers = desResult.balkedCustomers;
            ((DESResult) this.result).balkingProbability = desResult.balkingProbability;
            ((DESResult) this.result).retriedCustomers = desResult.retriedCustomers;
            ((DESResult) this.result).retrialDropped = desResult.retrialDropped;
            ((DESResult) this.result).avgOrbitSize = desResult.avgOrbitSize;
        }

        // Transfer FCR (Finite Capacity Region) metrics if available
        this.result.QNfcr = result.QNfcr;
        this.result.UNfcr = result.UNfcr;
        this.result.RNfcr = result.RNfcr;
        this.result.TNfcr = result.TNfcr;
        this.result.ANfcr = result.ANfcr;
        this.result.WNfcr = result.WNfcr;

        // Handle transient results if available
        if (result.QNt != null && result.QNt.length > 0) {
            // Create empty response time transient matrices (not computed by DES)
            Matrix[][] RNt = new Matrix[M][R];
            Matrix[][] CNt = new Matrix[1][R];
            Matrix[][] XNt = new Matrix[1][R];
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    RNt[i][r] = new Matrix(0, 0);
                }
            }
            for (int r = 0; r < R; r++) {
                CNt[0][r] = new Matrix(0, 0);
                XNt[0][r] = new Matrix(0, 0);
            }
            this.setTranAvgResults(result.QNt, result.UNt, RNt, result.TNt, CNt, XNt, runtime);
            if (result.t != null) {
                this.result.t = result.t;
            }
        }
    }

    /**
     * Runs the DES analyzer for LayeredNetwork models.
     */
    private void runLNAnalyzer() {
        this.resetRandomGeneratorSeed(options.seed);
        String method = options.method;
        line_debug(options.verbose, String.format("DES LN solver starting: method=%s, samples=%d, seed=%d",
                method, options.samples, options.seed));
        line_debug(options.verbose, "Running DES LN simulation, calling solver_des_ln_analyzer");

        LNDESResult lnResult;
        try {
            lnResult = solver_des_ln_analyzer(this.lsn, this.options, this);
        } catch (RuntimeException e) {
            e.printStackTrace();
            throw new RuntimeException("DES LN simulation failed.", e);
        }

        // Store results
        this.result = lnResult;
        this.result.method = lnResult.method;
        this.result.runtime = lnResult.runtime;
    }

    /**
     * Run the DES analyzer with current options including init_sol.
     * This method is called by SolverENV for transient analysis with initial conditions.
     *
     * @return SolverResult containing transient metrics (QNt, UNt, TNt, t)
     */
    public SolverResult runMethodSpecificAnalyzer() {
        // Re-initialize sn with current model state (may have been modified by SolverENV)
        this.sn = this.model.getStruct(false);

        // Ensure transient mode by setting finite timespan if not already set
        if (this.options.timespan == null || Double.isInfinite(this.options.timespan[1])) {
            // Calculate default timespan based on minimum rate
            double minRate = Double.MAX_VALUE;
            for (int i = 0; i < this.sn.rates.length(); i++) {
                double r = this.sn.rates.get(i);
                if (r > 0 && r < minRate) {
                    minRate = r;
                }
            }
            if (minRate == Double.MAX_VALUE) {
                minRate = 1.0;
            }
            this.options.timespan = new double[]{0.0, 30.0 / minRate};
        }

        // Call runAnalyzer which delegates to solver_des_analyzer
        try {
            runAnalyzer();
        } catch (Exception e) {
            throw new RuntimeException("DES transient analysis failed", e);
        }

        return this.result;
    }

    @Override
    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverDES.getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    /**
     * Checks if the solver supports the given LayeredNetwork model.
     *
     * @param model the LayeredNetwork model to check
     * @return true if the model is supported
     */
    public boolean supports(LayeredNetwork model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverDES.getLNFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    /**
     * Returns the feature set supported by the DES solver for LayeredNetwork models.
     *
     * @return the feature set supported for LQN models
     */
    public static FeatureSet getLNFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Host", "Processor",
                "Task", "Entry", "Activity",
                "SyncCall", "AsyncCall",
                "ActivityPrecedence_PRE_SEQ", "ActivityPrecedence_POST_SEQ",
                "ActivityPrecedence_PRE_AND", "ActivityPrecedence_POST_AND",
                "ActivityPrecedence_PRE_OR", "ActivityPrecedence_POST_OR",
                "SchedStrategy_REF", "SchedStrategy_FCFS", "SchedStrategy_PS", "SchedStrategy_INF",
                "Exp", "Erlang", "HyperExp", "PH", "APH", "Coxian", "Det", "Uniform", "Gamma"
        });
        return featSupported;
    }

    /**
     * Returns the LayeredNetwork model (if solving an LQN).
     *
     * @return the LayeredNetwork model, or null if solving a regular Network
     */
    public LayeredNetwork getLNModel() {
        return this.lnModel;
    }

    /**
     * Returns the LayeredNetworkStruct (if solving an LQN).
     *
     * @return the LayeredNetworkStruct, or null if solving a regular Network
     */
    public LayeredNetworkStruct getLNStruct() {
        return this.lsn;
    }

    /**
     * Returns average metrics for a LayeredNetwork model as a table.
     * This method should be called after getAvg() for LayeredNetwork models.
     *
     * @return LayeredNetworkAvgTable containing average metrics per LQN element
     * @throws RuntimeException if not solving a LayeredNetwork model or results not available
     */
    public LayeredNetworkAvgTable getLNAvgTable() {
        if (this.lnModel == null) {
            throw new RuntimeException("getLNAvgTable requires a LayeredNetwork model");
        }
        if (!(this.result instanceof LNDESResult)) {
            throw new RuntimeException("getLNAvgTable requires LNDESResult - call getAvg() first");
        }

        LNDESResult lnResult = (LNDESResult) this.result;

        List<Double> Qval = new ArrayList<Double>();
        List<Double> Uval = new ArrayList<Double>();
        List<Double> Rval = new ArrayList<Double>();
        List<Double> Wval = new ArrayList<Double>();
        List<Double> Aval = new ArrayList<Double>();
        List<Double> Tval = new ArrayList<Double>();

        List<String> nodeNames = new ArrayList<String>();
        List<String> nodeTypes = new ArrayList<String>();

        // Populate from lnResult matrices - iterate through all LQN elements
        for (int i = 1; i <= this.lsn.nidx; i++) {
            Qval.add(lnResult.QLN != null ? lnResult.QLN.get(0, i) : 0.0);
            Uval.add(lnResult.ULN != null ? lnResult.ULN.get(0, i) : 0.0);
            Rval.add(lnResult.RLN != null ? lnResult.RLN.get(0, i) : 0.0);
            Wval.add(lnResult.WLN != null ? lnResult.WLN.get(0, i) : 0.0);
            Aval.add(lnResult.ALN != null ? lnResult.ALN.get(0, i) : 0.0);
            Tval.add(lnResult.TLN != null ? lnResult.TLN.get(0, i) : 0.0);

            nodeNames.add(this.lsn.names.get(i));
            nodeTypes.add(getNodeTypeName((int) this.lsn.type.get(0, i)));
        }

        LayeredNetworkAvgTable table = new LayeredNetworkAvgTable(Qval, Uval, Rval, Wval, Aval, Tval);
        table.setNodeNames(nodeNames);
        table.setNodeTypes(nodeTypes);
        table.setOptions(this.options);

        return table;
    }

    /**
     * Returns the type name for a LayeredNetworkElement type constant.
     *
     * @param type the type constant
     * @return the type name string
     */
    private String getNodeTypeName(int type) {
        switch (type) {
            case LayeredNetworkElement.HOST:
                return "Host";
            case LayeredNetworkElement.TASK:
                return "Task";
            case LayeredNetworkElement.ENTRY:
                return "Entry";
            case LayeredNetworkElement.ACTIVITY:
                return "Activity";
            case LayeredNetworkElement.CALL:
                return "Call";
            default:
                return "Unknown";
        }
    }

    /**
     * Returns the default solver options for the DES solver.
     *
     * @return Default solver options with SolverType.DES
     */
    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.DES);
    }

    // =====================================================
    // Sample Methods - State Trajectory Generation
    // =====================================================

    /**
     * Generates a sample path (state trajectory) for a specific node using DES transient simulation.
     *
     * @param node The stateful node to sample
     * @param numSamples Number of time points to sample
     * @return SampleResult containing the state trajectory for the node
     */
    public jline.io.Ret.SampleResult sample(jline.lang.nodes.StatefulNode node, int numSamples) {
        // Configure for transient analysis using options.samples as time horizon
        DESOptions desOptions = (DESOptions) this.options;
        double[] originalTimespan = desOptions.timespan;
        desOptions.timespan = new double[]{0.0, (double) desOptions.samples};

        try {
            // Run transient simulation
            runAnalyzer();

            // Get results
            DESResult desResult = (DESResult) this.result;
            if (desResult == null || desResult.t == null) {
                return new jline.io.Ret.SampleResult();
            }

            NetworkStruct sn = this.sn;
            int nodeIdx = node.getNodeIndex();
            int isf = (int) sn.nodeToStateful.get(0, nodeIdx);

            // Build state matrix from transient queue lengths
            int numTimePoints = desResult.t.getNumRows();
            int numClasses = sn.nclasses;
            Matrix t = desResult.t.copy();
            Matrix state = new Matrix(numTimePoints, numClasses);

            if (desResult.QNt != null && isf < desResult.QNt.length) {
                for (int k = 0; k < numClasses && k < desResult.QNt[isf].length; k++) {
                    Matrix classData = desResult.QNt[isf][k];
                    if (classData != null) {
                        for (int ti = 0; ti < numTimePoints && ti < classData.getNumRows(); ti++) {
                            state.set(ti, k, classData.get(ti, 0));
                        }
                    }
                }
            }

            return new jline.io.Ret.SampleResult("des", t, state, new Matrix(0, 0), false, nodeIdx, numTimePoints);

        } catch (Exception e) {
            line_debug(options.verbose, "DES sample failed: " + e.getMessage());
            return new jline.io.Ret.SampleResult();
        } finally {
            // Restore original timespan
            desOptions.timespan = originalTimespan;
        }
    }

    /**
     * Generates an aggregated sample path for a specific node.
     *
     * For DES, the sample() method already returns per-class queue lengths (nir),
     * so this method simply marks the result as aggregated without additional
     * marginal computation. This differs from SSA which requires toMarginal
     * processing on raw state vectors.
     *
     * @param node The stateful node to sample
     * @param numSamples Number of time points to sample
     * @return SampleResult containing the aggregated state trajectory
     */
    public jline.io.Ret.SampleResult sampleAggr(jline.lang.nodes.StatefulNode node, int numSamples) {
        jline.io.Ret.SampleResult result = sample(node, numSamples);

        if (result != null && result.state != null && result.state instanceof Matrix) {
            // DES sample() already returns aggregated per-class queue lengths from QNt,
            // so we just need to mark the result as aggregated
            Matrix stateMatrix = (Matrix) result.state;
            return new jline.io.Ret.SampleResult("des", result.t, stateMatrix, result.event, true, result.nodeIndex, result.numSamples);
        }
        return result;
    }

    /**
     * Generates a system-wide sample path (state trajectory) for all stateful nodes.
     *
     * @param numSamples Number of time points to sample
     * @return SampleResult containing the state trajectory for all nodes
     */
    public jline.io.Ret.SampleResult sampleSys(int numSamples) {
        // Configure for transient analysis using options.samples as time horizon
        DESOptions desOptions = (DESOptions) this.options;
        double[] originalTimespan = desOptions.timespan;
        desOptions.timespan = new double[]{0.0, (double) desOptions.samples};

        try {
            // Run transient simulation
            runAnalyzer();

            DESResult desResult = (DESResult) this.result;
            if (desResult == null || desResult.t == null) {
                return new jline.io.Ret.SampleResult();
            }

            NetworkStruct sn = this.sn;
            int numTimePoints = desResult.t.getNumRows();
            int numClasses = sn.nclasses;
            Matrix t = desResult.t.copy();

            // Build state list for each stateful node
            List<Matrix> stateList = new ArrayList<>();
            for (int isf = 0; isf < sn.nstateful; isf++) {
                Matrix nodeState = new Matrix(numTimePoints, numClasses);

                if (desResult.QNt != null && isf < desResult.QNt.length) {
                    for (int k = 0; k < numClasses && k < desResult.QNt[isf].length; k++) {
                        Matrix classData = desResult.QNt[isf][k];
                        if (classData != null) {
                            for (int ti = 0; ti < numTimePoints && ti < classData.getNumRows(); ti++) {
                                nodeState.set(ti, k, classData.get(ti, 0));
                            }
                        }
                    }
                }
                stateList.add(nodeState);
            }

            return new jline.io.Ret.SampleResult("des", t, stateList, new Matrix(0, 0), false, numTimePoints);

        } catch (Exception e) {
            line_debug(options.verbose, "DES sampleSys failed: " + e.getMessage());
            return new jline.io.Ret.SampleResult();
        } finally {
            desOptions.timespan = originalTimespan;
        }
    }

    /**
     * Generates an aggregated system-wide sample path for all stateful nodes.
     *
     * @param numSamples Number of time points to sample
     * @return SampleResult containing the aggregated state trajectory for all nodes
     */
    @SuppressWarnings("unchecked")
    public jline.io.Ret.SampleResult sampleSysAggr(int numSamples) {
        jline.io.Ret.SampleResult result = sampleSys(numSamples);

        if (result != null && result.state != null && result.state instanceof List) {
            List<Matrix> stateList = (List<Matrix>) result.state;
            NetworkStruct sn = this.sn;

            List<Matrix> aggregatedStateList = new ArrayList<>();
            for (int isf = 0; isf < stateList.size() && isf < sn.nstateful; isf++) {
                Matrix nodeState = stateList.get(isf);
                int nodeIndex = ((Double) sn.statefulToNode.get(isf)).intValue();

                // Apply marginal aggregation
                jline.lang.state.State.StateMarginalStatistics marginal =
                    jline.lang.state.ToMarginal.toMarginal(sn, nodeIndex, nodeState, null, null, null, null, null);

                aggregatedStateList.add(marginal.nir);
            }

            return new jline.io.Ret.SampleResult("des", result.t, aggregatedStateList, result.event, true, result.numSamples);
        }
        return result;
    }

    /**
     * Generates a sample path with streaming metrics pushed via OTLP.
     *
     * @param node The stateful node to sample
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleResult containing the state trajectory for the node
     */
    public jline.io.Ret.SampleResult stream(jline.lang.nodes.StatefulNode node, StreamingOptions streamingOptions) {
        return stream(node, this.options.samples, streamingOptions);
    }

    /**
     * Generates a sample path with streaming metrics pushed via OTLP.
     *
     * @param node The stateful node to sample
     * @param numSamples Number of time points to sample
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleResult containing the state trajectory for the node
     */
    public jline.io.Ret.SampleResult stream(jline.lang.nodes.StatefulNode node, int numSamples, StreamingOptions streamingOptions) {
        try {
            // Initialize streaming collector
            this.sn = this.model.getStruct(true);
            this.stream = new Collector(streamingOptions, this.sn);

            // Run sample with streaming active
            jline.io.Ret.SampleResult result = sample(node, numSamples);

            // Flush any remaining metrics
            if (result != null && result.t != null && result.t.length() > 0) {
                double finalTime = result.t.get(result.t.length() - 1);
                this.stream.flush(finalTime);
            }

            return result;
        } finally {
            // Cleanup streaming collector
            if (this.stream != null) {
                this.stream.shutdown();
                this.stream = null;
            }
        }
    }

    /**
     * Generates an aggregated sample path with streaming metrics pushed via OTLP.
     *
     * @param node The stateful node to sample
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleResult containing the aggregated state trajectory
     */
    public jline.io.Ret.SampleResult streamAggr(jline.lang.nodes.StatefulNode node, StreamingOptions streamingOptions) {
        return streamAggr(node, this.options.samples, streamingOptions);
    }

    /**
     * Generates an aggregated sample path with streaming metrics pushed via OTLP.
     *
     * @param node The stateful node to sample
     * @param numSamples Number of time points to sample
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleResult containing the aggregated state trajectory
     */
    public jline.io.Ret.SampleResult streamAggr(jline.lang.nodes.StatefulNode node, int numSamples, StreamingOptions streamingOptions) {
        try {
            // Initialize streaming collector
            this.sn = this.model.getStruct(true);
            this.stream = new Collector(streamingOptions, this.sn);

            // Run sampleAggr with streaming active
            jline.io.Ret.SampleResult result = sampleAggr(node, numSamples);

            // Flush any remaining metrics
            if (result != null && result.t != null && result.t.length() > 0) {
                double finalTime = result.t.get(result.t.length() - 1);
                this.stream.flush(finalTime);
            }

            return result;
        } finally {
            // Cleanup streaming collector
            if (this.stream != null) {
                this.stream.shutdown();
                this.stream = null;
            }
        }
    }

    /**
     * Returns the current streaming collector, or null if not streaming.
     *
     * @return the streaming collector
     */
    public Collector getStream() {
        return this.stream;
    }

    // =====================================================
    // Probability Methods - State Probability Estimation
    // =====================================================

    /**
     * Estimates the steady-state probability of a specific state at a node via DES simulation.
     *
     * This method runs a DES simulation and estimates state probabilities by computing
     * the fraction of total simulation time spent in each state. The simulation
     * collects state trajectories and weights each state by its dwell time.
     *
     * @param node The stateful node to analyze
     * @param state The state vector to compute probability for. If null, uses current network state.
     * @return Estimated probability of the specified state (0 if state not observed)
     */
    public ProbabilityResult getProb(jline.lang.nodes.StatefulNode node, Matrix state) {
        // Get sample path via transient simulation
        jline.io.Ret.SampleResult sampleResult = sample(node, this.options.samples);

        if (sampleResult == null || sampleResult.t == null || !(sampleResult.state instanceof Matrix)) {
            return new ProbabilityResult(0.0);
        }

        Matrix t = sampleResult.t;
        Matrix stateMatrix = (Matrix) sampleResult.state;

        int numTimePoints = t.getNumRows();
        if (numTimePoints < 2) {
            return new ProbabilityResult(0.0);
        }

        // Get target state
        Matrix targetState = state;
        if (targetState == null) {
            NetworkStruct sn = this.sn;
            int nodeIdx = node.getNodeIndex();
            int isf = (int) sn.nodeToStateful.get(0, nodeIdx);
            targetState = sn.state.get(isf);
        }

        // Compute time-weighted probability: sum time spent in target state / total time
        double totalTime = t.get(numTimePoints - 1, 0) - t.get(0, 0);
        if (totalTime <= 0) {
            return new ProbabilityResult(0.0);
        }

        double timeInState = 0.0;
        int numClasses = stateMatrix.getNumCols();

        for (int ti = 0; ti < numTimePoints - 1; ti++) {
            double dt = t.get(ti + 1, 0) - t.get(ti, 0);

            // Check if current state matches target state
            boolean matches = true;
            for (int k = 0; k < numClasses && k < targetState.length(); k++) {
                if (Math.abs(stateMatrix.get(ti, k) - targetState.get(k)) > 1e-10) {
                    matches = false;
                    break;
                }
            }

            if (matches) {
                timeInState += dt;
            }
        }

        return new ProbabilityResult(timeInState / totalTime);
    }

    /**
     * Estimates state probability using node index.
     *
     * @param nodeIndex Index of the stateful node
     * @param state The state vector to compute probability for
     * @return Estimated probability of the specified state
     */
    @Override
    public ProbabilityResult getProb(int nodeIndex, Matrix state) {
        NetworkStruct sn = this.sn != null ? this.sn : this.model.getStruct(true);
        int isf = (int) sn.nodeToStateful.get(0, nodeIndex);
        jline.lang.nodes.StatefulNode node = this.model.getStatefulNodes().get(isf);
        return getProb(node, state);
    }

    /**
     * Estimates state probability using current network state.
     *
     * @param node The stateful node to analyze
     * @return Estimated probability of current network state at the node
     */
    public ProbabilityResult getProb(jline.lang.nodes.StatefulNode node) {
        return getProb(node, null);
    }

    /**
     * Estimates the steady-state probability of a specific aggregated (per-class) state at a node.
     *
     * This method estimates the probability of observing a specific per-class job distribution
     * (e.g., [2 jobs of class 1, 1 job of class 2]) at a station. States are aggregated over
     * service phases - only the number of jobs per class matters.
     *
     * @param node The stateful node to analyze
     * @param stateAggr The aggregated state vector (per-class job counts) to compute probability for.
     *                  If null, uses current network state aggregated over phases.
     * @return Estimated probability of the specified aggregated state (0 if state not observed)
     */
    public ProbabilityResult getProbAggr(jline.lang.nodes.StatefulNode node, Matrix stateAggr) {
        // Get aggregated sample path via transient simulation
        jline.io.Ret.SampleResult sampleResult = sampleAggr(node, this.options.samples);

        if (sampleResult == null || sampleResult.t == null || !(sampleResult.state instanceof Matrix)) {
            return new ProbabilityResult(0.0);
        }

        Matrix t = sampleResult.t;
        Matrix stateMatrix = (Matrix) sampleResult.state;

        int numTimePoints = t.getNumRows();
        if (numTimePoints < 2) {
            return new ProbabilityResult(0.0);
        }

        // Get target state
        Matrix targetState = stateAggr;
        if (targetState == null) {
            NetworkStruct sn = this.sn != null ? this.sn : this.model.getStruct(true);
            int nodeIdx = node.getNodeIndex();
            int isf = (int) sn.nodeToStateful.get(0, nodeIdx);
            Matrix nodeState = sn.state.get(isf);
            // Aggregate to per-class counts
            jline.lang.state.State.StateMarginalStatistics marginal =
                jline.lang.state.ToMarginal.toMarginal(sn, nodeIdx, nodeState, null, null, null, null, null);
            targetState = marginal.nir;
        }

        // Compute time-weighted probability: sum time spent in target state / total time
        double totalTime = t.get(numTimePoints - 1, 0) - t.get(0, 0);
        if (totalTime <= 0) {
            return new ProbabilityResult(0.0);
        }

        double timeInState = 0.0;
        int numClasses = stateMatrix.getNumCols();

        for (int ti = 0; ti < numTimePoints - 1; ti++) {
            double dt = t.get(ti + 1, 0) - t.get(ti, 0);

            // Check if current state matches target state
            boolean matches = true;
            for (int k = 0; k < numClasses && k < targetState.length(); k++) {
                if (Math.abs(stateMatrix.get(ti, k) - targetState.get(k)) > 1e-10) {
                    matches = false;
                    break;
                }
            }

            if (matches) {
                timeInState += dt;
            }
        }

        return new ProbabilityResult(timeInState / totalTime);
    }

    /**
     * Estimates aggregated state probability using node index.
     *
     * @param nodeIndex Index of the stateful node
     * @param stateAggr The aggregated state vector to compute probability for
     * @return Estimated probability of the specified aggregated state
     */
    @Override
    public ProbabilityResult getProbAggr(int nodeIndex, Matrix stateAggr) {
        NetworkStruct sn = this.sn != null ? this.sn : this.model.getStruct(true);
        int isf = (int) sn.nodeToStateful.get(0, nodeIndex);
        jline.lang.nodes.StatefulNode node = this.model.getStatefulNodes().get(isf);
        return getProbAggr(node, stateAggr);
    }

    /**
     * Estimates aggregated state probability using current network state.
     *
     * @param node The stateful node to analyze
     * @return Estimated probability of current aggregated network state at the node
     */
    public ProbabilityResult getProbAggr(jline.lang.nodes.StatefulNode node) {
        return getProbAggr(node, null);
    }

    /**
     * Estimates the joint steady-state probability of the entire system state via DES simulation.
     *
     * This method estimates the probability of observing the current system state
     * (combined state across all stateful nodes) using simulation-based estimation.
     * States include phase information from service distributions.
     *
     * @return Estimated joint probability of the current system state (0 if state not observed)
     */
    @SuppressWarnings("unchecked")
    @Override
    public ProbabilityResult getProbSys() {
        // Get system-wide sample path via transient simulation
        jline.io.Ret.SampleResult sampleResult = sampleSys(this.options.samples);

        if (sampleResult == null || sampleResult.t == null || !(sampleResult.state instanceof List)) {
            return new ProbabilityResult(0.0);
        }

        Matrix t = sampleResult.t;
        List<Matrix> stateList = (List<Matrix>) sampleResult.state;

        int numTimePoints = t.getNumRows();
        if (numTimePoints < 2 || stateList.isEmpty()) {
            return new ProbabilityResult(0.0);
        }

        NetworkStruct sn = this.sn != null ? this.sn : this.model.getStruct(true);

        // Build target state from current network state
        List<Matrix> targetStates = new ArrayList<>();
        for (int isf = 0; isf < sn.nstateful; isf++) {
            targetStates.add(sn.state.get(isf));
        }

        // Compute time-weighted probability: sum time spent in target state / total time
        double totalTime = t.get(numTimePoints - 1, 0) - t.get(0, 0);
        if (totalTime <= 0) {
            return new ProbabilityResult(0.0);
        }

        double timeInState = 0.0;

        for (int ti = 0; ti < numTimePoints - 1; ti++) {
            double dt = t.get(ti + 1, 0) - t.get(ti, 0);

            // Check if current joint state matches target joint state
            boolean matches = true;
            for (int isf = 0; isf < stateList.size() && isf < targetStates.size() && matches; isf++) {
                Matrix nodeState = stateList.get(isf);
                Matrix targetState = targetStates.get(isf);
                int numClasses = nodeState.getNumCols();

                for (int k = 0; k < numClasses && k < targetState.length(); k++) {
                    if (Math.abs(nodeState.get(ti, k) - targetState.get(k)) > 1e-10) {
                        matches = false;
                        break;
                    }
                }
            }

            if (matches) {
                timeInState += dt;
            }
        }

        return new ProbabilityResult(timeInState / totalTime);
    }

    /**
     * Estimates the joint steady-state probability of the entire aggregated system state.
     *
     * This method estimates the probability of observing the current system state
     * (combined per-class job counts across all stateful nodes) using simulation.
     * States are aggregated over service phases - only job counts per class matter.
     *
     * @return Estimated joint probability of the current aggregated system state (0 if state not observed)
     */
    @SuppressWarnings("unchecked")
    @Override
    public ProbabilityResult getProbSysAggr() {
        // Get aggregated system-wide sample path via transient simulation
        jline.io.Ret.SampleResult sampleResult = sampleSysAggr(this.options.samples);

        if (sampleResult == null || sampleResult.t == null || !(sampleResult.state instanceof List)) {
            return new ProbabilityResult(0.0);
        }

        Matrix t = sampleResult.t;
        List<Matrix> stateList = (List<Matrix>) sampleResult.state;

        int numTimePoints = t.getNumRows();
        if (numTimePoints < 2 || stateList.isEmpty()) {
            return new ProbabilityResult(0.0);
        }

        NetworkStruct sn = this.sn != null ? this.sn : this.model.getStruct(true);

        // Build target state from current network state (aggregated)
        List<Matrix> targetStates = new ArrayList<>();
        for (int isf = 0; isf < sn.nstateful; isf++) {
            int nodeIdx = ((Double) sn.statefulToNode.get(isf)).intValue();
            Matrix nodeState = sn.state.get(isf);
            jline.lang.state.State.StateMarginalStatistics marginal =
                jline.lang.state.ToMarginal.toMarginal(sn, nodeIdx, nodeState, null, null, null, null, null);
            targetStates.add(marginal.nir);
        }

        // Compute time-weighted probability: sum time spent in target state / total time
        double totalTime = t.get(numTimePoints - 1, 0) - t.get(0, 0);
        if (totalTime <= 0) {
            return new ProbabilityResult(0.0);
        }

        double timeInState = 0.0;

        for (int ti = 0; ti < numTimePoints - 1; ti++) {
            double dt = t.get(ti + 1, 0) - t.get(ti, 0);

            // Check if current joint state matches target joint state
            boolean matches = true;
            for (int isf = 0; isf < stateList.size() && isf < targetStates.size() && matches; isf++) {
                Matrix nodeState = stateList.get(isf);
                Matrix targetState = targetStates.get(isf);
                int numClasses = nodeState.getNumCols();

                for (int k = 0; k < numClasses && k < targetState.length(); k++) {
                    if (Math.abs(nodeState.get(ti, k) - targetState.get(k)) > 1e-10) {
                        matches = false;
                        break;
                    }
                }
            }

            if (matches) {
                timeInState += dt;
            }
        }

        return new ProbabilityResult(timeInState / totalTime);
    }

    // =====================================================
    // Transient CDF Methods
    // =====================================================

    /**
     * Returns cumulative distribution functions of response times during transient analysis.
     * Uses response time samples collected during DES simulation to compute empirical CDFs.
     *
     * @return DistributionResult containing empirical CDFs for each station-class pair
     */
    @Override
    public DistributionResult getTranCdfRespT() {
        return getTranCdfRespT(getAvgRespTHandles());
    }

    /**
     * Returns cumulative distribution functions of response times during transient analysis.
     *
     * @param R response time handles specifying which metrics to compute
     * @return DistributionResult containing empirical CDFs for each station-class pair
     */
    @Override
    public DistributionResult getTranCdfRespT(AvgHandle R) {
        // Ensure simulation has been run
        if (this.result == null || !(this.result instanceof DESResult)) {
            try {
                runAnalyzer();
            } catch (Exception e) {
                throw new RuntimeException("Failed to run DES simulation for transient CDF: " + e.getMessage(), e);
            }
        }

        DESResult desResult = (DESResult) this.result;
        NetworkStruct sn = this.sn != null ? this.sn : this.model.getStruct(true);

        DistributionResult distResult = new DistributionResult(sn.nstations, sn.nclasses, "response_time");

        if (desResult.respTimeSamples == null) {
            return distResult;
        }

        // Compute empirical CDF for each station-class pair
        for (int i = 0; i < sn.nstations; i++) {
            for (int r = 0; r < sn.nclasses; r++) {
                if (i < desResult.respTimeSamples.length && r < desResult.respTimeSamples[i].length) {
                    List<Double> samples = desResult.respTimeSamples[i][r];
                    if (samples != null && !samples.isEmpty()) {
                        Matrix cdf = computeEmpiricalCDF(samples);
                        distResult.setCdf(i, r, cdf);
                    }
                }
            }
        }

        return distResult;
    }

    /**
     * Returns cumulative distribution functions of passage times during transient analysis.
     * For DES, passage times are equivalent to response times in single-visit networks.
     *
     * @return DistributionResult containing empirical CDFs for passage times
     */
    @Override
    public DistributionResult getTranCdfPassT() {
        return getTranCdfPassT(getAvgRespTHandles());
    }

    /**
     * Returns cumulative distribution functions of passage times during transient analysis.
     *
     * @param R response time handles specifying which metrics to compute
     * @return DistributionResult containing empirical CDFs for passage times
     */
    @Override
    public DistributionResult getTranCdfPassT(AvgHandle R) {
        // For single-visit networks, passage time equals response time
        // Use the same implementation as getTranCdfRespT
        return getTranCdfRespT(R);
    }

    /**
     * Computes empirical CDF from a list of samples.
     *
     * @param samples list of response time observations
     * @return Matrix with columns [F, X] where F is CDF values and X is data values
     */
    private Matrix computeEmpiricalCDF(List<Double> samples) {
        if (samples == null || samples.isEmpty()) {
            return new Matrix(0, 2);
        }

        // Sort samples
        double[] sortedSamples = samples.stream().mapToDouble(Double::doubleValue).sorted().toArray();
        int n = sortedSamples.length;

        // Build unique values with their CDF
        List<Double> uniqueX = new ArrayList<>();
        List<Double> cdfF = new ArrayList<>();

        double prevVal = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n; i++) {
            double val = sortedSamples[i];
            if (val != prevVal) {
                uniqueX.add(val);
                cdfF.add((i + 1.0) / n);
                prevVal = val;
            } else {
                // Update the CDF value for this duplicate
                cdfF.set(cdfF.size() - 1, (i + 1.0) / n);
            }
        }

        // Create result matrix [F, X]
        Matrix result = new Matrix(uniqueX.size(), 2);
        for (int i = 0; i < uniqueX.size(); i++) {
            result.set(i, 0, cdfF.get(i));  // F (CDF value)
            result.set(i, 1, uniqueX.get(i));  // X (data value)
        }

        return result;
    }

    // =====================================================
    // Transient Probability Methods
    // =====================================================

    /**
     * Computes transient state probabilities at a specific node over time using DES simulation.
     *
     * @param node the stateful node to analyze
     * @return ProbabilityResult containing transient probability data
     */
    public ProbabilityResult getTranProb(StatefulNode node) {
        // Ensure finite timespan
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProb in SolverDES requires a finite timespan. " +
                    "Use: SolverDES(model, \"timespan\", new double[]{0, T})");
        }

        // Run transient simulation
        try {
            getTranAvg();
        } catch (Exception e) {
            throw new RuntimeException("Failed to run transient DES simulation: " + e.getMessage(), e);
        }

        // Get sample path for this node
        jline.io.Ret.SampleResult sampleResult = sample(node, this.options.samples);

        if (sampleResult == null || sampleResult.t == null) {
            return new ProbabilityResult();
        }

        // Return probability result with the sampled state trajectory stored in probability field
        ProbabilityResult result = new ProbabilityResult();
        if (sampleResult.state instanceof Matrix) {
            result.probability = (Matrix) sampleResult.state;
        }
        result.state = sampleResult.t;
        return result;
    }

    /**
     * Computes transient aggregated state probabilities at a specific node over time.
     *
     * @param node the stateful node to analyze
     * @return ProbabilityResult containing transient aggregated probability data
     */
    public ProbabilityResult getTranProbAggr(StatefulNode node) {
        // Ensure finite timespan
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProbAggr in SolverDES requires a finite timespan. " +
                    "Use: SolverDES(model, \"timespan\", new double[]{0, T})");
        }

        // Run transient simulation
        try {
            getTranAvg();
        } catch (Exception e) {
            throw new RuntimeException("Failed to run transient DES simulation: " + e.getMessage(), e);
        }

        // Get aggregated sample path for this node
        jline.io.Ret.SampleResult sampleResult = sampleAggr(node, this.options.samples);

        if (sampleResult == null || sampleResult.t == null) {
            return new ProbabilityResult();
        }

        // Return probability result with the sampled state trajectory
        ProbabilityResult result = new ProbabilityResult();
        if (sampleResult.state instanceof Matrix) {
            result.probability = (Matrix) sampleResult.state;
        }
        result.state = sampleResult.t;
        result.isAggregated = true;
        return result;
    }

    /**
     * Computes transient system-wide state probabilities over time using DES simulation.
     *
     * @return ProbabilityResult containing transient system probability data
     */
    public ProbabilityResult getTranProbSys() {
        // Ensure finite timespan
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProbSys in SolverDES requires a finite timespan. " +
                    "Use: SolverDES(model, \"timespan\", new double[]{0, T})");
        }

        // Run transient simulation
        try {
            getTranAvg();
        } catch (Exception e) {
            throw new RuntimeException("Failed to run transient DES simulation: " + e.getMessage(), e);
        }

        // Get system-wide sample path
        jline.io.Ret.SampleResult sampleResult = sampleSys(this.options.samples);

        if (sampleResult == null || sampleResult.t == null) {
            return new ProbabilityResult();
        }

        // Return probability result with time stored in state field
        ProbabilityResult result = new ProbabilityResult();
        result.state = sampleResult.t;
        return result;
    }

    /**
     * Computes transient system-wide aggregated state probabilities over time.
     *
     * @return ProbabilityResult containing transient system aggregated probability data
     */
    public ProbabilityResult getTranProbSysAggr() {
        // Ensure finite timespan
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProbSysAggr in SolverDES requires a finite timespan. " +
                    "Use: SolverDES(model, \"timespan\", new double[]{0, T})");
        }

        // Run transient simulation
        try {
            getTranAvg();
        } catch (Exception e) {
            throw new RuntimeException("Failed to run transient DES simulation: " + e.getMessage(), e);
        }

        // Get aggregated system-wide sample path
        jline.io.Ret.SampleResult sampleResult = sampleSysAggr(this.options.samples);

        if (sampleResult == null || sampleResult.t == null) {
            return new ProbabilityResult();
        }

        // Return probability result with time stored in state field
        ProbabilityResult result = new ProbabilityResult();
        result.state = sampleResult.t;
        result.isAggregated = true;
        return result;
    }
}
