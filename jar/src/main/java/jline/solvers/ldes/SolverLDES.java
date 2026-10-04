/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import jline.GlobalConstants;
import jline.io.M2M;
import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SolverType;
import jline.lang.constant.NodeType;
import jline.lang.nodes.Cache;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.LayeredNetworkElement;
import jline.lang.layered.LayeredNetworkStruct;
import jline.solvers.AvgHandle;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.NetworkAvgCacheTable;
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
import jline.lang.reward.RewardFunction;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import jline.io.LDESResultIO;
import jline.io.LineModelIO;

import javax.xml.parsers.ParserConfigurationException;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput;
import static jline.io.InputOutput.line_debug;
import static jline.solvers.ldes.analyzers.Solver_ldes_analyzer.solver_ldes_analyzer;
import static jline.solvers.ldes.analyzers.Solver_ldes_analyzer_parallel.solver_ldes_analyzer_parallel;
import static jline.solvers.ldes.analyzers.Solver_ldes_ln_analyzer.solver_ldes_ln_analyzer;

/**
 * @brief LINE Discrete Event Simulator (LDES) solver using SSJ library.
 *
 * @details SolverLDES implements a discrete-event simulation solver that uses the SSJ
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
 * - **Cache**: List-based caches (LRU/FIFO/RR/SFIFO), optionally with a delayed-hit
 *   retrieval system (misses fetched through retrieval queues; concurrent requests
 *   for an in-flight item become delayed hits). LDES measures the retrieval latency
 *   directly, which the CTMC and SSA solvers report as NaN.
 *
 * @subsection des_sched Scheduling Strategies
 * - **FCFS/FCFSPR/FCFSPI**: First-Come-First-Served (non-preemptive, preemptive-resume, preemptive-independent)
 * - **LCFS/LCFSPR/LCFSPI**: Last-Come-First-Served (non-preemptive, preemptive-resume, preemptive-independent)
 * - **PS/DPS/GPS**: Processor Sharing variants (standard, discriminatory, generalized)
 * - **LPS**: Limited Processor Sharing
 * - **HOL**: Head-of-Line priority scheduling
 * - **SIRO**: Service In Random Order
 * - **SJF/LJF**: Shortest/Longest Job First
 * - **SEPT/LEPT**: Shortest/Longest Expected Processing Time
 * - **SRPT/SRPTPRIO**: Shortest Remaining Processing Time (preemptive)
 * - **PSJF/FB/LRPT/SETF**: size/age-based disciplines
 * - **FSP**: Fair Sojourn Protocol (preemptive; ranks by virtual PS finish time, dominates PS per-job)
 * - **EDD/EDF**: Earliest Due Date / Earliest Deadline First
 * - **PAS**: Pass-and-swap / order-independent queue (state-dependent total rate mu(c) with swapping graph)
 * - **POLLING**: Polling server (GATED, EXHAUSTIVE, KLIMITED)
 * - **EXT**: External arrivals (Source section)
 * - **INF**: Infinite server (for Delay nodes)
 * - Priority variants: FCFSPRIO, FCFSPRPRIO, FCFSPIPRIO, LCFSPRIO, LCFSPRPRIO, LCFSPIPRIO, PSPRIO, DPSPRIO, GPSPRIO
 *
 * @subsection des_routing Routing Strategies
 * - **PROB**: Probabilistic routing based on routing matrix
 * - **RAND**: Uniform random selection among destinations
 * - **RROBIN**: Round-robin cycling through destinations
 * - **WRROBIN**: Weighted round-robin
 * - **JSQ**: Join the Shortest Queue
 * - **SQ**: Power-of-K-choices (select shortest queue among K random samples)
 *
 * @subsection des_dist Service Time Distributions
 * - **Exp**: Exponential distribution
 * - **Erlang**: Erlang distribution (sum of exponentials)
 * - **HyperExp**: Hyperexponential (mixture of exponentials)
 * - **PH/APH**: Phase-type distributions (general and acyclic)
 * - **Coxian/Cox2**: Coxian distribution (general and 2-phase)
 * - **Det**: Deterministic (constant) service time
 * - **Uniform**: Uniform distribution
 * - **Gamma**: Gamma distribution
 * - **Pareto**: Pareto distribution
 * - **Weibull**: Weibull distribution
 * - **Lognormal**: Lognormal distribution
 * - **MAP/DMAP/MMAP/BMAP**: Markovian Arrival Processes (continuous, discrete, marked, batch)
 * - **MMPP2**: 2-state Markov-Modulated Poisson Process
 * - **ME/RAP**: Matrix-Exponential / Rational Arrival Process
 * - **NHPP**: Non-homogeneous Poisson process, piecewise-constant intensity (cyclic or not)
 * - **Immediate**: Zero service time (instantaneous)
 * - **Disabled**: No service (jobs bypass station)
 * - **Replayer/Trace**: Trace-driven service times from file
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
 * LDES supports stochastic Petri nets with:
 * - **Place nodes**: Token storage with capacity constraints
 * - **Transition nodes**: Token firing with enabling/inhibiting conditions
 * - **Timed transitions**: Exponential or phase-type firing delays
 * - **Immediate transitions**: Zero-delay firing with priority and weight
 *
 * @section des_lqn LayeredNetwork (LQN) Support
 * LDES simulates LayeredNetwork models directly from the LayeredNetworkStruct
 * (hosts, tasks, entries, activities, sync/async/forwarding calls, seq/AND/OR
 * precedences, two-phase entries). It also supports layered cache-queueing
 * models: a CacheTask holds a segmented cache of ItemEntry items with a
 * replacement policy (RR/FIFO/SFIFO/LRU/HLRU/CLIMB/QLRU); a read is a
 * synchronous or asynchronous call to an ItemEntry whose bound activity carries
 * a POST_CACHE precedence. On each read the simulator draws an item from the
 * item popularity distribution, tests it against the live cache content and
 * updates the content per the replacement policy, then routes the request to the
 * hit or the miss continuation activity. Hits/misses are exact sample-path
 * events, in contrast to the characteristic-time approximation used by the
 * analytic LN(MVA)/LN(NC) cache decomposition.
 *
 * @section des_analysis Analysis Modes
 * - **Steady-state**: Long-run average performance metrics (default)
 * - **Transient**: Time-varying metrics over specified time horizon
 *
 * @section des_warmstart Warm Start
 * Passing an auxiliary solver to the constructor, e.g.
 * {@code new SolverLDES(model, new SolverMVA(model))}, or calling
 * {@link #initFromSolver(NetworkSolver)}, warm-starts the simulation from that
 * solver's steady-state solution: SolverCTMC yields the mode of the exact
 * aggregate stationary distribution, any other solver a rounded mean
 * queue-length placement conserving the closed populations. The placement is
 * stored in options.init_sol and the transient filter is disabled, removing
 * the initialization bias (see jline.examples.java.advanced.LDESWarmStartExample).
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
 * SolverLDES solver = new SolverLDES(model, "seed", 12345, "samples", 100000);
 * SolverResult result = solver.getAvg();
 * @endcode
 *
 * @see LDESResult Result container for LDES metrics
 * @see LDESOptions Configuration options for LDES solver
 * @see NetworkSolver Base class for network solvers
 * @since 1.0
 * @author QORE Lab, Imperial College London
 */
public class SolverLDES extends NetworkSolver {

    /** LayeredNetwork model (non-null if solving LQN) */
    private LayeredNetwork lnModel;

    /** LayeredNetwork structure (populated if solving LQN) */
    private LayeredNetworkStruct lsn;

    /** Streaming collector for pushing metrics via OTLP */
    private Collector stream;

    /** Thread pool for parallel replication execution */
    public ExecutorService threadPool;

    /** Number of parallel threads for replication execution */
    private int numThreads;

    /**
     * Constructs a SolverLDES with the given model using default options.
     *
     * @param model the queueing network model to solve
     */
    public SolverLDES(Network model) {
        this(model, new LDESOptions());
        this.result = new LDESResult();
    }

    /**
     * Constructs a SolverLDES with the given LayeredNetwork model using default options.
     *
     * @param model the layered queueing network model to solve
     */
    public SolverLDES(LayeredNetwork model) {
        this(model, new LDESOptions());
    }

    /**
     * Constructs a SolverLDES with the given LayeredNetwork model and variable arguments.
     *
     * @param model the layered queueing network model to solve
     * @param args variable arguments for solver options
     */
    public SolverLDES(LayeredNetwork model, Object... args) {
        super(null, "SolverLDES");
        this.lnModel = model;
        this.lsn = model.getStruct(true);
        this.setOptions(Solver.parseOptions(new LDESOptions(), args));
        this.result = new LNLDESResult();
    }

    /**
     * Constructs a SolverLDES with the given LayeredNetwork model and options.
     *
     * @param model the layered queueing network model to solve
     * @param options solver configuration options
     */
    public SolverLDES(LayeredNetwork model, SolverOptions options) {
        super(null, "SolverLDES", options);
        this.lnModel = model;
        this.lsn = model.getStruct(true);
        this.result = new LNLDESResult();
    }

    /**
     * Constructs a SolverLDES with the given model and variable arguments.
     *
     * @param model the queueing network model to solve
     * @param args variable arguments for solver options
     */
    public SolverLDES(Network model, Object... args) {
        super(model, "SolverLDES");
        this.setOptions(Solver.parseOptions(new LDESOptions(), args));
        this.result = new LDESResult();
    }

    /**
     * Constructs a SolverLDES with the given model and method.
     *
     * @param model the queueing network model to solve
     * @param method the solution method to use
     */
    public SolverLDES(Network model, String method) {
        super(model, "SolverLDES", new LDESOptions().method(method));
        this.result = new LDESResult();
    }

    /**
     * Constructs a SolverLDES with the given model and options.
     *
     * @param model the queueing network model to solve
     * @param options solver configuration options
     */
    public SolverLDES(Network model, SolverOptions options) {
        super(model, "SolverLDES", options);
        this.result = new LDESResult();
    }

    /**
     * Load a JMT model file (.jsimg/.jsim/.jsimw/.jmva) and simulate with LDES.
     *
     * @param filename path to the JMT model file
     */
    public SolverLDES(String filename) {
        this(new M2M().JMT2LINE(filename));
    }

    /**
     * Load a JMT model file (.jsimg/.jsim/.jsimw/.jmva) and simulate with LDES,
     * accepting variable option arguments (e.g., "seed", 23000, "samples", 50000).
     *
     * @param filename path to the JMT model file
     * @param args variable arguments for solver options
     */
    public SolverLDES(String filename, Object... args) {
        this(new M2M().JMT2LINE(filename), args);
    }

    /**
     * Constructs a SolverLDES that warm-starts the simulation from the steady-state
     * distribution computed by another solver. The auxiliary solver is used to find the
     * stationary distribution of the model; the resulting distribution decides the
     * initial state of the simulation (see {@link #initFromSolver(NetworkSolver)}).
     *
     * @param model the queueing network model to solve
     * @param initSolver auxiliary solver used to compute the steady-state distribution
     * @param args variable arguments for solver options
     */
    public SolverLDES(Network model, NetworkSolver initSolver, Object... args) {
        this(model, args);
        this.initFromSolver(initSolver);
    }

    /**
     * Warm-start the simulation from the steady-state solution of an auxiliary solver.
     *
     * If the auxiliary solver is a {@link SolverCTMC}, the exact stationary distribution
     * over the aggregate state space is computed and the initial state is set to the
     * mode of that distribution (the most probable aggregate state). For any other
     * network solver, the steady-state mean queue lengths are used instead and rounded
     * to an integer placement that conserves each closed-class population.
     *
     * Since the simulation starts (approximately) in steady state, the transient
     * removal filter is disabled (tranfilter = "fixed" with warmupfrac = 0), so every
     * simulated sample contributes to the estimators.
     *
     * @param initSolver auxiliary solver used to compute the steady-state distribution
     * @return this solver, for chaining
     */
    @Override
    public SolverLDES initFromSolver(NetworkSolver initSolver) {
        NetworkStruct snl = this.getStruct();
        int M = snl.nstations;
        int K = snl.nclasses;
        Matrix placement = jline.solvers.WarmStart.warmStartPlacement(initSolver, snl);

        // LDES consumes the placement through options.init_sol (station-major
        // vector) rather than the model initial state.
        Matrix initSol = new Matrix(1, M * K);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                initSol.set(0, i * K + r, placement.get(i, r));
            }
        }
        this.options.init_sol = initSol;
        if (this.options instanceof LDESOptions) {
            LDESOptions ldesOptions = (LDESOptions) this.options;
            ldesOptions.tranfilter = "fixed";
            ldesOptions.warmupfrac = 0.0;
        }
        return this;
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
     * Returns the feature set supported by the LDES solver.
     * Supports multiclass Jackson queueing networks with FCFS queues.
     *
     * @return the feature set supported by the LDES solver
     */
    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(new String[]{
                "Sink", "Source",
                "Queue", "Delay",
                "Fork", "Join", "Forker", "Joiner",  // Fork-Join node support
                "Place", "Transition",  // Stochastic Petri Net support
                "QueueingPlace",  // Queueing place (QPN embedded queue): FCFS/LCFS/SIRO/INF, renewal service
                "Linkage", "Enabling", "Inhibiting", "Timing", "Firing", "Storage",  // Petri net section support
                "Logger", "LogTunnel",  // Logger node support
                "Cache", "CacheClassSwitcher", "CacheRetrieval",  // Cache node support with hit/miss class switching
                "ReplacementStrategy_LRU", "ReplacementStrategy_FIFO", "ReplacementStrategy_RR", "ReplacementStrategy_SFIFO",  // Cache replacement strategies
                "ReplacementStrategy_HLRU", "ReplacementStrategy_CLIMB", "ReplacementStrategy_QLRU",
                "Buffer",  // Finite buffer capacity support
                "Region",  // Finite capacity region support
                "Exp", "Erlang", "HyperExp", "PH", "APH", "Coxian", "Cox2", "MAP", "DMAP", "MMAP", "BMAP", "MMPP2", "ME", "RAP", "Immediate", "Disabled", "Replayer", "Trace",  // Trace is an alias of Replayer
                "Det", "Uniform", "Gamma", "Pareto", "Weibull", "Lognormal",  // Additional continuous distributions
                "Geometric",  // Lattice-valued interarrival/service time on {1,2,...} (Geo/Geo/1 and slotted models)
                "Bernoulli", "Binomial", "Poisson",  // Counting distributions; their zero atom becomes an immediate interval (continuous mode only)
                "NHPP",  // Piecewise-constant-intensity non-homogeneous Poisson process
                "Server", "JobSink", "RandomSource",
                "InfiniteServer", "SharedServer", "ServiceTunnel", "DelayStation",  // internal station-section markers
                "SchedStrategy_FCFS", "SchedStrategy_INF",
                "SchedStrategy_HOL",  // Priority scheduling (FCFS with priorities)
                "SchedStrategy_FCFSPRIO",  // FCFS with priorities (non-preemptive), handled via priorityComparator
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
                "SchedStrategy_FSP",  // Fair Sojourn Protocol (virtual PS finish time ranking)
                "SchedStrategy_EDD",  // Earliest Due Date (non-preemptive, deadline ordering)
                "SchedStrategy_EDF",  // Earliest Deadline First (preemptive-resume by deadline)
                "SchedStrategy_SETF",  // Shortest Elapsed Time First (attained-service ordering)
                "SchedStrategy_FCFSPR",    // FCFS Preemptive Resume
                "SchedStrategy_FCFSPI",    // FCFS Preemptive Independent
                "SchedStrategy_FCFSPRPRIO",   // FCFSPR with priorities
                "SchedStrategy_FCFSPIPRIO",   // FCFSPI with priorities
                "SchedStrategy_LPS",       // Limited Processor Sharing
                "SchedStrategy_EXT",
                "SchedStrategy_POLLING",  // Polling scheduling (GATED, EXHAUSTIVE, KLIMITED)
                "SchedStrategy_PAS",      // Pass-and-swap (order-independent) queue
                "SchedStrategy_OI",       // Order-independent queue (PAS specialization, empty swap graph)
                "Router", "Dispatcher",  // Router node support (Dispatcher is the internal router section)
                "ClassSwitch", "StatelessClassSwitcher",  // Class switching node support
                "RoutingStrategy_PROB", "RoutingStrategy_RAND",
                "RoutingStrategy_RROBIN", "RoutingStrategy_WRROBIN",
                "RoutingStrategy_JSQ",  // Join Shortest Queue routing
                "RoutingStrategy_SQ",  // Power of K Choices routing
                "OpenClass",
                "ClosedClass",
                "SelfLoopingClass",
                "OpenSignal",           // G-network signal class in open networks
                "ClosedSignal",         // G-network signal class in closed networks
                "SignalType_NEGATIVE",
                "SignalType_REPLY",
                "SignalType_CATASTROPHE",
                "BatchArrival",          // Engine reads sn.arrivalbatch (Geo^X and other batch arrival streams)
                "SignalBatchRemoval",    // Engine reads sn.signalremdist
                "SignalRemovalPolicy",   // Engine reads sn.signalrempolicy
                "LoadDependence",  // Load-dependent service rates
                "ClassDependence", // Class-dependent service rate handles (setLimitedClassDependence)
                "JointDependence", // Joint-dependent (non-product-form) service rate handles (setJointDependence)
                "SetupDelayOff",   // Engine simulates the SETUP/DELAYOFF server states
                "Balking",         // Engine reads sn.balkingStrategy / balkingThresholds
                "Reneging",        // Engine collects renegingRate / avgRenegingWaitTime
                "Retrial"          // Engine collects retrialDropped and successful retries
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
     * The number of replications is determined by LDESOptions.replications.
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
                throw new RuntimeException("LDES LN analysis failed: " + e.getMessage(), e);
            }
            return this.result;
        }
        // For Network models, use parent implementation
        return super.getAvg();
    }

    /**
     * Returns the list of valid solution methods for LDES.
     *
     * @return list of valid method names
     */
    public List<String> listValidMethods() {
        return listValidMethods(null);
    }

    /**
     * Returns the list of valid solution methods for LDES.
     *
     * @param model the network model (unused, for interface compatibility)
     * @return list of valid method names
     */
    public List<String> listValidMethods(Network model) {
        return Arrays.asList("default", "parallel");
    }

    /**
     * Simulates the model on a remote LDES REST server instead of in this JVM.
     *
     * The exchange is the same pair of documents the command line uses: the
     * model is serialized with {@link LineModelIO} and POSTed to
     * &lt;restUrl&gt;/api/v1/solve, and the ldes-result document that comes back is
     * read with {@link LDESResultIO}. The network struct is not part of that
     * document, so the local one is attached to the result.
     *
     * @param restUrl base URL of the server, e.g. "http://localhost:8080"
     * @param ldesOptions the LDES options, or null if only base options are set
     * @return the deserialized simulation result
     * @throws IOException if the model cannot be serialized or the server fails
     */
    private LDESResult solveRemote(String restUrl, LDESOptions ldesOptions) throws IOException {
        String url = restUrl.replaceAll("/+$", "");
        if (!url.matches(".*/api/v\\d+/solve$")) {
            url = url + "/api/v1/solve";
        }

        java.nio.file.Path work = java.nio.file.Files.createTempDirectory("jline_ldes_rest_");
        try {
            java.nio.file.Path modelPath = work.resolve("model.json");
            LineModelIO.save(this.model, modelPath.toString());
            String modelText = new String(java.nio.file.Files.readAllBytes(modelPath),
                    java.nio.charset.StandardCharsets.UTF_8);

            JsonObject modelNode = new JsonObject();
            modelNode.addProperty("content", modelText);
            modelNode.addProperty("base64", false);
            JsonArray flagNode = new JsonArray();
            for (String flag : buildRemoteFlags(ldesOptions)) {
                flagNode.add(flag);
            }
            JsonObject request = new JsonObject();
            request.add("model", modelNode);
            request.add("flags", flagNode);

            line_debug(options.verbose, "LDES REST: POST " + url + " " + flagNode);
            JsonObject response = postJson(url, request.toString(), remoteTimeoutMillis());

            String status = response.has("status") ? response.get("status").getAsString() : "";
            if (!"ok".equals(status) || !response.has("result")) {
                String message = response.has("message") ? response.get("message").getAsString()
                        : "unspecified error";
                String stderr = response.has("stderr") ? response.get("stderr").getAsString() : "";
                throw new IOException("LDES REST solve failed: " + message
                        + (stderr.isEmpty() ? "" : (" Engine stderr: " + stderr.trim())));
            }

            java.nio.file.Path resultPath = work.resolve("result.json");
            java.nio.file.Files.write(resultPath,
                    response.get("result").toString().getBytes(java.nio.charset.StandardCharsets.UTF_8));
            LDESResult result = LDESResultIO.load(resultPath.toString());
            // The result document carries metrics only; the struct stays local.
            result.sn = this.sn;
            return result;
        } finally {
            deleteRecursively(work);
        }
    }

    /**
     * Maps the solver options onto the LDES long-form flags that follow
     * "solve &lt;model&gt; -o &lt;result&gt;" on the command line. The REST server accepts
     * the same tokens, so the remote and local runs are configured identically.
     *
     * @param o the LDES options, or null if only base options are set
     * @return the flag tokens, in command-line order
     */
    private List<String> buildRemoteFlags(LDESOptions o) {
        List<String> flags = new ArrayList<String>();
        flags.add("--samples");
        flags.add(String.valueOf(options.samples));
        // Always emitted: the CLI default is a random seed, so omitting it would
        // make the remote run irreproducible against the local one.
        flags.add("--seed");
        flags.add(String.valueOf(options.seed));
        if (options.method != null && !"default".equals(options.method)) {
            flags.add("--method");
            flags.add(options.method);
        }
        if (o == null) {
            return flags;
        }
        if (o.cnvgon) {
            flags.add("--cnvgon");
            flags.add("--cnvgtol");
            flags.add(String.valueOf(o.cnvgtol));
        }
        if (!"mser5".equals(o.tranfilter)) {
            flags.add("--tranfilter");
            flags.add(o.tranfilter);
        }
        if (o.warmupfrac != LDESOptions.DEFAULT_WARMUP_FRAC) {
            flags.add("--warmupfrac");
            flags.add(String.valueOf(o.warmupfrac));
        }
        if (!"obm".equals(o.cimethod)) {
            flags.add("--cimethod");
            flags.add(o.cimethod);
        }
        if (o.slotted) {
            flags.add("--slotted");
            if (o.slotLength != LDESOptions.DEFAULT_SLOT_LENGTH) {
                flags.add("--slotlength");
                flags.add(String.valueOf(o.slotLength));
            }
        }
        if (o.replications > 1) {
            flags.add("--replications");
            flags.add(String.valueOf(o.replications));
            flags.add("--numthreads");
            flags.add(String.valueOf(o.numThreads));
        }
        if (o.maxSimEvents > 0) {
            flags.add("--maxevents");
            flags.add(String.valueOf(o.maxSimEvents));
        }
        if (!Double.isInfinite(o.maxTime) && o.maxTime > 0) {
            flags.add("--maxtime");
            flags.add(String.valueOf(o.maxTime));
        }
        return flags;
    }

    /**
     * Socket timeout for a remote solve: the cooperative --maxtime budget plus
     * grace for the server to serialize its result, or ten minutes if unbounded.
     *
     * @return the timeout in milliseconds
     */
    private int remoteTimeoutMillis() {
        double budget = options.timeout;
        if (!Double.isInfinite(budget) && budget > 0) {
            return (int) Math.min(Integer.MAX_VALUE, (budget + 30.0) * 1000.0);
        }
        return 600000;
    }

    /**
     * POSTs a JSON document and returns the parsed JSON response, whether the
     * server answered 2xx or an error status (the error body is itself JSON and
     * carries the engine diagnostics).
     *
     * @param url the endpoint
     * @param body the request document
     * @param timeoutMillis connect and read timeout
     * @return the parsed response document
     * @throws IOException if the exchange fails or the response is not JSON
     */
    private static JsonObject postJson(String url, String body, int timeoutMillis) throws IOException {
        HttpURLConnection conn = (HttpURLConnection) new java.net.URL(url).openConnection();
        try {
            conn.setRequestMethod("POST");
            conn.setDoOutput(true);
            conn.setConnectTimeout(Math.min(timeoutMillis, 30000));
            conn.setReadTimeout(timeoutMillis);
            conn.setRequestProperty("Content-Type", "application/json; charset=utf-8");
            byte[] payload = body.getBytes(java.nio.charset.StandardCharsets.UTF_8);
            conn.setFixedLengthStreamingMode(payload.length);
            OutputStream os = conn.getOutputStream();
            try {
                os.write(payload);
            } finally {
                os.close();
            }

            InputStream in = (conn.getResponseCode() >= 400) ? conn.getErrorStream() : conn.getInputStream();
            if (in == null) {
                throw new IOException("LDES REST server returned HTTP " + conn.getResponseCode()
                        + " with no body");
            }
            StringBuilder sb = new StringBuilder();
            BufferedReader reader = new BufferedReader(
                    new InputStreamReader(in, java.nio.charset.StandardCharsets.UTF_8));
            try {
                String line;
                while ((line = reader.readLine()) != null) {
                    sb.append(line);
                }
            } finally {
                reader.close();
            }
            try {
                return JsonParser.parseString(sb.toString()).getAsJsonObject();
            } catch (RuntimeException e) {
                throw new IOException("LDES REST server returned a non-JSON body: " + sb);
            }
        } finally {
            conn.disconnect();
        }
    }

    /**
     * Removes a temporary directory and its contents, ignoring failures: a
     * leaked scratch directory is not worth failing a completed solve.
     *
     * @param root the directory to remove
     */
    private static void deleteRecursively(java.nio.file.Path root) {
        if (root == null) {
            return;
        }
        try {
            java.util.stream.Stream<java.nio.file.Path> walk = java.nio.file.Files.walk(root);
            try {
                walk.sorted(java.util.Comparator.reverseOrder()).forEach(p -> {
                    try {
                        java.nio.file.Files.deleteIfExists(p);
                    } catch (IOException ignored) {
                    }
                });
            } finally {
                walk.close();
            }
        } catch (IOException ignored) {
        }
    }

    @Override
    public void runAnalyzer() throws IllegalAccessException, ParserConfigurationException, IOException {
        long T0 = System.nanoTime();
        if (this.options == null) {
            this.options = new SolverOptions(SolverType.LDES);
        }
        // options.events (DES event budget) overrides options.samples when set;
        // samples remains accepted as a deprecated alias for the event budget.
        if (this.options.events > 0) {
            this.options.samples = this.options.events;
        }
        // Map the generic wall-clock budget (options.timeout, seconds) onto the
        // LDES-specific maxTime (enforced by the SSJ event loop) unless maxTime was
        // set explicitly. Handles the case where options is a base SolverOptions
        // (timeout only) rather than an LDESOptions.
        if (this.options instanceof LDESOptions) {
            LDESOptions _lo = (LDESOptions) this.options;
            if (Double.isInfinite(_lo.maxTime) && Double.isFinite(_lo.timeout) && _lo.timeout > 0) {
                _lo.maxTime = _lo.timeout;
            }
        }
        // Propagate solver verbose level to global
        GlobalConstants.Verbose = options.verbose;

        // Check if we're solving a LayeredNetwork model
        if (this.lnModel != null) {
            runLNAnalyzer();
            return;
        }

        if (this.enableChecks && !supports(this.model)) {
            throw new RuntimeException("This model is not supported by the LDES solver.");
        }
        this.resetRandomGeneratorSeed(options.seed);

        // Clear any cache hit/miss/delayed/latency split left on the cache nodes by
        // a previously run solver. The struct rebuild bakes a populated actualHitProb
        // into the cache routing (Network.getRoutingMatrix), which the simulator then
        // reads; a stale value from another solver (e.g. SSA) would otherwise make the
        // cache node metrics order-dependent. The simulator repopulates its own split.
        // Mirrors the equivalent clear in SolverCTMC.
        for (int ci = 0; ci < this.model.getNodes().size(); ci++) {
            if (this.model.getNodes().get(ci) instanceof Cache) {
                Cache cacheNode = (Cache) this.model.getNodes().get(ci);
                cacheNode.setResultHitProb(new Matrix(0, 0));
                cacheNode.setResultMissProb(new Matrix(0, 0));
                cacheNode.setResultDelayedHitProb(new Matrix(0, 0));
                cacheNode.setResultResidT(new Matrix(0, 0));
                cacheNode.setResultHitProbList(new Matrix(0, 0));
                cacheNode.setResultItemProb(new Matrix(0, 0));
            }
        }
        String method = options.method;

        // Check for parallel replication mode (only if options is LDESOptions)
        LDESOptions ldesOptions = null;
        boolean useParallel = false;
        if (this.options instanceof LDESOptions) {
            ldesOptions = (LDESOptions) this.options;
            useParallel = ldesOptions.replications > 1;
        }

        if (useParallel) {
            // Initialize thread pool for parallel replications
            this.numThreads = Math.min(ldesOptions.numThreads, ldesOptions.replications);
            this.threadPool = Executors.newFixedThreadPool(this.numThreads);
            line_debug(options.verbose, String.format("LDES solver starting (parallel): replications=%d, threads=%d, samples=%d, seed=%d",
                    ldesOptions.replications, this.numThreads, options.samples, options.seed));
        } else {
            line_debug(options.verbose, String.format("LDES solver starting: method=%s, samples=%d, seed=%d",
                    method, options.samples, options.seed));
        }
        line_debug(options.verbose, "Running LDES simulation, calling solver_ldes_analyzer");

        // Remote engine, if options.restUrl names an LDES REST server: the
        // container simulates the model and returns the same ldes-result
        // document the CLI writes, so the metrics are identical for a fixed
        // seed and event budget.
        String restUrl = (ldesOptions != null) ? ldesOptions.restUrl : null;

        LDESResult result;
        try {
            if (restUrl != null && !restUrl.isEmpty()) {
                result = solveRemote(restUrl, ldesOptions);
            } else if (useParallel) {
                result = solver_ldes_analyzer_parallel(this.sn, ldesOptions, this);
            } else {
                result = solver_ldes_analyzer(this.sn, this.options, this);
            }
        } catch (RuntimeException e) {
            // Carry the cause's message forward: the engine rejects unsupported
            // configurations (e.g. heterogeneous-server NHPP service, non-cyclic
            // NHPP in steady state) with an actionable diagnostic, and the
            // outer getAvg wrapper reports only getMessage(), so a bare
            // "LDES simulation failed." would discard it.
            String detail = (e.getMessage() != null) ? e.getMessage() : e.toString();
            throw new RuntimeException("LDES simulation failed: " + detail, e);
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

        // Propagate cache hit/miss ratios (set on the Cache nodes by the
        // simulator) into this.sn.nodeparam so the node-level hit/miss
        // throughput table is computed from the actual probabilities rather
        // than the uniform nodevisits fallback.
        if (this.sn != null && this.sn.nodeparam != null) {
            for (int ind = 0; ind < this.sn.nnodes; ind++) {
                if (this.sn.nodetype.get(ind) == NodeType.Cache) {
                    Cache cacheNode = (Cache) this.model.getNodes().get(ind);
                    CacheNodeParam cnp = (CacheNodeParam) this.sn.nodeparam.get(cacheNode);
                    if (cnp != null) {
                        Matrix hp = cacheNode.getHitRatio();
                        Matrix mp = cacheNode.getMissRatio();
                        if (hp != null && !hp.isEmpty()) cnp.actualhitprob = hp;
                        if (mp != null && !mp.isEmpty()) cnp.actualmissprob = mp;
                    }
                }
            }
        }

        double runtime = result.runtime;
        int M = sn.nstations;
        int R = sn.nclasses;
        // Use arrival rate from LDES result (includes dropped jobs), fallback to computed if null
        Matrix AN = result.AN != null ? result.AN : snGetArvRFromTput(sn, TN, getAvgTputHandles());
        Matrix WN = new Matrix(0, 0);

        this.result.method = result.method;
        this.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN, runtime, this.result.method, options.samples);

        // Set tardiness results
        this.result.TardN = TardN;
        this.result.SysTardN = SysTardN;

        // Transfer confidence interval and convergence data if available
        if (result instanceof LDESResult) {
            LDESResult ldesResult = (LDESResult) result;
            ((LDESResult) this.result).QNCI = ldesResult.QNCI;
            ((LDESResult) this.result).UNCI = ldesResult.UNCI;
            ((LDESResult) this.result).RNCI = ldesResult.RNCI;
            ((LDESResult) this.result).TNCI = ldesResult.TNCI;
            ((LDESResult) this.result).ANCI = ldesResult.ANCI;
            ((LDESResult) this.result).WNCI = ldesResult.WNCI;
            // Transfer convergence data
            ((LDESResult) this.result).converged = ldesResult.converged;
            ((LDESResult) this.result).stoppingReason = ldesResult.stoppingReason;
            ((LDESResult) this.result).convergenceBatches = ldesResult.convergenceBatches;
            ((LDESResult) this.result).QNRelPrec = ldesResult.QNRelPrec;
            ((LDESResult) this.result).UNRelPrec = ldesResult.UNRelPrec;
            ((LDESResult) this.result).RNRelPrec = ldesResult.RNRelPrec;
            ((LDESResult) this.result).TNRelPrec = ldesResult.TNRelPrec;
            // Transfer impatience statistics
            ((LDESResult) this.result).renegedCustomers = ldesResult.renegedCustomers;
            ((LDESResult) this.result).avgRenegingWaitTime = ldesResult.avgRenegingWaitTime;
            ((LDESResult) this.result).renegingRate = ldesResult.renegingRate;
            ((LDESResult) this.result).balkedCustomers = ldesResult.balkedCustomers;
            ((LDESResult) this.result).balkingProbability = ldesResult.balkingProbability;
            ((LDESResult) this.result).retriedCustomers = ldesResult.retriedCustomers;
            ((LDESResult) this.result).retrialDropped = ldesResult.retrialDropped;
            ((LDESResult) this.result).avgOrbitSize = ldesResult.avgOrbitSize;
            // Transfer simulation statistics
            ((LDESResult) this.result).totalSimulatedEvents = ldesResult.totalSimulatedEvents;
            ((LDESResult) this.result).QNSamples = ldesResult.QNSamples;
            ((LDESResult) this.result).UNSamples = ldesResult.UNSamples;
            ((LDESResult) this.result).RNSamples = ldesResult.RNSamples;
            ((LDESResult) this.result).TNSamples = ldesResult.TNSamples;
            // Transfer Markov reward metrics
            ((LDESResult) this.result).avgReward = ldesResult.avgReward;
            ((LDESResult) this.result).tranReward = ldesResult.tranReward;
            ((LDESResult) this.result).rewardTime = ldesResult.rewardTime;
            ((LDESResult) this.result).rewardNames = ldesResult.rewardNames;
            ((LDESResult) this.result).stateHistogramSpace = ldesResult.stateHistogramSpace;
            ((LDESResult) this.result).stateHistogramTime = ldesResult.stateHistogramTime;
            // The per-job response time samples are the empirical CDF input of
            // getCdfRespT and of the --respt-samples JSON export. runAnalyzer
            // publishes a fresh LDESResult and copies fields onto it one by
            // one, so omitting this line left the samples in the handler's
            // result and both consumers saw null.
            ((LDESResult) this.result).respTimeSamples = ldesResult.respTimeSamples;
            ((LDESResult) this.result).stateTrajectorySpace = ldesResult.stateTrajectorySpace;
            ((LDESResult) this.result).stateTrajectoryTime = ldesResult.stateTrajectoryTime;
        }

        // Transfer FCR (Finite Capacity Region) metrics if available
        this.result.QNfcr = result.QNfcr;
        this.result.UNfcr = result.UNfcr;
        this.result.RNfcr = result.RNfcr;
        this.result.TNfcr = result.TNfcr;
        this.result.ANfcr = result.ANfcr;
        this.result.WNfcr = result.WNfcr;
        this.result.WeightNfcr = result.WeightNfcr;
        this.result.MemOccNfcr = result.MemOccNfcr;
        this.result.DropRateNfcr = result.DropRateNfcr;
        this.result.DropRateJoin = result.DropRateJoin;

        // Handle transient results if available
        if (result.QNt != null && result.QNt.length > 0) {
            // Create empty response time transient matrices (not computed by LDES)
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
     * Runs the LDES analyzer for LayeredNetwork models.
     */
    private void runLNAnalyzer() {
        this.resetRandomGeneratorSeed(options.seed);
        String method = options.method;
        line_debug(options.verbose, String.format("LDES LN solver starting: method=%s, samples=%d, seed=%d",
                method, options.samples, options.seed));
        line_debug(options.verbose, "Running LDES LN simulation, calling solver_ldes_ln_analyzer");

        LNLDESResult lnResult;
        try {
            lnResult = solver_ldes_ln_analyzer(this.lsn, this.options, this);
        } catch (RuntimeException e) {
            throw new RuntimeException("LDES LN simulation failed.", e);
        }

        // Store results
        this.result = lnResult;
        this.result.method = lnResult.method;
        this.result.runtime = lnResult.runtime;
    }

    /**
     * Run the LDES analyzer with current options including init_sol.
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

        // Call runAnalyzer which delegates to solver_ldes_analyzer
        try {
            runAnalyzer();
        } catch (Exception e) {
            throw new RuntimeException("LDES transient analysis failed", e);
        }

        return this.result;
    }

    @Override
    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverLDES.getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    /**
     * LDES is a discrete-event simulator; all methods are stochastic.
     *
     * @param method the method name to classify
     * @return true always
     */
    @Override
    public boolean isStochasticMethod(String method) {
        return true;
    }

    /**
     * Checks if the solver supports the given LayeredNetwork model.
     *
     * @param model the LayeredNetwork model to check
     * @return true if the model is supported
     */
    public boolean supports(LayeredNetwork model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverLDES.getLNFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    /**
     * Returns the feature set supported by the LDES solver for LayeredNetwork models.
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
                // Layered cache-queueing models (CacheTask + ItemEntry): the read is a
                // sync/async call to an ItemEntry resolved to a hit or miss branch by
                // simulating the replacement policy directly over the item popularity.
                "CacheTask", "ItemEntry", "Cache", "ActivityPrecedence_POST_CACHE",
                "ReplacementStrategy_RR", "ReplacementStrategy_FIFO", "ReplacementStrategy_SFIFO",
                "ReplacementStrategy_LRU", "ReplacementStrategy_HLRU", "ReplacementStrategy_CLIMB",
                "ReplacementStrategy_QLRU",
                "Exp", "Erlang", "HyperExp", "PH", "APH", "Coxian", "Det", "Uniform", "Gamma",
                "Lognormal", "Weibull", "Pareto"  // i.i.d. renewal host-demand distributions (generic sampler)
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
     * Cache performance table for a layered cache-queueing (LCQ) model: one row per
     * CacheTask read (ItemEntry), reporting hit / delayed-hit / miss probabilities and
     * rates. Mirrors the flat-network {@link jline.solvers.NetworkSolver#getAvgCacheTable}
     * convention (delayed hits reported in their own column; reads = hits + misses +
     * delayed hits). Returns an empty table for LQNs without a CacheTask.
     */
    public NetworkAvgCacheTable getAvgCacheTable() {
        if (this.lnModel == null) {
            throw new RuntimeException("getAvgCacheTable(LN) requires a LayeredNetwork model");
        }
        if (!(this.result instanceof LNLDESResult) || ((LNLDESResult) this.result).TLN == null) {
            getAvg();
        }
        LNLDESResult r = (LNLDESResult) this.result;

        List<Double> List_ = new ArrayList<Double>(), ListCap = new ArrayList<Double>(), Items = new ArrayList<Double>();
        List<Double> HitProb = new ArrayList<Double>(), DelayedHitProb = new ArrayList<Double>(), MissProb = new ArrayList<Double>();
        List<Double> HitRate = new ArrayList<Double>(), DelayedHitRate = new ArrayList<Double>(), MissRate = new ArrayList<Double>();
        List<Double> ArvR = new ArrayList<Double>(), ResidT = new ArrayList<Double>();
        List<String> nodeNames = new ArrayList<String>(), classNames = new ArrayList<String>();

        for (int k = 0; k < r.cacheTaskIdx.size(); k++) {
            int tidx = r.cacheTaskIdx.get(k);
            int eidx = r.cacheItemEntryIdx.get(k);
            double ph = r.cacheHitProb.get(k), pd = r.cacheDelayedProb.get(k), pm = r.cacheMissProb.get(k);
            double rate = r.cacheReadRate.get(k);
            int[] cap = (lsn.itemcap != null) ? lsn.itemcap.get(tidx) : null;
            double totcap = 0.0;
            if (cap != null) for (int c : cap) totcap += c;
            List_.add(0.0);
            ListCap.add(totcap);
            Items.add(lsn.nitems.get(0, tidx));
            HitProb.add(ph); DelayedHitProb.add(pd); MissProb.add(pm);
            HitRate.add(rate * ph); DelayedHitRate.add(rate * pd); MissRate.add(rate * pm);
            ArvR.add(rate * (pm + pd));   // rate of requests entering the retrieval system
            ResidT.add(0.0);
            nodeNames.add(lsn.names.get(tidx));
            classNames.add(lsn.names.get(eidx));
        }

        NetworkAvgCacheTable table = new NetworkAvgCacheTable(List_, ListCap, Items,
                HitProb, DelayedHitProb, MissProb, HitRate, DelayedHitRate, MissRate, ArvR, ResidT);
        table.setNodeNames(nodeNames);
        table.setClassNames(classNames);
        table.setOptions(this.options);
        return table;
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
        if (!(this.result instanceof LNLDESResult)) {
            throw new RuntimeException("getLNAvgTable requires LNLDESResult - call getAvg() first");
        }

        LNLDESResult lnResult = (LNLDESResult) this.result;

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
            // LINE convention: utilization of a non-infinite-server station is a
            // per-server fraction in [0,1]. The simulator accumulates mean busy
            // servers (0..multiplicity); rescale by the host processor multiplicity
            // to match the LN/LQNS convention (see SolverLQNS.hostProcessorMult).
            // Infinite-server processors return mult 1.0 here and keep mean-busy-servers.
            double uRaw = lnResult.ULN != null ? lnResult.ULN.get(0, i) : 0.0;
            double m = hostProcessorMult(this.lsn, i);
            Uval.add(m > 1.0 ? uRaw / m : uRaw);
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
     * Returns the multiplicity of the host processor on which the given LQN
     * element resides, walking up the parent chain (element -> ... -> processor).
     * Returns 1.0 for infinite-server processors (infinite multiplicity), so that
     * their utilization keeps the mean-busy-servers value (which may exceed 1),
     * matching the LINE convention and SolverLQNS normalization.
     *
     * @param lqn the LayeredNetworkStruct
     * @param idx the absolute element index
     * @return the host processor multiplicity, or 1.0 if none / infinite
     */
    private static double hostProcessorMult(LayeredNetworkStruct lqn, int idx) {
        int cur = idx;
        for (int hops = 0; hops <= lqn.nidx; hops++) {
            if (cur < 1 || cur > lqn.nidx) return 1.0;
            int t = (int) lqn.type.get(cur);
            if (t == LayeredNetworkElement.PROCESSOR) {
                double m = lqn.mult.get(cur);
                return (m > 0 && !Double.isInfinite(m)) ? m : 1.0;
            }
            int p = (int) lqn.parent.get(cur);
            if (p <= 0 || p == cur) return 1.0;
            cur = p;
        }
        return 1.0;
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
     * Returns the default solver options for the LDES solver.
     *
     * @return Default solver options with SolverType.LDES
     */
    public static SolverOptions defaultOptions() {
        return new LDESOptions();
    }

    // =====================================================
    // Sample Methods - State Trajectory Generation
    // =====================================================

    /**
     * Generates a sample path (state trajectory) for a specific node using LDES transient simulation.
     *
     * @param node The stateful node to sample
     * @param numEvents Number of time points to sample
     * @return SampleResult containing the state trajectory for the node
     */
    public jline.io.Ret.SampleResult sample(jline.lang.nodes.StatefulNode node, int numEvents) {
        // Configure for transient analysis; the horizon is passed explicitly via
        // timespan and options.samples (the steady-state event budget) is left
        // untouched, since the engine ignores it in transient mode.
        LDESOptions ldesOptions = (LDESOptions) this.options;
        double[] originalTimespan = ldesOptions.timespan;
        double horizon = (numEvents > 0) ? numEvents : ldesOptions.samples;
        ldesOptions.timespan = new double[]{0.0, horizon};

        try {
            // Run transient simulation
            runAnalyzer();

            // Get results
            LDESResult ldesResult = (LDESResult) this.result;
            if (ldesResult == null || ldesResult.t == null) {
                return new jline.io.Ret.SampleResult();
            }

            NetworkStruct sn = this.sn;
            int nodeIdx = node.getNodeIndex();
            int isf = (int) sn.nodeToStateful.get(0, nodeIdx);

            // Build state matrix from transient queue lengths
            int numTimePoints = ldesResult.t.getNumRows();
            int numClasses = sn.nclasses;
            Matrix t = ldesResult.t.copy();
            Matrix state = new Matrix(numTimePoints, numClasses);

            if (ldesResult.QNt != null && isf < ldesResult.QNt.length) {
                for (int k = 0; k < numClasses && k < ldesResult.QNt[isf].length; k++) {
                    Matrix classData = ldesResult.QNt[isf][k];
                    if (classData != null) {
                        for (int ti = 0; ti < numTimePoints && ti < classData.getNumRows(); ti++) {
                            state.set(ti, k, classData.get(ti, 0));
                        }
                    }
                }
            }

            return new jline.io.Ret.SampleResult("des", t, state, new Matrix(0, 0), false, nodeIdx, numTimePoints);

        } catch (Exception e) {
            line_debug(options.verbose, "LDES sample failed: " + e.getMessage());
            return new jline.io.Ret.SampleResult();
        } finally {
            // Restore original timespan
            ldesOptions.timespan = originalTimespan;
        }
    }

    /**
     * Generates an aggregated sample path for a specific node.
     *
     * For LDES, the sample() method already returns per-class queue lengths (nir),
     * so this method simply marks the result as aggregated without additional
     * marginal computation. This differs from SSA which requires toMarginal
     * processing on raw state vectors.
     *
     * @param node The stateful node to sample
     * @param numEvents Number of time points to sample
     * @return SampleResult containing the aggregated state trajectory
     */
    public jline.io.Ret.SampleResult sampleAggr(jline.lang.nodes.StatefulNode node, int numEvents) {
        jline.io.Ret.SampleResult result = sample(node, numEvents);

        if (result != null && result.state != null && result.state instanceof Matrix) {
            // LDES sample() already returns aggregated per-class queue lengths from QNt,
            // so we just need to mark the result as aggregated
            Matrix stateMatrix = (Matrix) result.state;
            return new jline.io.Ret.SampleResult("des", result.t, stateMatrix, result.event, true, result.nodeIndex, result.numEvents);
        }
        return result;
    }

    /**
     * Generates a system-wide sample path (state trajectory) for all stateful nodes.
     *
     * @param numEvents Number of time points to sample
     * @return SampleResult containing the state trajectory for all nodes
     */
    public jline.io.Ret.SampleResult sampleSys(int numEvents) {
        // Configure for transient analysis; the horizon is passed explicitly via
        // timespan and options.samples (the steady-state event budget) is left
        // untouched, since the engine ignores it in transient mode.
        LDESOptions ldesOptions = (LDESOptions) this.options;
        double[] originalTimespan = ldesOptions.timespan;
        double horizon = (numEvents > 0) ? numEvents : ldesOptions.samples;
        ldesOptions.timespan = new double[]{0.0, horizon};

        try {
            // Run transient simulation
            runAnalyzer();

            LDESResult ldesResult = (LDESResult) this.result;
            if (ldesResult == null || ldesResult.t == null) {
                return new jline.io.Ret.SampleResult();
            }

            NetworkStruct sn = this.sn;
            int numTimePoints = ldesResult.t.getNumRows();
            int numClasses = sn.nclasses;
            Matrix t = ldesResult.t.copy();

            // Build state list for each stateful node
            List<Matrix> stateList = new ArrayList<>();
            for (int isf = 0; isf < sn.nstateful; isf++) {
                Matrix nodeState = new Matrix(numTimePoints, numClasses);

                if (ldesResult.QNt != null && isf < ldesResult.QNt.length) {
                    for (int k = 0; k < numClasses && k < ldesResult.QNt[isf].length; k++) {
                        Matrix classData = ldesResult.QNt[isf][k];
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
            line_debug(options.verbose, "LDES sampleSys failed: " + e.getMessage());
            return new jline.io.Ret.SampleResult();
        } finally {
            ldesOptions.timespan = originalTimespan;
        }
    }

    /**
     * Generates an aggregated system-wide sample path for all stateful nodes.
     *
     * @param numEvents Number of time points to sample
     * @return SampleResult containing the aggregated state trajectory for all nodes
     */
    @SuppressWarnings("unchecked")
    public jline.io.Ret.SampleResult sampleSysAggr(int numEvents) {
        jline.io.Ret.SampleResult result = sampleSys(numEvents);

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

            return new jline.io.Ret.SampleResult("des", result.t, aggregatedStateList, result.event, true, result.numEvents);
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
     * @param numEvents Number of time points to sample
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleResult containing the state trajectory for the node
     */
    public jline.io.Ret.SampleResult stream(jline.lang.nodes.StatefulNode node, int numEvents, StreamingOptions streamingOptions) {
        try {
            // Initialize streaming collector
            this.sn = this.model.getStruct(true);
            this.stream = new Collector(streamingOptions, this.sn);

            // Run sample with streaming active
            jline.io.Ret.SampleResult result = sample(node, numEvents);

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
     * @param numEvents Number of time points to sample
     * @param streamingOptions Configuration for streaming (endpoint, mode, frequency)
     * @return SampleResult containing the aggregated state trajectory
     */
    public jline.io.Ret.SampleResult streamAggr(jline.lang.nodes.StatefulNode node, int numEvents, StreamingOptions streamingOptions) {
        try {
            // Initialize streaming collector
            this.sn = this.model.getStruct(true);
            this.stream = new Collector(streamingOptions, this.sn);

            // Run sampleAggr with streaming active
            jline.io.Ret.SampleResult result = sampleAggr(node, numEvents);

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
     * Estimates the steady-state probability of a specific state at a node via LDES simulation.
     *
     * This method runs a LDES simulation and estimates state probabilities by computing
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
     * Estimates the joint steady-state probability of the entire system state via LDES simulation.
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
    // Markov reward methods
    // =====================================================

    /**
     * Steady-state expected reward for every reward function defined on the model.
     *
     * The reward integral E[r] = (1/T) integral_0^T r(X(t)) dt is accumulated inside
     * the discrete-event simulation on the exact integer joint state, so it is correct
     * for nonlinear rewards (e.g. E[n^2]) and converges to the stationary expectation
     * sum_s pi(s) r(s) computed exactly by the CTMC solver.
     *
     * @return map from reward name to expected reward value in steady state
     * @throws IllegalStateException if no rewards are defined on the model
     */
    public Map<String, Double> getAvgReward() {
        // Refresh the struct so the reward functions are visible to the simulation.
        this.sn = this.model.getStruct(true);
        NetworkStruct sn = this.sn;
        if (sn.reward == null || sn.reward.isEmpty()) {
            throw new IllegalStateException(
                "No rewards defined. Use model.setReward(name, rewardFn) before calling reward analysis.");
        }
        try {
            runAnalyzer();
        } catch (Exception e) {
            throw new RuntimeException("LDES reward analysis failed: " + e.getMessage(), e);
        }
        if (this.result instanceof LDESResult) {
            LDESResult r = (LDESResult) this.result;
            if (r.avgReward != null) {
                return new LinkedHashMap<String, Double>(r.avgReward);
            }
        }
        throw new RuntimeException("LDES reward computation produced no result.");
    }

    /**
     * Steady-state expected reward for a single named reward function.
     *
     * @param rewardName the reward name
     * @return the expected reward value in steady state
     * @throws IllegalArgumentException if the reward name is not defined
     */
    public double getAvgReward(String rewardName) {
        Map<String, Double> all = getAvgReward();
        Double value = all.get(rewardName);
        if (value == null) {
            throw new IllegalArgumentException("Reward '" + rewardName + "' not found.");
        }
        return value;
    }

    /**
     * Single-run transient reward trajectory r(X(t)) recorded along the simulated path.
     *
     * The associated time vector is available via {@link #getRewardTimeVector()}.
     *
     * @param rewardName optional reward name to filter; if null, all rewards are returned
     * @return map from reward name to the reward time series
     * @throws IllegalStateException if no rewards are defined on the model
     */
    public Map<String, double[]> getTranReward(String rewardName) {
        this.sn = this.model.getStruct(true);
        NetworkStruct sn = this.sn;
        if (sn.reward == null || sn.reward.isEmpty()) {
            throw new IllegalStateException(
                "No rewards defined. Use model.setReward(name, rewardFn) before calling getTranReward.");
        }
        try {
            runAnalyzer();
        } catch (Exception e) {
            throw new RuntimeException("LDES transient reward analysis failed: " + e.getMessage(), e);
        }
        Map<String, double[]> result = new LinkedHashMap<String, double[]>();
        if (this.result instanceof LDESResult) {
            LDESResult r = (LDESResult) this.result;
            if (r.tranReward != null) {
                for (Map.Entry<String, double[]> e : r.tranReward.entrySet()) {
                    if (rewardName == null || rewardName.equals(e.getKey())) {
                        result.put(e.getKey(), e.getValue());
                    }
                }
                return result;
            }
        }
        throw new RuntimeException("LDES transient reward computation produced no result.");
    }

    /**
     * Transient reward trajectory for all reward functions.
     *
     * @return map from reward name to the reward time series
     */
    public Map<String, double[]> getTranReward() {
        return getTranReward(null);
    }

    /**
     * Time vector associated with {@link #getTranReward()}.
     *
     * @return the time points, or null if not available
     */
    public Matrix getRewardTimeVector() {
        if (this.result instanceof LDESResult) {
            return ((LDESResult) this.result).rewardTime;
        }
        return null;
    }

    /**
     * Names of the reward functions defined on the model.
     *
     * @return list of reward names (empty if none are defined)
     */
    public List<String> getRewardNames() {
        NetworkStruct sn = this.sn != null ? this.sn : this.model.getStruct(true);
        if (sn.reward == null) {
            return new ArrayList<String>();
        }
        return new ArrayList<String>(sn.reward.keySet());
    }

    /**
     * Runs the simulation with state-histogram export enabled and returns the exact
     * joint-state residence-time histogram. Row {@code i} of the returned space matrix
     * is an aggregated system state (CTMC {@code stateSpaceAggr} layout, station-major)
     * and entry {@code i} of the time vector is its residence time. Host languages that
     * cannot pass a Java reward function into the simulation use this to evaluate
     * arbitrary Markov rewards: E[r] = sum_i (time_i / sum(time)) * r(state_i).
     *
     * @return a two-element array {space, time}; both null if the histogram is empty
     */
    public Matrix[] getStateHistogram() {
        if (this.options instanceof LDESOptions) {
            ((LDESOptions) this.options).exportStateHistogram = true;
        }
        this.sn = this.model.getStruct(true);
        try {
            runAnalyzer();
        } catch (Exception e) {
            throw new RuntimeException("LDES state-histogram analysis failed: " + e.getMessage(), e);
        }
        if (this.result instanceof LDESResult) {
            LDESResult r = (LDESResult) this.result;
            return new Matrix[] { r.stateHistogramSpace, r.stateHistogramTime };
        }
        return new Matrix[] { null, null };
    }

    /**
     * Runs the simulation with state export enabled and returns the integer joint-state
     * trajectory along the sampled path: {@code {space, time}} where row {@code i} of
     * space is the aggregated state (CTMC {@code stateSpaceAggr} layout) at {@code time[i]}.
     * Host languages evaluate their reward functions on each state to obtain the
     * transient reward trajectory r(X(t)).
     *
     * @return a two-element array {space, time}; both null if unavailable
     */
    public Matrix[] getStateTrajectory() {
        if (this.options instanceof LDESOptions) {
            ((LDESOptions) this.options).exportStateHistogram = true;
        }
        this.sn = this.model.getStruct(true);
        try {
            runAnalyzer();
        } catch (Exception e) {
            throw new RuntimeException("LDES state-trajectory analysis failed: " + e.getMessage(), e);
        }
        if (this.result instanceof LDESResult) {
            LDESResult r = (LDESResult) this.result;
            return new Matrix[] { r.stateTrajectorySpace, r.stateTrajectoryTime };
        }
        return new Matrix[] { null, null };
    }

    // =====================================================
    // Transient CDF Methods
    // =====================================================

    /**
     * Returns cumulative distribution functions of response times during transient analysis.
     * Uses response time samples collected during LDES simulation to compute empirical CDFs.
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
        if (this.result == null || !(this.result instanceof LDESResult)) {
            try {
                runAnalyzer();
            } catch (Exception e) {
                throw new RuntimeException("Failed to run LDES simulation for transient CDF: " + e.getMessage(), e);
            }
        }

        LDESResult ldesResult = (LDESResult) this.result;
        NetworkStruct sn = this.sn != null ? this.sn : this.model.getStruct(true);

        DistributionResult distResult = new DistributionResult(sn.nstations, sn.nclasses, "response_time");

        if (ldesResult.respTimeSamples == null) {
            return distResult;
        }

        // Compute empirical CDF for each station-class pair
        for (int i = 0; i < sn.nstations; i++) {
            for (int r = 0; r < sn.nclasses; r++) {
                if (i < ldesResult.respTimeSamples.length && r < ldesResult.respTimeSamples[i].length) {
                    List<Double> samples = ldesResult.respTimeSamples[i][r];
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
     * For LDES, passage times are equivalent to response times in single-visit networks.
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
     * Computes transient state probabilities at a specific node over time using LDES simulation.
     *
     * @param node the stateful node to analyze
     * @return ProbabilityResult containing transient probability data
     */
    public ProbabilityResult getTranProb(StatefulNode node) {
        // Ensure finite timespan
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProb in SolverLDES requires a finite timespan. " +
                    "Use: SolverLDES(model, \"timespan\", new double[]{0, T})");
        }

        // Run transient simulation
        try {
            getTranAvg();
        } catch (Exception e) {
            throw new RuntimeException("Failed to run transient LDES simulation: " + e.getMessage(), e);
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
            throw new RuntimeException("getTranProbAggr in SolverLDES requires a finite timespan. " +
                    "Use: SolverLDES(model, \"timespan\", new double[]{0, T})");
        }

        // Run transient simulation
        try {
            getTranAvg();
        } catch (Exception e) {
            throw new RuntimeException("Failed to run transient LDES simulation: " + e.getMessage(), e);
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
     * Computes transient system-wide state probabilities over time using LDES simulation.
     *
     * @return ProbabilityResult containing transient system probability data
     */
    public ProbabilityResult getTranProbSys() {
        // Ensure finite timespan
        if (this.options.timespan == null || !Double.isFinite(this.options.timespan[1])) {
            throw new RuntimeException("getTranProbSys in SolverLDES requires a finite timespan. " +
                    "Use: SolverLDES(model, \"timespan\", new double[]{0, T})");
        }

        // Run transient simulation
        try {
            getTranAvg();
        } catch (Exception e) {
            throw new RuntimeException("Failed to run transient LDES simulation: " + e.getMessage(), e);
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
            throw new RuntimeException("getTranProbSysAggr in SolverLDES requires a finite timespan. " +
                    "Use: SolverLDES(model, \"timespan\", new double[]{0, T})");
        }

        // Run transient simulation
        try {
            getTranAvg();
        } catch (Exception e) {
            throw new RuntimeException("Failed to run transient LDES simulation: " + e.getMessage(), e);
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
