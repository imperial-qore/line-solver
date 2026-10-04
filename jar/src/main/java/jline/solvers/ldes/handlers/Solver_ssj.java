/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * @file Solver_ssj.java
 * @brief SSJ-based discrete event simulation engine for queueing networks and Petri nets.
 *
 * Java port of Solver_ssj.kt. See the original Kotlin file for full documentation.
 */
package jline.solvers.ldes.handlers;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.mam.Map_sample;
import jline.api.mam.Me_sample;
import jline.api.mam.Rap_sample;
import jline.api.mam.Mmap_sample;
import jline.api.mam.Dmap_sample;
import static jline.io.SysUtils.lineTempName;
import static jline.io.InputOutput.line_warning;
import jline.lang.NodeParam;
import java.io.FileWriter;
import jline.lang.NetworkStruct;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.HeteroSchedPolicy;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.NodeType;
import jline.lang.constant.PollingType;
import jline.lang.constant.RemovalPolicy;
import jline.lang.constant.ProcessType;
import jline.lang.constant.ReplacementStrategy;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.constant.TimingStrategy;
import jline.lang.JobClass;
import jline.lang.nodes.Cache;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodeparam.QueueNodeParam;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.processes.DiscreteDistribution;
import jline.lang.processes.Replayer;
import jline.solvers.SolverOptions;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.LDESResult;
import jline.streaming.Collector;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.simevents.Event;
import umontreal.ssj.simevents.Simulator;
import umontreal.ssj.stat.Tally;

import java.io.File;
import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Set;

/**
 * Top-level Java translation of Solver_ssj.kt.
 *
 * Hosts the public top-level functions (solver_ssj, solver_ssj_transient) as
 * static methods and the package-internal {@code SSJSimulator} class as a
 * nested static class.
 */
public final class Solver_ssj {

    private Solver_ssj() {
        // Static utility container
    }

    /**
     * Custom random variate generator for deterministic (constant) values.
     * Used when simulating deterministic distributions in SSJ.
     */
    public static class ConstantGen extends RandomVariateGen {
        private final double value;

        public ConstantGen(RandomStream stream, double value) {
            super(stream, null);
            this.value = value;
        }

        @Override
        public double nextDouble() {
            return value;
        }
    }

    /**
     * Random variate generator for a geometric distribution supported on
     * {1,2,...}, i.e. the trial index of the first success.
     *
     * <p>SSJ's {@code GeometricGen} counts failures and is therefore supported on
     * {0,1,...}, which is one slot short of {@link jline.lang.processes.Geometric}
     * and admits zero-length intervals. Shifting by one restores both the support
     * and the moments: mean 1/p and SCV 1-p.
     */
    public static class ShiftedGeometricGen extends RandomVariateGen {
        private final double p;

        public ShiftedGeometricGen(RandomStream stream, double p) {
            super(stream, null);
            if (!(p > 0.0) || p > 1.0) {
                throw new IllegalArgumentException("Geometric success probability p=" + p
                        + " must lie in (0,1]");
            }
            this.p = p;
        }

        @Override
        public double nextDouble() {
            if (p >= 1.0) {
                // Degenerate case: the first trial always succeeds. Drawing from
                // the stream anyway keeps the substream consumption independent
                // of p, so a parameter sweep stays synchronized.
                this.stream.nextDouble();
                return 1.0;
            }
            double u = this.stream.nextDouble();
            return FastMath.ceil(FastMath.log(1.0 - u) / FastMath.log(1.0 - p));
        }
    }

    /**
     * Steady-state queueing network simulation using the SSJ library.
     */
    public static LDESResult solver_ssj(NetworkStruct sn, SolverOptions options, Collector stream) {
        assertNoNonCyclicScheduleInSteadyState(sn);
        Matrix initSol = options.init_sol;
        SSJSimulator simulator = new SSJSimulator(sn, options, initSol, stream);
        double maxEvents = (double) options.samples;
        simulator.simulate(maxEvents);
        return simulator.getLDESResult();
    }

    /**
     * Rejects a steady-state run of a model carrying a non-cyclic rate schedule.
     *
     * <p>A non-cyclic NHPP has zero intensity past its horizon, so its total
     * intensity is finite: the process emits a bounded number of events and then
     * falls silent forever. The steady state is therefore the empty system, and
     * a sample-count-driven run cannot reach its target -- it averages over
     * whatever horizon it happened to stop at and reports a plausible-looking
     * but meaningless table (throughput exceeding the arrival rate, for
     * instance). A non-cyclic NHPP is a transient construct; failing here is
     * deliberate, the alternative being a silently wrong result.
     */
    private static void assertNoNonCyclicScheduleInSteadyState(NetworkStruct sn) {
        if (sn == null || sn.procid == null) {
            return;
        }
        for (Map.Entry<jline.lang.nodes.Station, Map<jline.lang.JobClass, ProcessType>> entry
                : sn.procid.entrySet()) {
            jline.lang.nodes.Station station = entry.getKey();
            for (Map.Entry<jline.lang.JobClass, ProcessType> inner : entry.getValue().entrySet()) {
                if (inner.getValue() != ProcessType.NHPP) {
                    continue;
                }
                Map<jline.lang.JobClass, MatrixCell> procRow =
                        (sn.proc != null) ? sn.proc.get(station) : null;
                MatrixCell proc = (procRow != null) ? procRow.get(inner.getKey()) : null;
                if (proc != null && proc.get(2) != null && proc.get(2).get(0, 0) == 0.0) {
                    throw new RuntimeException(
                            "LDES: station '" + station.getName() + "' class '"
                            + inner.getKey().getName() + "' carries a non-cyclic NHPP,"
                            + " whose intensity is zero past its horizon, so the steady"
                            + " state is the empty system and a steady-state run is not"
                            + " meaningful. Use transient analysis (getTranAvg with"
                            + " options.timespan inside the schedule horizon), or make"
                            + " the schedule cyclic.");
                }
            }
        }
    }

    public static LDESResult solver_ssj(NetworkStruct sn, SolverOptions options) {
        return solver_ssj(sn, options, null);
    }

    /**
     * Transient analysis of queueing network using the SSJ library.
     */
    public static LDESResult solver_ssj_transient(NetworkStruct sn, SolverOptions options, Collector stream) {
        Matrix initSol = options.init_sol;
        SSJSimulator simulator = new SSJSimulator(sn, options, initSol, stream);
        double[] timespan = options.timespan;
        if (timespan != null && timespan.length >= 2) {
            simulator.simulateTransient(timespan[1]);
        }
        return simulator.getTransientLDESResult();
    }

    /**
     * Average the transient time series (QNt/UNt/TNt value columns) across
     * independent replications, preserving the (shared) time column. Yields the
     * Monte-Carlo estimate of the transient mean E[.](t): a single sample path
     * is a random realization, not the transient mean, so the LN transient
     * oracle averages over replications. The bucket grid is deterministic
     * (timespan[1]/N) so bucket k has the same time across runs. Steady-state
     * fields are taken from the first run. Used by the parallel replication
     * analyzer; public so it can drive the ensemble transient there.
     */
    public static LDESResult averageTransientResults(List<LDESResult> runs) {
        LDESResult base = runs.get(0);
        int R = runs.size();
        if (base.QNt == null) {
            return base;
        }
        int nS = base.QNt.length;
        for (int i = 0; i < nS; i++) {
            int nK = base.QNt[i].length;
            for (int k = 0; k < nK; k++) {
                averageSeriesInto(base.QNt[i][k], runs, R, i, k, 0);
                averageSeriesInto(base.UNt[i][k], runs, R, i, k, 1);
                averageSeriesInto(base.TNt[i][k], runs, R, i, k, 2);
            }
        }
        return base;
    }

    // which: 0=QNt, 1=UNt, 2=TNt. Accumulates value column (col 0) across runs
    // in place on `dst` (the first run's matrix) then divides by R.
    private static void averageSeriesInto(Matrix dst, List<LDESResult> runs, int R, int i, int k, int which) {
        if (dst == null || dst.getNumRows() == 0) {
            return;
        }
        int rows = dst.getNumRows();
        for (int r = 1; r < R; r++) {
            LDESResult run = runs.get(r);
            Matrix src = (which == 0) ? run.QNt[i][k] : (which == 1) ? run.UNt[i][k] : run.TNt[i][k];
            if (src == null) {
                continue;
            }
            int rr = Math.min(rows, src.getNumRows());
            for (int t = 0; t < rr; t++) {
                dst.set(t, 0, dst.get(t, 0) + src.get(t, 0));
            }
        }
        for (int t = 0; t < rows; t++) {
            dst.set(t, 0, dst.get(t, 0) / R);
        }
    }

    public static LDESResult solver_ssj_transient(NetworkStruct sn, SolverOptions options) {
        return solver_ssj_transient(sn, options, null);
    }

    // =========================================================================
    // SSJSimulator - Java port of the Kotlin internal class.
    // =========================================================================

    /**
     * Core discrete event simulation engine using SSJ library.
     *
     * <p>Java port of the Kotlin {@code internal class SSJSimulator}.
     */
    static final class SSJSimulator {

        // Event list and clock for this run. An instance is used rather than the
        // static Sim facade, whose default Simulator is shared JVM-wide: parallel
        // replications would otherwise interleave on one event list and one clock.
        // Every event is bound to this instance through SimEvent below.
        private final Simulator ssjSim = new Simulator();

        /**
         * Base class of every event scheduled by this simulator, binding the event
         * to the enclosing run's {@link #ssjSim} rather than to SSJ's JVM-wide
         * default simulator. All event classes below must extend this, not Event.
         */
        private abstract class SimEvent extends Event {
            SimEvent() {
                super(SSJSimulator.this.ssjSim);
                if (SSJSimulator.this.slotted) {
                    setPriority(slotPhase());
                }
            }

            /**
             * Position of this event within a slot boundary, used only in slotted
             * mode. SSJ orders the event list by (time, priority) ascending, so a
             * lower phase fires first. Without this, events falling on the same
             * integer instant would be ordered by insertion, which is an emergent
             * rather than a modelled convention.
             *
             * <p>The default is the completion phase because most events in the
             * engine are service completions; the arrival, routing and bookkeeping
             * classes override it.
             */
            protected double slotPhase() {
                return SLOT_PHASE_COMPLETION;
            }
        }

        /** Service completions and other departures resolve first within a slot. */
        static final double SLOT_PHASE_COMPLETION = 1.0;
        /** Internal movements triggered by a completion resolve next. */
        static final double SLOT_PHASE_INTERNAL = 2.0;
        /** External arrivals see the post-completion state of the slot. */
        static final double SLOT_PHASE_ARRIVAL = 3.0;
        /** Sampling and control events observe the settled state of the slot. */
        static final double SLOT_PHASE_BOOKKEEPING = 9.0;

        // Network parameters derived from sn
        private final NetworkStruct sn;
        private final SolverOptions options;
        private final Matrix initSol;
        private final Collector stream;

        private final int numClasses;
        private final int numStations;
        private final int numNodes;
        private final long seed;

        /** Discrete time scale: every interval must land on the slot lattice. */
        private final boolean slotted;
        /** Slot length in model time units; meaningful only when slotted. */
        private final double slotLength;

        /**
         * Server state for setup and delayoff support.
         */
        private enum ServerState {
            OFF,
            SETUP,
            ACTIVE,
            DELAYOFF
        }

        // Node classification
        private final List<Integer> sourceNodes = new ArrayList<Integer>();
        private final List<Integer> sourceStations = new ArrayList<Integer>();
        private final List<Integer> serviceNodes = new ArrayList<Integer>();
        private final List<Integer> serviceStations = new ArrayList<Integer>();
        private final List<Boolean> isDelayNode = new ArrayList<Boolean>();
        private final List<Integer> sinkNodes = new ArrayList<Integer>();
        private final List<Integer> loggerNodes = new ArrayList<Integer>();
        private final List<Integer> routerNodes = new ArrayList<Integer>();
        private final List<Integer> classSwitchNodes = new ArrayList<Integer>();
        private final List<Integer> forkNodes = new ArrayList<Integer>();
        private final List<Integer> joinNodes = new ArrayList<Integer>();
        private final List<Integer> joinStations = new ArrayList<Integer>();
        private final List<Integer> placeNodes = new ArrayList<Integer>();
        private final List<Integer> transitionNodes = new ArrayList<Integer>();
        private final List<Integer> cacheNodes = new ArrayList<Integer>();

        /** A request parked at a cache as a delayed hit, awaiting fetch completion. */
        private static final class HeldRequest {
            final int accessClass;       // class with which the request accessed the cache
            final double holdTime;       // simulation time at which it was parked
            final long jobId;
            HeldRequest(int accessClass, double holdTime, long jobId) {
                this.accessClass = accessClass;
                this.holdTime = holdTime;
                this.jobId = jobId;
            }
        }

        /** Cache state information */
        private static final class CacheStateInfo {
            final int nodeIdx;
            final int numItems;
            final int[] levelCapacities;
            final ReplacementStrategy replacementStrategy;
            double qlru = 1.0;                // q-LRU admission probability on a miss
            final LinkedList<Integer>[] levels;
            DiscreteAccessSampler[] accessSamplers; // [class] -> item sampler; null entry = class has no popularity (bound at runtime)
            Matrix[][] accost;                // [class][item] -> (h+1)x(h+1) list-access probabilities (row 0 = miss row)
            final int[] hitClass;
            final int[] missClass;
            long[] totalHits;
            long[] totalMisses;
            long[][] hitsPerList;             // [accessClass][list] -> hit count in that list
            double[][] itemLevelTime;         // [item][list] -> time-integrated presence in that list
            double lastContentUpdateTime;     // sim time of the last itemLevelTime accumulation
            double occupancyStartTime;        // sim time occupancy/count statistics were (re)started
            // --- retrieval system (delayed hits) ---
            boolean hasRetrieval;
            int[][] retrievalClass;           // [item][accessClass] -> retrieval class index, or -1
            int[] retrievalClassToItem;       // [class] -> item if class is a retrieval class, else -1
            boolean[] inFlight;               // [item] -> a fetch for this item is in progress
            List<HeldRequest>[] heldDelayedHits;  // [item] -> requests parked while item is in flight
            long[] totalDelayedHits;          // [accessClass] -> delayed-hit count
            double totalDelayedHitWait;       // accumulated waiting time of released delayed hits
            double[] fetchStartTime;          // [item] -> sim time the in-flight fetch began
            double totalFetchTime;            // accumulated fetch sojourn of completed fetches
            long completedFetches;            // number of completed fetches

            CacheStateInfo(int nodeIdx, int numItems, int[] levelCapacities,
                           ReplacementStrategy replacementStrategy,
                           LinkedList<Integer>[] levels,
                           int[] hitClass, int[] missClass,
                           long[] totalHits, long[] totalMisses) {
                this.nodeIdx = nodeIdx;
                this.numItems = numItems;
                this.levelCapacities = levelCapacities;
                this.replacementStrategy = replacementStrategy;
                this.levels = levels;
                this.accessSamplers = null;   // deferred to initializeCacheSamplers
                this.hitClass = hitClass;
                this.missClass = missClass;
                this.totalHits = totalHits;
                this.totalMisses = totalMisses;
                this.hitsPerList = new long[totalHits.length][levels.length];
                this.itemLevelTime = new double[numItems][levels.length];
                this.lastContentUpdateTime = 0.0;
                this.occupancyStartTime = 0.0;
                this.hasRetrieval = false;
            }
        }
        private CacheStateInfo[] cacheStates;

        /** Sentinel returned by processCacheAccess when a request is parked as a delayed hit. */
        private static final int CACHE_HELD = Integer.MIN_VALUE;

        // Place token storage: [placeListIdx] -> token counts per class.
        // For an ordinary place this is the marking and also the tokens available to
        // output transitions. For a queueing place (QPN embedded queue) this is the total
        // marking (queued + in-service + depository); the tokens actually available to
        // output transitions are held separately in placeDepository.
        private int[][] placeTokens;

        // --- Queueing place (QPN embedded queue) state, indexed by placeListIdx ---
        /** True if the place has an embedded queue (service assigned via Place.setService). */
        private boolean[] isQueueingPlace;
        /** Number of servers of the embedded queue (Integer.MAX_VALUE for INF). */
        private int[] placeNumServers;
        /** Tokens in the depository, available to output transitions. */
        private int[][] placeDepository;
        /** Number of tokens currently in service at the embedded queue. */
        private int[] placeBusy;
        /** Number of servers currently busy per class (for per-class utilization). */
        private int[][] currentPlaceInService;
        /** Time integral of busy servers per class (for utilization). */
        private double[][] totalPlaceBusyTime;
        /** Last time the busy-server integral was flushed, per place. */
        private double[] lastPlaceBusyUpdateTime;
        /** FIFO waiting line of token classes not yet in service. */
        private java.util.ArrayDeque<Integer>[] placeWaiting;
        /** Embedded-queue service process type per class. */
        private ProcessType[][] placeSvcType;
        /** Embedded-queue service generators (EXP/DET) per class. */
        private RandomVariateGen[][] placeSvcGen;
        /** Embedded-queue service process matrices (renewal PH families) per class. */
        private MatrixCell[][] placeSvcProc;
        /** RNG for renewal phase-type embedded-queue service sampling per class. */
        private java.util.Random[][] placeSvcRng;
        /** Inverse-CDF samplers for ME embedded-queue service processes per class. */
        private Me_sample.MeSampler[][] placeSvcMeSampler;

        /** Transition mode information */
        private static final class TransitionModeInfo {
            final int modeIdx;
            final String modeName;
            final TimingStrategy timingStrategy;
            final int priority;
            final double weight;
            final int numServers;
            final int[][] enablingConditions;
            final int[][] inhibitingConditions;
            final int[][] firingOutcomes;
            // Marking-dependent firing-rate multiplier g(marking); null == unit.
            final SerializableFunction<Matrix, Double> firingDep;

            TransitionModeInfo(int modeIdx, String modeName, TimingStrategy timingStrategy,
                               int priority, double weight, int numServers,
                               int[][] enablingConditions, int[][] inhibitingConditions,
                               int[][] firingOutcomes,
                               SerializableFunction<Matrix, Double> firingDep) {
                this.modeIdx = modeIdx;
                this.modeName = modeName;
                this.timingStrategy = timingStrategy;
                this.priority = priority;
                this.weight = weight;
                this.numServers = numServers;
                this.enablingConditions = enablingConditions;
                this.inhibitingConditions = inhibitingConditions;
                this.firingOutcomes = firingOutcomes;
                this.firingDep = firingDep;
            }
        }

        // Transition parameters: [transitionNodeIdx] -> list of modes
        private List<TransitionModeInfo>[] transitionModes;

        // Transition firing distributions: [transitionNodeIdx][modeIdx] -> distribution generator
        private RandomVariateGen[][] transitionFiringGens;

        // Tokens currently being processed by transitions (for multi-server)
        private int[][] transitionInService;

        // In-flight firing clocks of marking-dependent (dependent) modes. Their
        // rate changes with the marking, so on every marking change they are
        // canceled and resampled at the new rate (exact for exponential firing).
        private final List<TransitionFiring> inflightDependent = new ArrayList<TransitionFiring>();

        // Tokens in transit at transitions: [transListIdx][classIdx] -> token count
        private int[][] tokensInTransit;
        // Time-weighted tokens in transit
        private double[][] totalTransitTokenTime;
        private double[] lastTransitUpdateTime;

        // Place statistics tracking
        private double[][] totalPlaceTokenTime;
        private int[][] placeCompletions;
        private double[] lastPlaceUpdateTime;
        private double[][] placeTransitTokenTime;
        private int[][] placeTokensInTransit;
        private double[] lastPlaceTransitUpdateTime;

        /** Fork tracking record (parent -> outstanding tasks). */
        private static final class ForkJobInfo {
            final long parentJobId;
            final int parentClassId;
            final double parentSystemArrivalTime;
            final int forkNodeIdx;
            final int totalTasks;
            int completedTasks;
            double firstJoinArrivalTime;
            final List<Double> forkedJobJoinArrivalTimes;
            final List<Integer> forkedJobJoinClasses;

            ForkJobInfo(long parentJobId, int parentClassId, double parentSystemArrivalTime,
                        int forkNodeIdx, int totalTasks) {
                this.parentJobId = parentJobId;
                this.parentClassId = parentClassId;
                this.parentSystemArrivalTime = parentSystemArrivalTime;
                this.forkNodeIdx = forkNodeIdx;
                this.totalTasks = totalTasks;
                this.completedTasks = 0;
                this.firstJoinArrivalTime = -1.0;
                this.forkedJobJoinArrivalTimes = new ArrayList<Double>();
                this.forkedJobJoinClasses = new ArrayList<Integer>();
            }
        }

        /** A single forked sibling job. */
        private static final class ForkedJob {
            final long forkJobId;
            final long parentJobId;
            final int classId;
            final int priority;
            final double systemArrivalTime;
            final double queueArrivalTime;
            final double randomRank;

            ForkedJob(long forkJobId, long parentJobId, int classId, int priority,
                      double systemArrivalTime, double queueArrivalTime, double randomRank) {
                this.forkJobId = forkJobId;
                this.parentJobId = parentJobId;
                this.classId = classId;
                this.priority = priority;
                this.systemArrivalTime = systemArrivalTime;
                this.queueArrivalTime = queueArrivalTime;
                this.randomRank = randomRank;
            }
        }

        // Maps parent job ID to fork info (for join synchronization)
        private final Map<Long, ForkJobInfo> forkJobInfoMap = new HashMap<Long, ForkJobInfo>();

        // Maps forked job ID to parent job ID
        private final Map<Long, Long> forkedJobParentMap = new HashMap<Long, Long>();

        // Fork parameters
        private int[] forkFanOut;

        // Join parameters
        private int[] joinToForkMap;
        private int[] forkToJoinMap;

        // Join strategies
        private JoinStrategy[][] joinStrategies;

        // Join required counts
        private int[][] joinRequired;

        // Join statistics tracking
        private double[][] totalJoinQueueTime;
        private int[][] joinCompletions;
        private double[][] lastJoinUpdateTime;
        private int[][] currentJoinQueueLength;
        private Tally[][] joinResponseTimeTally;
        private int[][] arrivedAtJoin;
        // Siblings discarded at a quorum/PARTIAL Join: a task arriving after the
        // quorum already fired finds its parent record removed and is dropped.
        private int[][] droppedByJoin;

        // Counter for generating unique forked job IDs
        private long nextForkedJobId = 0L;

        // Class switch matrix storage
        private double[][][] classSwitchMatrices; // [nodeIdx][fromClass][toClass]

        // Routing strategy tracking
        private RoutingStrategy[][] nodeRoutingStrategies;
        private int[][] roundRobinCounters;
        private double[][][] wrrobinWeights;
        private int[][][] rroutlinks;
        private int sqD = 2;
        private int[][] sqDByNodeClass;

        // Derived arrays
        private double[][] lambdas;
        private double[][] mus;
        private int[] numServers;
        private int[] bufferCapacities;
        private int[][] classCapacities;

        // Finite capacity regions
        private final int numRegions;
        private final List<Integer> fcRegionIndices;
        private final List<Integer> fcRegionGlobalMax;
        private final double[] fcRegionGlobalMaxMem;  // per-region memory budget (-1 = unbounded)
        private final Matrix fcRegionClassMax;
        private final boolean[][] fcRegionDropRule;
        private final double[][][] fcRegionLinConA;
        private final double[][] fcRegionLinConb;

        // Blocking policy support (BAS and BBS)
        private int[][] stationDropRule;

        /** BBS blocked-server record. */
        private static final class BBSBlockedServer {
            final Customer customer;
            final int destQueueIdx;
            final int destClassId;
            final int sourceQueueIdx;
            final int serverId;
            final int sourceClassId;
            final double blockStartTime;

            BBSBlockedServer(Customer customer, int destQueueIdx, int destClassId,
                             int sourceQueueIdx, int serverId, int sourceClassId,
                             double blockStartTime) {
                this.customer = customer;
                this.destQueueIdx = destQueueIdx;
                this.destClassId = destClassId;
                this.sourceQueueIdx = sourceQueueIdx;
                this.serverId = serverId;
                this.sourceClassId = sourceClassId;
                this.blockStartTime = blockStartTime;
            }
        }

        private final Map<Integer, List<BBSBlockedServer>> bbsBlockedServers = new HashMap<Integer, List<BBSBlockedServer>>();
        private final Map<Integer, Set<Integer>> bbsDestinationToSources = new HashMap<Integer, Set<Integer>>();

        /** BAS waiting-job record. */
        private static final class BASWaitingJob {
            final Customer customer;
            final int destQueueIdx;
            final int destClassId;
            final int sourceQueueIdx;
            final int serverId;
            final int sourceClassId;
            final double arrivalTime;

            BASWaitingJob(Customer customer, int destQueueIdx, int destClassId,
                          int sourceQueueIdx, int serverId, int sourceClassId,
                          double arrivalTime) {
                this.customer = customer;
                this.destQueueIdx = destQueueIdx;
                this.destClassId = destClassId;
                this.sourceQueueIdx = sourceQueueIdx;
                this.serverId = serverId;
                this.sourceClassId = sourceClassId;
                this.arrivalTime = arrivalTime;
            }
        }

        private final Map<Integer, List<BASWaitingJob>> basOutgoingBuffer = new HashMap<Integer, List<BASWaitingJob>>();
        private final Map<Integer, Set<Integer>> basDestinationToSources = new HashMap<Integer, Set<Integer>>();

        // Counts
        private int numServiceNodes = 0;
        private int numSources = 0;

        // Tracking for max queue length (used internally)
        private int maxQueueLengthReached = 0;

        // Class priorities (from sn.classprio)
        private final int[] classPrio;

        // Class deadlines (from sn.classdeadline)
        private final double[] classDeadline;

        // Scheduling strategies per service node
        private SchedStrategy[] schedStrategies;

        // Polling scheduling state
        private boolean[] isPollingStation;
        private PollingType[] pollingType;
        private int[] pollingK;
        private int[] pollingCurrentClass;
        private int[] pollingJobsServedInRound;
        private int[] pollingGateSize;
        private boolean[] pollingInSwitchover;
        private LinkedList<Customer>[][] pollingQueues;
        private RandomVariateGen[][] pollingSwitchoverGens;

        // Class type classification
        private boolean[] isOpenClass;
        private boolean[] isClosedClass;
        private int[] closedClassPopulation;
        private int[] referenceStation;

        // Signal class detection
        private boolean[] isSignalClass;
        private boolean[] isNegativeSignal;
        private boolean hasNegativeSignals = false;
        private DiscreteDistribution[] signalRemovalDist;
        private RemovalPolicy[] signalRemovalPolicy;
        private boolean[] isCatastropheSignal;
        private boolean hasCatastropheSignals = false;
        // True when any signal class removes jobs (negative or catastrophe).
        // In-service and delay job tracking must be active for both: a model
        // whose only signal is a catastrophe has hasNegativeSignals == false,
        // so gating the tracking on that flag alone left in-service jobs
        // invisible to the removal routine and under-removed at every station.
        private boolean hasRemovalSignals = false;

        // REPLY signal detection
        private boolean[] isReplySignal;
        private boolean hasReplySignals = false;

        private int[] synchCallReplyClass;

        // Spawn-on-completion mapping (sn.classspawn): completing a job of
        // class k injects a fresh job of class spawnClassOf[k] at the same
        // station (LQN phase-2 continuations), -1 if none
        private int[] spawnClassOf;

        // Lazy per-class cache of whether a spawn target class routes into a
        // Join node: 0 unknown, 1 no, 2 yes. Such a continuation stands in
        // for a forked sibling at the Join, so it inherits the fork identity
        // of the job that triggered the spawn.
        private byte[] spawnJoinReach;

        // Server blocking state
        private boolean[][] serverBlocked;

        /** Pending reply tracking for synchronous calls. */
        private static final class PendingReply {
            final long jobId;
            final int originalClassId;
            final int queueIdx;
            final int serverId;
            final double blockStartTime;

            PendingReply(long jobId, int originalClassId, int queueIdx,
                         int serverId, double blockStartTime) {
                this.jobId = jobId;
                this.originalClassId = originalClassId;
                this.queueIdx = queueIdx;
                this.serverId = serverId;
                this.blockStartTime = blockStartTime;
            }
        }
        private final Map<Long, PendingReply> pendingReplyMap = new HashMap<Long, PendingReply>();
        /**
         * Outer call a nested synchronous call was issued from, keyed by the
         * job id of the inner call. A customer carries one job id, so a call
         * made while another is outstanding (A calls B, B calls C, as in a
         * layered network) needs a fresh id, and the outer id must be handed
         * back to the token when the inner reply arrives, or the outer reply
         * finds no pending record and its caller stays blocked forever.
         */
        private final Map<Long, Long> replyParentJobId = new HashMap<Long, Long>();

        // Job ID counter for reply tracking
        private long nextJobId = 0L;

        // Comparators
        private final Comparator<Customer> fcfsComparator = new Comparator<Customer>() {
            @Override public int compare(Customer c1, Customer c2) {
                return Double.compare(c1.orderTime, c2.orderTime);
            }
        };

        private final Comparator<Customer> lcfsComparator = new Comparator<Customer>() {
            @Override public int compare(Customer c1, Customer c2) {
                return Double.compare(c2.orderTime, c1.orderTime);
            }
        };

        private final Comparator<Customer> siroComparator = new Comparator<Customer>() {
            @Override public int compare(Customer c1, Customer c2) {
                return Double.compare(c1.randomRank, c2.randomRank);
            }
        };

        // LINE orders class priorities ascending: lower value = more urgent, 0 highest.
        // (SaveHandlers inverts this when exporting to JMT, which orders them the other
        // way.) Every priority comparator here must follow the ascending convention.
        private final Comparator<Customer> priorityComparator = new Comparator<Customer>() {
            @Override public int compare(Customer c1, Customer c2) {
                if (c1.priority != c2.priority) {
                    return Integer.compare(c1.priority, c2.priority);
                }
                return Double.compare(c1.orderTime, c2.orderTime);
            }
        };

        private final Comparator<Customer> priorityLcfsComparator = new Comparator<Customer>() {
            @Override public int compare(Customer c1, Customer c2) {
                if (c1.priority != c2.priority) {
                    return Integer.compare(c1.priority, c2.priority);
                }
                return Double.compare(c2.orderTime, c1.orderTime);
            }
        };

        private final Comparator<Customer> sjfComparator = new Comparator<Customer>() {
            @Override public int compare(Customer c1, Customer c2) {
                if (c1.serviceTime != c2.serviceTime) {
                    return Double.compare(c1.serviceTime, c2.serviceTime);
                }
                return Double.compare(c1.orderTime, c2.orderTime);
            }
        };

        private final Comparator<Customer> ljfComparator = new Comparator<Customer>() {
            @Override public int compare(Customer c1, Customer c2) {
                if (c1.serviceTime != c2.serviceTime) {
                    return Double.compare(c2.serviceTime, c1.serviceTime);
                }
                return Double.compare(c1.orderTime, c2.orderTime);
            }
        };

        private final Comparator<Customer> eddComparator = new Comparator<Customer>() {
            @Override public int compare(Customer c1, Customer c2) {
                if (c1.absoluteDeadline != c2.absoluteDeadline) {
                    return Double.compare(c1.absoluteDeadline, c2.absoluteDeadline);
                }
                return Double.compare(c1.orderTime, c2.orderTime);
            }
        };

        private final Comparator<Customer> edfComparator = eddComparator;

        /** Triple key used for preemption history lookup (queueIdx, sysArrival, queueArrival). */
        private static final class PreemptionKey {
            final int queueIdx;
            final double systemArrivalTime;
            final double queueArrivalTime;

            PreemptionKey(int queueIdx, double systemArrivalTime, double queueArrivalTime) {
                this.queueIdx = queueIdx;
                this.systemArrivalTime = systemArrivalTime;
                this.queueArrivalTime = queueArrivalTime;
            }

            @Override
            public boolean equals(Object o) {
                if (this == o) return true;
                if (!(o instanceof PreemptionKey)) return false;
                PreemptionKey k = (PreemptionKey) o;
                return queueIdx == k.queueIdx
                        && Double.doubleToLongBits(systemArrivalTime) == Double.doubleToLongBits(k.systemArrivalTime)
                        && Double.doubleToLongBits(queueArrivalTime) == Double.doubleToLongBits(k.queueArrivalTime);
            }

            @Override
            public int hashCode() {
                int h = queueIdx;
                long b1 = Double.doubleToLongBits(systemArrivalTime);
                long b2 = Double.doubleToLongBits(queueArrivalTime);
                h = 31 * h + (int) (b1 ^ (b1 >>> 32));
                h = 31 * h + (int) (b2 ^ (b2 >>> 32));
                return h;
            }
        }

        // SRPT comparator
        private Comparator<Customer> createSRPTComparator(final int queueIdx) {
            return new Comparator<Customer>() {
                @Override public int compare(Customer c1, Customer c2) {
                    PreemptionKey key1 = new PreemptionKey(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime);
                    PreemptionKey key2 = new PreemptionKey(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime);
                    PreemptionRecord r1 = preemptedJobHistory.get(key1);
                    PreemptionRecord r2 = preemptedJobHistory.get(key2);
                    double remaining1 = (r1 != null) ? r1.remainingWork : c1.serviceTime;
                    double remaining2 = (r2 != null) ? r2.remainingWork : c2.serviceTime;
                    if (remaining1 != remaining2) {
                        return Double.compare(remaining1, remaining2);
                    }
                    return Double.compare(c1.queueArrivalTime, c2.queueArrivalTime);
                }
            };
        }

        /**
         * Virtual PS finish time for a single target customer at queueIdx.
         * The PS shadow includes all jobs currently waiting at queueIdx, all jobs
         * currently in service at queueIdx, and the target itself if not already
         * accounted for (so the comparator works during PriorityQueue.offer when
         * the new arrival is not yet visible to iteration).
         */
        private double computeFSPVirtualFinishTime(int queueIdx, Customer target) {
            PriorityQueue<Customer> waiting = waitQueues[queueIdx];
            List<PreemptiveCustomer> inService = preemptiveJobsInService[queueIdx];
            int c = (numServers[queueIdx] > 0) ? numServers[queueIdx] : 1;
            double now = ssjSim.time();

            PreemptionKey targetKey = new PreemptionKey(queueIdx,
                    target.systemArrivalTime, target.queueArrivalTime);
            int estCount = waiting.size() + ((inService != null) ? inService.size() : 0) + 1;
            List<Double> works = new ArrayList<Double>(estCount);
            boolean targetIncluded = false;

            for (Customer cust : waiting) {
                double w = residualWorkForCustomer(queueIdx, cust);
                works.add(w);
                PreemptionKey k = new PreemptionKey(queueIdx,
                        cust.systemArrivalTime, cust.queueArrivalTime);
                if (k.equals(targetKey)) {
                    targetIncluded = true;
                }
            }
            if (inService != null) {
                for (PreemptiveCustomer pc : inService) {
                    double residual = pc.remainingServiceWork - serviceWorkBetween(queueIdx, pc.classId, pc.serviceStartTime, now);
                    if (residual < 0.0) residual = 0.0;
                    works.add(residual);
                    PreemptionKey k = new PreemptionKey(queueIdx,
                            pc.systemArrivalTime, pc.queueArrivalTime);
                    if (k.equals(targetKey)) {
                        targetIncluded = true;
                    }
                }
            }
            double targetWork = residualWorkForCustomer(queueIdx, target);
            if (!targetIncluded) {
                works.add(targetWork);
            }

            Collections.sort(works);
            int n = works.size();
            double cumulative = 0.0;
            double prev = 0.0;
            for (int rank = 0; rank < n; rank++) {
                double w = works.get(rank);
                int jobsRemaining = n - rank;
                double rate = Math.min(1.0, (double) c / jobsRemaining);
                if (rate > 0.0) {
                    cumulative += (w - prev) / rate;
                }
                if (w >= targetWork) {
                    return now + cumulative;
                }
                prev = w;
            }
            return now + cumulative;
        }

        private double residualWorkForCustomer(int queueIdx, Customer cust) {
            PreemptionKey k = new PreemptionKey(queueIdx,
                    cust.systemArrivalTime, cust.queueArrivalTime);
            PreemptionRecord rec = preemptedJobHistory.get(k);
            return (rec != null) ? rec.remainingWork : cust.serviceTime;
        }

        /**
         * Residual virtual finish time for an in-service job at queueIdx.
         */
        private double computeFSPVirtualFinishTimeInService(int queueIdx, PreemptiveCustomer target) {
            PriorityQueue<Customer> waiting = waitQueues[queueIdx];
            List<PreemptiveCustomer> inService = preemptiveJobsInService[queueIdx];
            int c = (numServers[queueIdx] > 0) ? numServers[queueIdx] : 1;
            double now = ssjSim.time();

            int estCount = waiting.size() + ((inService != null) ? inService.size() : 0);
            List<Double> works = new ArrayList<Double>(estCount);
            for (Customer cust : waiting) {
                works.add(residualWorkForCustomer(queueIdx, cust));
            }
            if (inService != null) {
                for (PreemptiveCustomer pc : inService) {
                    double residual = pc.remainingServiceWork - serviceWorkBetween(queueIdx, pc.classId, pc.serviceStartTime, now);
                    if (residual < 0.0) residual = 0.0;
                    works.add(residual);
                }
            }
            double targetResidual = target.remainingServiceWork - serviceWorkBetween(queueIdx, target.classId, target.serviceStartTime, now);
            if (targetResidual < 0.0) targetResidual = 0.0;

            Collections.sort(works);
            int n = works.size();
            double cumulative = 0.0;
            double prev = 0.0;
            for (int rank = 0; rank < n; rank++) {
                double w = works.get(rank);
                int jobsRemaining = n - rank;
                double rate = Math.min(1.0, (double) c / jobsRemaining);
                if (rate > 0.0) {
                    cumulative += (w - prev) / rate;
                }
                if (w >= targetResidual) {
                    return now + cumulative;
                }
                prev = w;
            }
            return now + cumulative;
        }

        // FSP comparator (Fair Sojourn Protocol): rank by virtual PS finish time
        private Comparator<Customer> createFSPComparator(final int queueIdx) {
            return new Comparator<Customer>() {
                @Override public int compare(Customer c1, Customer c2) {
                    double vft1 = computeFSPVirtualFinishTime(queueIdx, c1);
                    double vft2 = computeFSPVirtualFinishTime(queueIdx, c2);
                    if (vft1 != vft2) {
                        return Double.compare(vft1, vft2);
                    }
                    return Double.compare(c1.queueArrivalTime, c2.queueArrivalTime);
                }
            };
        }

        // SRPTPRIO comparator
        private Comparator<Customer> createSRPTPRIOComparator(final int queueIdx) {
            return new Comparator<Customer>() {
                @Override public int compare(Customer c1, Customer c2) {
                    PreemptionKey key1 = new PreemptionKey(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime);
                    PreemptionKey key2 = new PreemptionKey(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime);
                    PreemptionRecord r1 = preemptedJobHistory.get(key1);
                    PreemptionRecord r2 = preemptedJobHistory.get(key2);
                    double remaining1 = (r1 != null) ? r1.remainingWork : c1.serviceTime;
                    double remaining2 = (r2 != null) ? r2.remainingWork : c2.serviceTime;
                    if (c1.priority != c2.priority) {
                        return Integer.compare(c1.priority, c2.priority);
                    }
                    if (remaining1 != remaining2) {
                        return Double.compare(remaining1, remaining2);
                    }
                    return Double.compare(c1.queueArrivalTime, c2.queueArrivalTime);
                }
            };
        }

        private Comparator<Customer> getComparatorForStrategy(SchedStrategy strategy, int queueIdx) {
            switch (strategy) {
                case SIRO:
                    return siroComparator;
                case LCFS:
                case LCFSPR:
                case LCFSPI:
                    return lcfsComparator;
                case HOL:
                case FCFSPRIO:
                case FCFSPR:
                case FCFSPI:
                case FCFSPRPRIO:
                case FCFSPIPRIO:
                case PSPRIO:
                case DPSPRIO:
                case GPSPRIO:
                    return priorityComparator;
                case LCFSPRIO:
                case LCFSPRPRIO:
                case LCFSPIPRIO:
                    return priorityLcfsComparator;
                case SEPT:
                    return createSEPTComparator(queueIdx);
                case LEPT:
                    return createLEPTComparator(queueIdx);
                case SJF:
                    return sjfComparator;
                case LJF:
                    return ljfComparator;
                case EDD:
                    return eddComparator;
                case EDF:
                    return edfComparator;
                case SRPT:
                    return createSRPTComparator(queueIdx);
                case SRPTPRIO:
                    return createSRPTPRIOComparator(queueIdx);
                case PSJF:
                    return createPSJFComparator(queueIdx);
                case FB:
                    return createFBComparator(queueIdx);
                case LRPT:
                    return createLRPTComparator(queueIdx);
                case SETF:
                    return createSETFComparator(queueIdx);
                case FSP:
                    return createFSPComparator(queueIdx);
                default:
                    return fcfsComparator;
            }
        }

        // SETF comparator
        private Comparator<Customer> createSETFComparator(final int queueIdx) {
            return new Comparator<Customer>() {
                @Override public int compare(Customer c1, Customer c2) {
                    PreemptionKey key1 = new PreemptionKey(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime);
                    PreemptionKey key2 = new PreemptionKey(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime);
                    PreemptionRecord r1 = preemptedJobHistory.get(key1);
                    PreemptionRecord r2 = preemptedJobHistory.get(key2);
                    double attained1 = (r1 != null) ? r1.elapsedTime : 0.0;
                    double attained2 = (r2 != null) ? r2.elapsedTime : 0.0;
                    if (attained1 != attained2) {
                        return Double.compare(attained1, attained2);
                    }
                    return Double.compare(c1.queueArrivalTime, c2.queueArrivalTime);
                }
            };
        }

        // PSJF comparator
        private Comparator<Customer> createPSJFComparator(final int queueIdx) {
            return new Comparator<Customer>() {
                @Override public int compare(Customer c1, Customer c2) {
                    PreemptionKey key1 = new PreemptionKey(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime);
                    PreemptionKey key2 = new PreemptionKey(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime);
                    PreemptionRecord r1 = preemptedJobHistory.get(key1);
                    PreemptionRecord r2 = preemptedJobHistory.get(key2);
                    double original1 = (r1 != null) ? r1.originalTotal : c1.serviceTime;
                    double original2 = (r2 != null) ? r2.originalTotal : c2.serviceTime;
                    if (original1 != original2) {
                        return Double.compare(original1, original2);
                    }
                    return Double.compare(c1.queueArrivalTime, c2.queueArrivalTime);
                }
            };
        }

        // LRPT comparator
        private Comparator<Customer> createLRPTComparator(final int queueIdx) {
            return new Comparator<Customer>() {
                @Override public int compare(Customer c1, Customer c2) {
                    PreemptionKey key1 = new PreemptionKey(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime);
                    PreemptionKey key2 = new PreemptionKey(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime);
                    PreemptionRecord r1 = preemptedJobHistory.get(key1);
                    PreemptionRecord r2 = preemptedJobHistory.get(key2);
                    double remaining1 = (r1 != null) ? r1.remainingWork : c1.serviceTime;
                    double remaining2 = (r2 != null) ? r2.remainingWork : c2.serviceTime;
                    if (remaining1 != remaining2) {
                        return Double.compare(remaining2, remaining1);
                    }
                    double original1 = (r1 != null) ? r1.originalTotal : c1.serviceTime;
                    double original2 = (r2 != null) ? r2.originalTotal : c2.serviceTime;
                    return Double.compare(original1, original2);
                }
            };
        }

        // FB / LAS comparator
        private Comparator<Customer> createFBComparator(final int queueIdx) {
            return new Comparator<Customer>() {
                @Override public int compare(Customer c1, Customer c2) {
                    PreemptionKey key1 = new PreemptionKey(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime);
                    PreemptionKey key2 = new PreemptionKey(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime);
                    PreemptionRecord r1 = preemptedJobHistory.get(key1);
                    PreemptionRecord r2 = preemptedJobHistory.get(key2);
                    double attained1 = (r1 != null) ? r1.elapsedTime : 0.0;
                    double attained2 = (r2 != null) ? r2.elapsedTime : 0.0;
                    if (attained1 != attained2) {
                        return Double.compare(attained1, attained2);
                    }
                    return Double.compare(c1.queueArrivalTime, c2.queueArrivalTime);
                }
            };
        }

        // SEPT comparator
        private Comparator<Customer> createSEPTComparator(final int queueIdx) {
            return new Comparator<Customer>() {
                @Override public int compare(Customer c1, Customer c2) {
                    double mu1 = mus[queueIdx][c1.classId];
                    double mu2 = mus[queueIdx][c2.classId];
                    if (mu1 != mu2) {
                        return Double.compare(mu2, mu1);
                    }
                    return Double.compare(c1.queueArrivalTime, c2.queueArrivalTime);
                }
            };
        }

        // LEPT comparator
        private Comparator<Customer> createLEPTComparator(final int queueIdx) {
            return new Comparator<Customer>() {
                @Override public int compare(Customer c1, Customer c2) {
                    double mu1 = mus[queueIdx][c1.classId];
                    double mu2 = mus[queueIdx][c2.classId];
                    if (mu1 != mu2) {
                        return Double.compare(mu1, mu2);
                    }
                    return Double.compare(c1.queueArrivalTime, c2.queueArrivalTime);
                }
            };
        }

        private boolean isPSScheduling(SchedStrategy strategy) {
            switch (strategy) {
                case PS:
                case DPS:
                case GPS:
                case PSPRIO:
                case DPSPRIO:
                case GPSPRIO:
                case LPS:
                    return true;
                default:
                    return false;
            }
        }

        private boolean isPSWithPriority(SchedStrategy strategy) {
            switch (strategy) {
                case PSPRIO:
                case DPSPRIO:
                case GPSPRIO:
                    return true;
                default:
                    return false;
            }
        }

        private boolean isPSWithWeights(SchedStrategy strategy) {
            switch (strategy) {
                case DPS:
                case GPS:
                case DPSPRIO:
                case GPSPRIO:
                    return true;
                default:
                    return false;
            }
        }

        private boolean isPreemptiveLCFSScheduling(SchedStrategy strategy) {
            switch (strategy) {
                case LCFSPR:
                case LCFSPI:
                case LCFSPRPRIO:
                case LCFSPIPRIO:
                    return true;
                default:
                    return false;
            }
        }

        private boolean isPreemptiveFCFSScheduling(SchedStrategy strategy) {
            switch (strategy) {
                case FCFSPR:
                case FCFSPI:
                case FCFSPRPRIO:
                case FCFSPIPRIO:
                    return true;
                default:
                    return false;
            }
        }

        private boolean isLCFSFamily(SchedStrategy strategy) {
            switch (strategy) {
                case LCFS:
                case LCFSPR:
                case LCFSPI:
                case LCFSPRIO:
                case LCFSPRPRIO:
                case LCFSPIPRIO:
                    return true;
                default:
                    return false;
            }
        }

        private boolean isSRPTScheduling(SchedStrategy strategy) {
            switch (strategy) {
                case SRPT:
                case SRPTPRIO:
                    return true;
                default:
                    return false;
            }
        }

        private boolean isPSJFScheduling(SchedStrategy strategy) {
            return strategy == SchedStrategy.PSJF;
        }

        private boolean isFBScheduling(SchedStrategy strategy) {
            return strategy == SchedStrategy.FB;
        }

        private boolean isLRPTScheduling(SchedStrategy strategy) {
            return strategy == SchedStrategy.LRPT;
        }

        private boolean isFSPScheduling(SchedStrategy strategy) {
            return strategy == SchedStrategy.FSP;
        }

        private boolean isEDFScheduling(SchedStrategy strategy) {
            return strategy == SchedStrategy.EDF;
        }

        // Comparator-based preemptive-resume policies (size-based plus deadline-based EDF).
        // EDD is deliberately excluded: it is the non-preemptive earliest-due-date variant.
        private boolean isSizeBasedPreemptiveScheduling(SchedStrategy strategy) {
            return isSRPTScheduling(strategy) || isPSJFScheduling(strategy)
                    || isFBScheduling(strategy) || isLRPTScheduling(strategy)
                    || isFSPScheduling(strategy) || isEDFScheduling(strategy);
        }

        // Random number generators
        private RandomVariateGen[][] arrivalGens;
        private RandomVariateGen[][] serviceGens;
        private RandomVariateGen[][] setupGens;
        private RandomVariateGen[][] delayoffGens;
        private MRG32k3a routingRng;
        private Random siroRng;

        // PH distribution support for service
        private ProcessType[][] serviceProcessType;
        private MatrixCell[][] serviceProc;
        private Random[][] serviceRng;
        // Stateful phase-carrying samplers for correlated (MAP/MMPP2/RAP) service, keyed
        // by [serviceNode][class]; created lazily so that autocorrelation is preserved
        // across successive services (renewal map_sample would reset the phase each call).
        private Map_sample.MapSampler[][] serviceMapSampler;
        // Inverse-CDF samplers for ME service processes, keyed by [serviceNode][class].
        // ME is not a phase-type process, so map_sample's CTMC walk cannot be used;
        // the sampler is cached because it builds an inversion table on construction.
        private Me_sample.MeSampler[][] serviceMeSampler;
        // Conditional-vector samplers for RAP service processes. They carry the
        // conditional phase vector across services, which is what preserves the
        // autocorrelation of a RAP (the role MapSampler plays for a MAP).
        private Rap_sample.RapSampler[][] serviceRapSampler;

        // Batch Markovian Service Process (BMSP) support: a BMAP assigned as a
        // station's service process is interpreted as a bulk (batch) server.
        // Scope: single-server, single-class, FCFS stations (the canonical
        // M/BMSP/1 bulk queue). One station-level clock removes min(k,N) jobs
        // FCFS at each firing (partial-batch truncation), matching the MAM
        // BMAP/BMSP/1 boundary convention (solver_mam_map_bmap_1). The clock is
        // frozen while the station is empty and resumes (phase retained) on the
        // next arrival, mirroring the MAP/MAP/1 service convention.
        private boolean[] isBatchServiceStation;      // [serviceNode]
        private int[] batchServiceClass;              // served class of a batch station (-1 if none)
        private Map_sample.BmapSampler[] batchServiceSampler; // [serviceNode]
        private Random[] batchServiceRng;             // [serviceNode]
        private int[] batchPendingSize;               // batch size k of the scheduled firing

        // PAS (pass-and-swap) RNG: exponential completion times + position/swap selection
        private MRG32k3a pasRng;

        // PH distribution support for arrivals
        private ProcessType[][] arrivalProcessType;
        private MatrixCell[][] arrivalProc;
        private Random[][] arrivalRng;
        // Stateful phase-carrying samplers for correlated (MAP/MMPP2/RAP) arrivals.
        private Map_sample.MapSampler[][] arrivalMapSampler;
        // Inverse-CDF samplers for ME arrival processes (renewal, table cached).
        private Me_sample.MeSampler[][] arrivalMeSampler;
        // Conditional-vector samplers for RAP arrival processes (correlated).
        private Rap_sample.RapSampler[][] arrivalRapSampler;

        // BMAP support
        private int[][] arrivalBatchSize;
        /** Explicit per-class batch-size law from sn.arrivalbatch; null = single arrivals. */
        private DiscreteDistribution[][] arrivalBatchDist;
        private Random[][] arrivalBatchRng;
        // Stateful phase-carrying samplers for correlated (BMAP) batch arrivals.
        private Map_sample.BmapSampler[][] arrivalBmapSampler;

        // MMAP (marked MAP) support: stateful phase-carrying marked samplers,
        // the 1-based mark sampled for the next scheduled arrival (0 = none),
        // and the mark->class binding of each source (sn.markidx):
        // markedGroupClass[srcIdx][m] = classId emitted by mark m (index 0 unused),
        // markedCarrierClass[srcIdx] = class whose ExternalArrival events drive
        // the shared sampler (-1 when the source has no marked group).
        private jline.api.mam.Mmap_sample.MmapSampler[][] arrivalMmapSampler;
        private int[][] arrivalPendingMark;
        private int[][] markedGroupClass;
        private int[] markedCarrierClass;

        // Rate schedule (NHPP) per source/class in breakpoint
        // form; null when the class carries no schedule. Layout per entry:
        // [0] breakpoints (n+1), [1] rates (n), [2] {cyclic ? 1 : 0}.
        private double[][][][] arrivalSchedule;
        // Rate schedule per service-node/class in the same layout; null when none.
        private double[][][][] serviceSchedule;

        // Service node state
        private PriorityQueue<Customer>[] waitQueues;
        private boolean[][] serverBusy;
        private int[] customersInService;

        // Setup and delayoff state tracking
        private ServerState[][] serverState;
        private boolean[] hasSetupDelayoff;
        private int[][] serverLastClass;
        private Event[][] pendingDelayoffEvents;

        /** In-service job tracking record (signal-based removal). */
        private static final class InServiceJob {
            final Customer customer;
            final Event departureEvent;

            InServiceJob(Customer customer, Event departureEvent) {
                this.customer = customer;
                this.departureEvent = departureEvent;
            }
        }

        /** A (queueIdx, serverId) pair used as map key. */
        private static final class IntPair {
            final int a;
            final int b;
            IntPair(int a, int b) { this.a = a; this.b = b; }
            @Override public boolean equals(Object o) {
                if (this == o) return true;
                if (!(o instanceof IntPair)) return false;
                IntPair p = (IntPair) o;
                return a == p.a && b == p.b;
            }
            @Override public int hashCode() {
                return 31 * a + b;
            }
        }

        private final Map<IntPair, InServiceJob> inServiceJobs = new HashMap<IntPair, InServiceJob>();

        /** Delay-node job tracking record. */
        private static final class DelayJob {
            final int queueIdx;
            final Customer customer;
            final Event departureEvent;

            DelayJob(int queueIdx, Customer customer, Event departureEvent) {
                this.queueIdx = queueIdx;
                this.customer = customer;
                this.departureEvent = departureEvent;
            }
        }

        private final Map<Long, DelayJob> delayJobs = new HashMap<Long, DelayJob>();
        private long nextDelayJobId = 0L;

        // ==================== Heterogeneous Server Support ====================
        private int[] numServerTypes;
        private int[][] serversPerType;
        private boolean[][][] serverCompat;
        private int[][] busyCountPerType;
        private int[][] serverToType;
        private HeteroSchedPolicy[] heteroSchedPolicy;
        private RandomVariateGen[][][] heteroServiceGens;
        private double[][][] heteroMus;
        private ProcessType[][][] heteroServiceProcType;
        private MatrixCell[][][] heteroServiceProc;
        private Random[][][] heteroServiceRng;
        // Stateful phase-carrying samplers for correlated heterogeneous service.
        private Map_sample.MapSampler[][][] heteroServiceMapSampler;
        // Inverse-CDF samplers for ME heterogeneous service processes.
        private Me_sample.MeSampler[][][] heteroServiceMeSampler;
        // Conditional-vector samplers for RAP heterogeneous service processes.
        private Rap_sample.RapSampler[][][] heteroServiceRapSampler;
        private List<Integer>[] serverTypeOrder;
        private int[][] alfsOrder;

        // Statistics per service node and class
        private Tally[][] responseTimeTally;
        private List<Double>[][] responseTimeSamples;
        private int[][] completedCustomers;
        private double[][] totalQueueTime;
        private double[][] lastQueueUpdateTime;
        private int[][] currentQueueLength;

        // Markov reward accumulation (active only when the model defines rewards
        // via model.setReward). The reward integral is evaluated on the exact
        // integer joint state, so it is correct for nonlinear rewards (e.g. E[n^2]),
        // unlike the transient-mean sample path exposed by sampleSys.
        private boolean hasReward = false;
        private List<String> rewardNames;
        private jline.lang.reward.RewardFunction[] rewardFnArr;
        private double[] rewardArea;                 // time-integral of r(X(t)) per reward
        private double rewardLastUpdateTime;         // global last flush time
        private double rewardTotalTime;              // total observation window
        private int rewardRowCols;                   // nstations * numClasses
        private Map<String, double[]> rewardStateCache;   // joint-state key -> reward values
        private List<double[]> rewardTranSeries;     // per flush: [t, r0, r1, ...]
        // State residence-time histogram export: lets host languages (MATLAB, native
        // Python) that cannot pass a Java reward function into the simulation evaluate
        // arbitrary rewards on the exact joint-state distribution. key -> [residenceTime, v0, v1, ...].
        private boolean exportStateHistogram = false;
        private java.util.LinkedHashMap<String, double[]> stateHistogram;
        // Integer joint-state trajectory [t, v0, v1, ...] per flush, used by host
        // languages to evaluate the transient reward r(X(t)) along the sample path.
        private List<double[]> stateTranSeries;

        // Utilization tracking per service node and class
        private double[][] totalBusyTime;
        private double[][] lastBusyUpdateTime;
        private int[][] currentBusyServers;

        // Blocking time tracking
        private double[][] totalBlockingTime;
        private int[][] currentBlockedServers;

        // Blocking policy tracking
        private int[][] basBlockedAtDest;

        // DEBUG: BAS event counters
        private int basBlockCount = 0;
        private int basUnblockCount = 0;
        private double q1QueueTimeAtBlock = 0.0;
        private double q1QueueTimeAfterBlock = 0.0;
        private int[][] bbsBlockedAtDest;
        private int[][] fcrBlockedAtDest;

        // Setup time tracking
        private double[][] totalSetupTime;
        private double[][] lastSetupUpdateTime;
        private int[][] currentServersInSetup;

        // Delayoff time tracking
        private double[][] totalDelayoffTime;
        private double[][] lastDelayoffUpdateTime;
        private int[][] currentServersInDelayoff;

        // System-level statistics
        private Tally[] systemResponseTimeTally;
        private Tally[] systemTardinessTally;
        private Tally[][] tardinessTally;
        private int[] systemCompletedCustomers;

        // Dropped customers due to finite buffer
        private int[][] droppedCustomers;

        // Arrived customers at each station (including dropped)
        private int[][] arrivedCustomers;

        // ==================== Impatience Statistics ====================
        private int[][] renegedCustomers;
        private double[][] totalRenegingWaitTime;
        private int[][] balkedCustomers;
        private List<OrbitJob>[] orbitJobs;
        private int[][] retriedCustomers;
        private int[][] maxRetriesExceeded;
        private int[][] currentOrbitSize;
        private double[][] totalOrbitTime;
        private double[][] lastOrbitUpdateTime;

        // Impatience configuration per station-class
        private boolean[][] hasPatienceConfig;
        private RandomVariateGen[][] patienceGens;
        private boolean[][] hasBalkingConfig;
        private boolean[][] hasRetrialConfig;
        private RandomVariateGen[][] retrialGens;
        private int[][] retrialMaxAttemptsConfig;

        /** Triple key for waiting impatient customers (queueIdx, sysArrivalBits, classId). */
        private static final class ImpatientKey {
            final int queueIdx;
            final long systemArrivalBits;
            final int classId;

            ImpatientKey(int queueIdx, long systemArrivalBits, int classId) {
                this.queueIdx = queueIdx;
                this.systemArrivalBits = systemArrivalBits;
                this.classId = classId;
            }

            @Override public boolean equals(Object o) {
                if (this == o) return true;
                if (!(o instanceof ImpatientKey)) return false;
                ImpatientKey k = (ImpatientKey) o;
                return queueIdx == k.queueIdx
                        && systemArrivalBits == k.systemArrivalBits
                        && classId == k.classId;
            }

            @Override public int hashCode() {
                int h = queueIdx;
                h = 31 * h + (int) (systemArrivalBits ^ (systemArrivalBits >>> 32));
                h = 31 * h + classId;
                return h;
            }
        }
        private final Map<ImpatientKey, ImpatientCustomer> waitingImpatientCustomers =
                new HashMap<ImpatientKey, ImpatientCustomer>();

        // PS scheduling state
        private List<PSCustomer>[] psJobsInService;
        private double[] psLastUpdateTime;
        private double[] psLastBusyUpdateTime;

        // PAS (pass-and-swap / order-independent) scheduling state
        private boolean[] isPASStation;
        private List<Customer>[] pasList;          // ordered class list, oldest first
        private Event[] pasDepartureEvent;         // single aggregate completion event
        private double[] pasLastBusyUpdateTime;
        private jline.util.SerializableFunction<jline.util.matrix.Matrix, Double>[] pasSvcRateFun;
        private jline.util.matrix.Matrix[] pasSwapGraph;

        // Preemptive LCFS scheduling state
        private List<PreemptiveCustomer>[] preemptiveJobsInService;
        private boolean[] isPreemptiveScheduling;

        /**
         * Record for storing complete preemption state.
         */
        private static final class PreemptionRecord {
            final double remainingWork;
            final double originalTotal;
            final double elapsedTime;
            final ProcessType distType;
            final Integer phaseParam;

            PreemptionRecord(double remainingWork, double originalTotal, double elapsedTime,
                             ProcessType distType, Integer phaseParam) {
                this.remainingWork = remainingWork;
                this.originalTotal = originalTotal;
                this.elapsedTime = elapsedTime;
                this.distType = distType;
                this.phaseParam = phaseParam;
            }

            PreemptionRecord(double remainingWork, double originalTotal, double elapsedTime) {
                this(remainingWork, originalTotal, elapsedTime, null, null);
            }
        }
        private Map<PreemptionKey, PreemptionRecord> preemptedJobHistory;

        // DPS/GPS weights
        private double[][] schedWeights;

        // LPS limits
        private int[] lpsLimits;

        // Load-dependent service support
        private double[][] lldScaling;
        private boolean[] isLoadDependent;
        private boolean hasLld = false;

        // Class-dependence support: beta_{i,r}(n) per service node, a function of
        // the per-class population vector at the station returning either a scalar
        // (chain-independent) or a length-R vector of per-class RATES.
        private SerializableFunction<Matrix, Matrix>[] cdFunctions = null;
        private boolean hasCd = false;

        // State-dependent departure event tracking (class- and load-dependence)
        private Map<Integer, Event>[] sdDepartureEvents = null;
        private Map<Integer, Customer>[] sdInServiceCustomers = null;

        // Finite capacity region state tracking
        private int[][] currentJobsInRegion;
        /** Memory currently occupied per region: incrementally maintained as
         *  sum over classes of currentJobsInRegion * fcRegionClassSize, so the
         *  per-arrival admission check stays O(1) (numClasses can be large:
         *  cache retrieval algorithms add one class per item). */
        private double[] currentMemInRegion;
        private int[][] droppedByRegion;

        // FCR time-weighted metrics tracking
        private double[][] totalRegionJobTime;        // time-integral of raw in-region job count (per class)
        private double[][] totalRegionWeightTime;     // time-integral of weighted occupation (jobs * classWeight)
        private double[][] totalRegionMemTime;        // time-integral of memory occupation (jobs * classSize)
        private double[] lastRegionUpdateTime;
        private int[][] regionCompletions;
        private Tally[][] regionResponseTimeTally;

        // FCR class weights
        private double[][] fcRegionClassWeights;
        private double[][] fcRegionClassSize;   // per-class memory footprint for the region memory constraint (JMT classSize)

        // FCR arrival rate tracking
        private double[][] lastRegionArrivalTime;
        private int[][] regionArrivalCount;
        private double[][] regionInterArrivalTimeSum;

        /** Record for an FCR-blocked customer along with its destination queue. */
        private static final class BlockedCustomer {
            final Customer customer;
            final int destQueueIdx;

            BlockedCustomer(Customer customer, int destQueueIdx) {
                this.customer = customer;
                this.destQueueIdx = destQueueIdx;
            }
        }

        private LinkedList<BlockedCustomer>[] fcRegionBlockedQueue;
        private int[][] blockedInRegion;

        // Trace writer for DEBUG verbose level
        private PrintWriter traceWriter = null;
        private boolean traceEnabled = false;

        // Logger node support
        private final Map<Integer, BufferedWriter> loggerWriters = new HashMap<Integer, BufferedWriter>();
        private final Map<Integer, LoggerConfig> loggerConfigs = new HashMap<Integer, LoggerConfig>();
        private final Map<Integer, double[]> loggerLastJobTimePerClass = new HashMap<Integer, double[]>();
        private final Map<Integer, Double> loggerLastJobTimeAny = new HashMap<Integer, Double>();
        private String simulationStartTime = "";

        /** Logger configuration record. */
        private static final class LoggerConfig {
            final int nodeIdx;
            final String fileName;
            final String filePath;
            final boolean logLoggerName;
            final boolean logTimestamp;
            final boolean logJobID;
            final boolean logJobClass;
            final boolean logTimeSameClass;
            final boolean logTimeAnyClass;
            final boolean logStartTime;
            final String loggerName;
            final String delimiter;
            final String decimalSeparator;

            LoggerConfig(int nodeIdx, String fileName, String filePath,
                         boolean logLoggerName, boolean logTimestamp, boolean logJobID,
                         boolean logJobClass, boolean logTimeSameClass,
                         boolean logTimeAnyClass, boolean logStartTime,
                         String loggerName, String delimiter, String decimalSeparator) {
                this.nodeIdx = nodeIdx;
                this.fileName = fileName;
                this.filePath = filePath;
                this.logLoggerName = logLoggerName;
                this.logTimestamp = logTimestamp;
                this.logJobID = logJobID;
                this.logJobClass = logJobClass;
                this.logTimeSameClass = logTimeSameClass;
                this.logTimeAnyClass = logTimeAnyClass;
                this.logStartTime = logStartTime;
                this.loggerName = loggerName;
                this.delimiter = delimiter;
                this.decimalSeparator = decimalSeparator;
            }

            LoggerConfig(int nodeIdx, String fileName, String filePath,
                         boolean logLoggerName, boolean logTimestamp, boolean logJobID,
                         boolean logJobClass, boolean logTimeSameClass,
                         boolean logTimeAnyClass, boolean logStartTime,
                         String loggerName) {
                this(nodeIdx, fileName, filePath, logLoggerName, logTimestamp, logJobID,
                        logJobClass, logTimeSameClass, logTimeAnyClass, logStartTime,
                        loggerName, ",", ".");
            }
        }

        /** Trace sampler for Replayer distribution - cycles through trace data sequentially. */
        private static final class TraceSampler {
            private final double[] data;
            private int index = 0;

            TraceSampler(double[] data) {
                this.data = data;
            }

            double nextSample() {
                if (data.length == 0) return 1e-9;
                double value = data[index];
                index = (index + 1) % data.length;
                return (value <= 0) ? 1e-9 : value;
            }

            void reset() {
                index = 0;
            }
        }

        // Trace samplers
        private TraceSampler[][] arrivalTraceSamplers;
        private TraceSampler[][] serviceTraceSamplers;

        /** Routing result: destination node and (possibly switched) class. */
        private static final class RoutingResult {
            final int destNode;
            final int destClassId;

            RoutingResult(int destNode, int destClassId) {
                this.destNode = destNode;
                this.destClassId = destClassId;
            }
        }

        /** Customer data class. */
        private static final class Customer {
            final int classId;
            final int priority;
            final double systemArrivalTime;
            final double queueArrivalTime;
            final double randomRank;
            double serviceTime;
            final long jobId;
            double absoluteDeadline;
            int assignedServerType;
            ForkedJob forkedJob;
            // Buffer-ordering key: equals queueArrivalTime except for customers
            // released from an FCR blocked FIFO, whose queueArrivalTime is kept
            // at the original arrival instant for response-time statistics but
            // who join the destination buffer at their release instant (JMT
            // semantics: an admission-gate overtaker is served first)
            double orderTime;

            Customer(int classId, int priority, double systemArrivalTime,
                     double queueArrivalTime, double randomRank, double serviceTime,
                     long jobId, double absoluteDeadline, int assignedServerType,
                     ForkedJob forkedJob) {
                this.classId = classId;
                this.priority = priority;
                this.systemArrivalTime = systemArrivalTime;
                this.queueArrivalTime = queueArrivalTime;
                this.randomRank = randomRank;
                this.serviceTime = serviceTime;
                this.jobId = jobId;
                this.absoluteDeadline = absoluteDeadline;
                this.assignedServerType = assignedServerType;
                this.forkedJob = forkedJob;
                this.orderTime = queueArrivalTime;
            }

            Customer(int classId, int priority, double systemArrivalTime,
                     double queueArrivalTime, double randomRank) {
                this(classId, priority, systemArrivalTime, queueArrivalTime, randomRank,
                        -1.0, -1L, Double.POSITIVE_INFINITY, -1, null);
            }

            Customer(int classId, int priority, double systemArrivalTime,
                     double queueArrivalTime, double randomRank, double serviceTime) {
                this(classId, priority, systemArrivalTime, queueArrivalTime, randomRank,
                        serviceTime, -1L, Double.POSITIVE_INFINITY, -1, null);
            }
        }

        /** PS Customer data class for Processor Sharing scheduling. */
        private static final class PSCustomer {
            final int classId;
            final int priority;
            final double systemArrivalTime;
            final double queueArrivalTime;
            final double totalServiceRequirement;
            double remainingServiceWork;
            Event scheduledDepartureEvent;
            final ForkedJob forkedJob;
            int assignedServerType;

            PSCustomer(int classId, int priority, double systemArrivalTime,
                       double queueArrivalTime, double totalServiceRequirement,
                       double remainingServiceWork, Event scheduledDepartureEvent,
                       ForkedJob forkedJob, int assignedServerType) {
                this.classId = classId;
                this.priority = priority;
                this.systemArrivalTime = systemArrivalTime;
                this.queueArrivalTime = queueArrivalTime;
                this.totalServiceRequirement = totalServiceRequirement;
                this.remainingServiceWork = remainingServiceWork;
                this.scheduledDepartureEvent = scheduledDepartureEvent;
                this.forkedJob = forkedJob;
                this.assignedServerType = assignedServerType;
            }

            PSCustomer(int classId, int priority, double systemArrivalTime,
                       double queueArrivalTime, double totalServiceRequirement,
                       double remainingServiceWork, Event scheduledDepartureEvent) {
                this(classId, priority, systemArrivalTime, queueArrivalTime,
                        totalServiceRequirement, remainingServiceWork,
                        scheduledDepartureEvent, null, -1);
            }
        }

        /** Preemptive LCFS Customer data class. */
        private static final class PreemptiveCustomer {
            final int classId;
            final int priority;
            final double systemArrivalTime;
            final double queueArrivalTime;
            final double randomRank;
            final double totalServiceRequirement;
            double remainingServiceWork;
            double elapsedServiceTime;
            double serviceStartTime;
            Event scheduledDepartureEvent;
            final int serverId;
            int assignedServerType;
            final double absoluteDeadline;

            PreemptiveCustomer(int classId, int priority, double systemArrivalTime,
                               double queueArrivalTime, double randomRank,
                               double totalServiceRequirement, double remainingServiceWork,
                               double elapsedServiceTime, double serviceStartTime,
                               Event scheduledDepartureEvent, int serverId,
                               int assignedServerType, double absoluteDeadline) {
                this.classId = classId;
                this.priority = priority;
                this.systemArrivalTime = systemArrivalTime;
                this.queueArrivalTime = queueArrivalTime;
                this.randomRank = randomRank;
                this.totalServiceRequirement = totalServiceRequirement;
                this.remainingServiceWork = remainingServiceWork;
                this.elapsedServiceTime = elapsedServiceTime;
                this.serviceStartTime = serviceStartTime;
                this.scheduledDepartureEvent = scheduledDepartureEvent;
                this.serverId = serverId;
                this.assignedServerType = assignedServerType;
                this.absoluteDeadline = absoluteDeadline;
            }

            PreemptiveCustomer(int classId, int priority, double systemArrivalTime,
                               double queueArrivalTime, double randomRank,
                               double totalServiceRequirement, double remainingServiceWork,
                               double elapsedServiceTime, double serviceStartTime,
                               Event scheduledDepartureEvent, int serverId) {
                this(classId, priority, systemArrivalTime, queueArrivalTime, randomRank,
                        totalServiceRequirement, remainingServiceWork, elapsedServiceTime,
                        serviceStartTime, scheduledDepartureEvent, serverId, -1,
                        Double.POSITIVE_INFINITY);
            }
        }

        // ==================== Impatience Data Classes ====================

        /**
         * Impatient customer tracking for reneging (timer-based abandonment).
         */
        private static final class ImpatientCustomer {
            final Customer customer;
            final double patienceTime;
            final double patienceDeadline;
            Event renegingEvent;

            ImpatientCustomer(Customer customer, double patienceTime,
                              double patienceDeadline, Event renegingEvent) {
                this.customer = customer;
                this.patienceTime = patienceTime;
                this.patienceDeadline = patienceDeadline;
                this.renegingEvent = renegingEvent;
            }

            ImpatientCustomer(Customer customer, double patienceTime, double patienceDeadline) {
                this(customer, patienceTime, patienceDeadline, null);
            }
        }

        /**
         * Orbit job tracking for retrial (customer retry after rejection).
         */
        private static final class OrbitJob {
            final Customer customer;
            final double orbitEntryTime;
            final int retrialAttempts;
            final int maxAttempts;
            final int sourceQueueIdx;
            Event retrialEvent;

            OrbitJob(Customer customer, double orbitEntryTime, int retrialAttempts,
                     int maxAttempts, int sourceQueueIdx, Event retrialEvent) {
                this.customer = customer;
                this.orbitEntryTime = orbitEntryTime;
                this.retrialAttempts = retrialAttempts;
                this.maxAttempts = maxAttempts;
                this.sourceQueueIdx = sourceQueueIdx;
                this.retrialEvent = retrialEvent;
            }

            OrbitJob(Customer customer, double orbitEntryTime, int retrialAttempts,
                     int maxAttempts, int sourceQueueIdx) {
                this(customer, orbitEntryTime, retrialAttempts, maxAttempts, sourceQueueIdx, null);
            }
        }

        /** Warm-up period as fraction of total simulation time (fallback if MSER-5 not used). */
        private double effectiveWarmupFraction = LDESOptions.DEFAULT_WARMUP_FRAC;

        /** Flag indicating whether warmup period is complete. */
        private boolean warmupDone = false;
        private int preemptDebugCount = 0;

        /** Time when warmup ended (for calculating actual simulation time). */
        private double warmupEndTime = 0.0;

        // MSER-5 transient detection
        /** Batch size for MSER (configurable via options.mserbatch). */
        private int effectiveMserBatchSize = LDESOptions.DEFAULT_MSER_BATCH;

        // Configurable CI parameters (read from LDESOptions)
        /** CI method: "obm", "bm", or "none". */
        private String effectiveCiMethod = "obm";

        /** OBM overlap fraction (0.0-1.0). */
        private double effectiveObmOverlap = LDESOptions.DEFAULT_OBM_OVERLAP;

        /** Minimum batch size for CI. */
        private int effectiveCiMinBatch = LDESOptions.DEFAULT_CI_MIN_BATCH;

        /** Minimum observations for CI. */
        private int effectiveCiMinObs = LDESOptions.DEFAULT_CI_MIN_OBS;

        /** Fraction of lowest frequencies for spectral log-periodogram regression. */
        private double effectiveSpectralLowFreqFrac = LDESOptions.DEFAULT_SPECTRAL_LOW_FREQ_FRAC;

        /** Sampling interval for collecting queue length observations. */
        private double mserSamplingInterval = 0.0;

        /** Queue length observations for MSER-5: [svcIdx][classIdx]. */
        private List<Double>[][] queueLengthObservations;

        /** Throughput observations for MSER-5. */
        private List<Integer>[][] throughputObservations;

        /** Busy time observations for MSER-5. */
        private List<Double>[][] busyTimeObservations;

        /** Cumulative queue-time (integral of QLen) observations for MSER-5. */
        private List<Double>[][] queueTimeObservations;

        /** Blocking time observations for MSER-5. */
        private List<Double>[][] blockingTimeObservations;

        /** Place completion observations for MSER-5. */
        private List<Integer>[][] placeCompletionObservations;

        /** Time stamps for observations. */
        private List<Double> observationTimes;

        /** MSER-5 truncation point (batch index) determined at end of simulation. */
        private int mserTruncationBatch = 0;

        /** Flag indicating if MSER-5 is enabled. */
        private boolean mserEnabled = true;

        // Event-count based stopping
        /** Total number of service completions (event count). */
        private long totalEventCount = 0L;

        /**
         * True while the initial closed-class population is being placed (and,
         * for queueing-Petri-net models, while immediate transitions fire to
         * settle the initial marking). Firings during this phase are model
         * setup, not simulation events: they must not increment the event
         * budget nor trigger warmup/MSER/convergence sampling, whose data
         * structures (e.g. observationTimes) are allocated only afterwards.
         */
        private boolean initializing = false;

        /** Maximum service completion events to simulate (from options.samples). */
        private long maxEvents = 0L;

        // Total simulation event stopping
        /** Total number of all simulation events. */
        private long totalSimEvents = 0L;

        /** Maximum total simulation events (-1 = no limit). */
        private long maxSimEventsLimit = -1L;

        /** Wall-clock time budget in seconds (Double.POSITIVE_INFINITY = no limit). */
        private double maxTimeLimit = Double.POSITIVE_INFINITY;

        /** System.nanoTime() captured when the simulation event loop started. */
        private long simStartNanos = System.nanoTime();

        /** Event count threshold for warmup completion. */
        private long warmupEventThreshold = 0L;

        /** Event count at last MSER sample. */
        private long lastMserEventCount = 0L;

        /** Event interval between MSER samples. */
        private long mserEventInterval = 0L;

        /** Event count at last streaming push. */
        private long lastStreamEventCount = 0L;

        /** Last simulation time when streaming was pushed. */
        private double lastStreamTime = 0.0;

        // ==================== Variance Reduction Support ====================

        /** Flag indicating if antithetic variates are enabled. */
        private boolean useAntitheticVariates = false;

        /** Flag indicating if control variates are enabled. */
        private boolean useControlVariates = false;

        /** Antithetic arrival generators: [sourceIdx][classIdx]. */
        private RandomVariateGen[][] antitheticArrivalGens;

        /** Antithetic service generators: [svcIdx][classIdx]. */
        private RandomVariateGen[][] antitheticServiceGens;

        /** Antithetic RNG for PH/MAP arrivals. */
        private Random[][] antitheticArrivalRng;

        /** Antithetic RNG for PH/MAP services. */
        private Random[][] antitheticServiceRng;

        // Control variates tracking
        /** Sum of sampled arrival times. */
        private double[][] arrivalSampleSum;

        /** Count of arrival samples. */
        private long[][] arrivalSampleCount;

        /** Expected mean arrival time. */
        private double[][] arrivalExpectedMean;

        /** Sum of sampled service times. */
        private double[][] serviceSampleSum;

        /** Count of service samples. */
        private long[][] serviceSampleCount;

        /** Expected mean service time. */
        private double[][] serviceExpectedMean;

        // Convergence detection state
        /** Flag indicating if convergence checking is enabled. */
        private boolean convergenceEnabled = true;

        /** Convergence tolerance. */
        private double convergenceTolerance = 0.05;

        /** Minimum number of batches before checking convergence. */
        private int convergenceMinBatches = 20;

        /** Number of events between convergence checks. */
        private long convergenceCheckInterval = 0L;

        /** Event count at last convergence check. */
        private long lastConvergenceCheckEventCount = 0L;

        /** Batch statistics for queue lengths. */
        private List<Double>[][] queueBatchMeans;

        /** Batch statistics for utilizations. */
        private List<Double>[][] utilBatchMeans;

        /** Batch statistics for response times. */
        private List<Double>[][] respTimeBatchMeans;

        /** Batch statistics for throughputs. */
        private List<Double>[][] throughputBatchMeans;

        /** Current batch observation count within a batch. */
        private int currentBatchObservations = 0;

        /** Cumulative queue time at start of current batch. */
        private double[][] batchStartQueueTime;

        /** Cumulative busy time at start of current batch. */
        private double[][] batchStartBusyTime;

        /** Completed customers at start of current batch. */
        private int[][] batchStartCompletions;

        /** Response time sum within current batch. */
        private double[][] currentBatchRespTimeSum;

        /** Response time count within current batch. */
        private int[][] currentBatchRespTimeCount;

        /** Time at start of current batch. */
        private double batchStartTime = 0.0;

        /** Flag indicating if simulation has converged. */
        private boolean hasConverged = false;

        /** Stopping reason. */
        private String stoppingReason = "max_events";

        /** Final CI half-widths for queue lengths. */
        private double[][] finalQNCI;

        /** Final CI half-widths for utilizations. */
        private double[][] finalUNCI;

        /** Final CI half-widths for response times. */
        private double[][] finalRNCI;

        /** Final CI half-widths for throughputs. */
        private double[][] finalTNCI;

        /** Final relative precision for queue lengths. */
        private double[][] finalQNRelPrec;

        /** Final relative precision for utilizations. */
        private double[][] finalUNRelPrec;

        /** Final relative precision for response times. */
        private double[][] finalRNRelPrec;

        /** Final relative precision for throughputs. */
        private double[][] finalTNRelPrec;

        // =====================================================================
        // Constructor — corresponds to the Kotlin `init {}` block in part.
        // The full `init` block continues to translate node classification logic
        // in part 2.
        // =====================================================================

        SSJSimulator(NetworkStruct sn, SolverOptions options, Matrix initSol, Collector stream) {
            this.sn = sn;
            this.options = options;
            this.initSol = initSol;
            this.stream = stream;

            this.numClasses = sn.nclasses;
            this.numStations = sn.nstations;
            this.numNodes = sn.nnodes;
            this.seed = (long) options.seed;

            // Discrete-time configuration. Read here rather than in simulateBody
            // because SimEvent assigns its slot phase at construction time and the
            // generators are validated against the lattice as they sample.
            if (options instanceof LDESOptions) {
                LDESOptions ldesOpts = (LDESOptions) options;
                this.slotted = ldesOpts.slotted;
                this.slotLength = ldesOpts.slotLength;
            } else {
                this.slotted = false;
                this.slotLength = LDESOptions.DEFAULT_SLOT_LENGTH;
            }
            if (this.slotted) {
                assertSlottedModelIsSupported(sn);
            }

            // Class priorities (from sn.classprio)
            this.classPrio = new int[numClasses];
            for (int k = 0; k < numClasses; k++) {
                this.classPrio[k] = (int) sn.classprio.get(k);
            }

            // Class deadlines (from sn.classdeadline)
            this.classDeadline = new double[numClasses];
            for (int k = 0; k < numClasses; k++) {
                if (sn.classdeadline != null) {
                    this.classDeadline[k] = sn.classdeadline.get(k);
                } else {
                    this.classDeadline[k] = Double.POSITIVE_INFINITY;
                }
            }

            // Finite capacity regions
            this.numRegions = sn.nregions;
            this.fcRegionIndices = new ArrayList<Integer>(numStations);
            for (int i = 0; i < numStations; i++) {
                this.fcRegionIndices.add(-1);
            }
            this.fcRegionGlobalMax = new ArrayList<Integer>(numRegions);
            for (int i = 0; i < numRegions; i++) {
                this.fcRegionGlobalMax.add(Integer.MAX_VALUE);
            }
            this.fcRegionGlobalMaxMem = new double[numRegions];
            for (int i = 0; i < numRegions; i++) {
                this.fcRegionGlobalMaxMem[i] = -1.0;  // unbounded until a member station sets it
            }
            this.fcRegionClassMax = new Matrix(numRegions, numClasses);
            this.fcRegionDropRule = new boolean[numRegions][numClasses];
            for (int f = 0; f < numRegions; f++) {
                for (int r = 0; r < numClasses; r++) {
                    this.fcRegionDropRule[f][r] = true;
                }
            }
            this.fcRegionLinConA = new double[numRegions][0][];
            this.fcRegionLinConb = new double[numRegions][0];

            initFromNetworkStructPart2();
        }

        // =====================================================================
        // Public surface methods (stubs to be filled in by subsequent chunks).
        // The real implementations live in chunks 2-9.
        // =====================================================================

        void simulate(double maxEvents) {
            simulateBody(maxEvents);
        }

        void simulateTransient(double endTime) {
            simulateTransientBody(endTime);
        }

        @SuppressWarnings("unchecked")
        LDESResult getLDESResult() {
            Matrix QN = new Matrix(numStations, numClasses);
            QN.fill(0.0);
            Matrix UN = new Matrix(numStations, numClasses);
            UN.fill(0.0);
            Matrix RN = new Matrix(numStations, numClasses);
            RN.fill(0.0);
            Matrix TN = new Matrix(numStations, numClasses);
            TN.fill(0.0);
            Matrix AN = new Matrix(numStations, numClasses);
            AN.fill(0.0);
            Matrix TardN = new Matrix(numStations, numClasses);
            TardN.fill(0.0);
            Matrix SysTardN = new Matrix(1, numClasses);
            SysTardN.fill(0.0);
            Matrix CN = new Matrix(1, numClasses);
            CN.fill(0.0);
            Matrix XN = new Matrix(1, numClasses);
            XN.fill(0.0);

            // Build set of classes that can receive jobs via class-switching
            Set<Integer> classSwitchClasses = new HashSet<Integer>();
            for (int c = 0; c < sn.nchains; c++) {
                List<Integer> classesInChain = new ArrayList<Integer>();
                for (int k = 0; k < numClasses; k++) {
                    if (sn.chains.get(c, k) > 0) {
                        classesInChain.add(k);
                    }
                }
                if (classesInChain.size() > 1) {
                    boolean chainHasJobs = false;
                    for (Integer k : classesInChain) {
                        if (sn.njobs.get(k) > 0) {
                            chainHasJobs = true;
                            break;
                        }
                    }
                    if (chainHasJobs) {
                        classSwitchClasses.addAll(classesInChain);
                    }
                }
            }

            // A spawn target receives jobs at its station even inside a
            // zero-population chain (LQN phase-2 continuations): unmask its
            // whole chain, which also covers its downstream switches.
            if (spawnClassOf != null) {
                for (int k = 0; k < spawnClassOf.length; k++) {
                    int sc = spawnClassOf[k];
                    if (sc < 0) continue;
                    for (int c = 0; c < sn.nchains; c++) {
                        if (sn.chains.get(c, sc) > 0) {
                            for (int k2 = 0; k2 < numClasses; k2++) {
                                if (sn.chains.get(c, k2) > 0) {
                                    classSwitchClasses.add(k2);
                                }
                            }
                        }
                    }
                }
            }

            // Finalize time-weighted queue length statistics for all service nodes
            double finalTime = ssjSim.time();
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    double elapsed = finalTime - lastQueueUpdateTime[qIdx][k];
                    if (elapsed > 0) {
                        double effectiveQueueLength = currentQueueLength[qIdx][k]
                                + currentBlockedServers[qIdx][k]
                                + basBlockedAtDest[qIdx][k]
                                + bbsBlockedAtDest[qIdx][k]
                                + fcrBlockedAtDest[qIdx][k];
                        totalQueueTime[qIdx][k] += effectiveQueueLength * elapsed;
                        lastQueueUpdateTime[qIdx][k] = finalTime;
                    }
                    double busyElapsed = finalTime - lastBusyUpdateTime[qIdx][k];
                    if (busyElapsed > 0) {
                        totalBusyTime[qIdx][k] += currentBusyServers[qIdx][k] * busyElapsed;
                        lastBusyUpdateTime[qIdx][k] = finalTime;
                    }
                }
            }

            // Finalize the Markov reward integral over the last interval.
            if (hasReward) {
                updateRewardStats();
            }

            // Extract results for each service node (Queue or Delay)
            for (int svcIdx = 0; svcIdx < serviceStations.size(); svcIdx++) {
                int serviceStation = serviceStations.get(svcIdx);
                for (int k = 0; k < numClasses; k++) {
                    double njobs = sn.njobs.get(k);
                    boolean isClosedWithZeroPopulation = !isOpenClass[k]
                            && Double.isFinite(njobs) && njobs == 0.0;
                    boolean canReceiveViaClassSwitch = classSwitchClasses.contains(k);
                    if ((!isClosedWithZeroPopulation || canReceiveViaClassSwitch)
                            && mus[svcIdx][k] < Double.MAX_VALUE) {
                        QN.set(serviceStation, k, getAvgQueueLength(svcIdx, k));
                        UN.set(serviceStation, k, getUtilization(svcIdx, k));
                        RN.set(serviceStation, k, getAvgResponseTime(svcIdx, k));
                        TN.set(serviceStation, k, getThroughput(svcIdx, k));
                        AN.set(serviceStation, k, getArrivalRate(svcIdx, k));
                        double tard = (tardinessTally[svcIdx][k].numberObs() > 0)
                                ? tardinessTally[svcIdx][k].average() : 0.0;
                        TardN.set(serviceStation, k, tard);
                    }
                }
            }

            // Extract results for Place nodes (Petri nets)
            double simTime = getActualSimTime();
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                int placeNodeIdx = placeNodes.get(placeListIdx);
                int placeStationIdx = (int) sn.nodeToStation.get(placeNodeIdx);
                if (placeStationIdx >= 0 && placeStationIdx < numStations) {
                    double currentTime = ssjSim.time();
                    double elapsed = currentTime - lastPlaceUpdateTime[placeListIdx];
                    if (elapsed > 0) {
                        for (int k = 0; k < numClasses; k++) {
                            totalPlaceTokenTime[placeListIdx][k] +=
                                    placeTokens[placeListIdx][k] * elapsed;
                        }
                    }
                    double transitElapsed = currentTime - lastPlaceTransitUpdateTime[placeListIdx];
                    if (transitElapsed > 0) {
                        for (int k = 0; k < numClasses; k++) {
                            placeTransitTokenTime[placeListIdx][k] +=
                                    placeTokensInTransit[placeListIdx][k] * transitElapsed;
                        }
                    }
                    for (int k = 0; k < numClasses; k++) {
                        double qlen = (simTime > 0)
                                ? (totalPlaceTokenTime[placeListIdx][k]
                                   + placeTransitTokenTime[placeListIdx][k]) / simTime
                                : 0.0;
                        QN.set(placeStationIdx, k, qlen);

                        double tput;
                        if (mserEnabled && placeCompletionObservations != null
                                && placeListIdx < placeCompletionObservations.length
                                && k < placeCompletionObservations[placeListIdx].length
                                && observationTimes != null) {
                            List<Integer> compObs = placeCompletionObservations[placeListIdx][k];
                            int truncationIdxLoc = mserTruncationBatch * effectiveMserBatchSize;
                            if (truncationIdxLoc < compObs.size()
                                    && truncationIdxLoc < observationTimes.size()) {
                                int startCompletions = compObs.get(truncationIdxLoc);
                                long endCompletions = placeCompletions[placeListIdx][k];
                                double startTime = observationTimes.get(truncationIdxLoc);
                                double endTime = ssjSim.time();
                                double elapsed2 = endTime - startTime;
                                tput = (elapsed2 > 0)
                                        ? (endCompletions - startCompletions) / elapsed2 : 0.0;
                            } else if (simTime > 0) {
                                tput = (double) placeCompletions[placeListIdx][k] / simTime;
                            } else {
                                tput = 0.0;
                            }
                        } else if (simTime > 0) {
                            tput = (double) placeCompletions[placeListIdx][k] / simTime;
                        } else {
                            tput = 0.0;
                        }
                        TN.set(placeStationIdx, k, tput);
                        // Mean token sojourn time at the place via Little's law
                        // (E[N]/X), matching JMT's per-place response time.
                        double respTime = (tput > 0) ? qlen / tput : 0.0;
                        RN.set(placeStationIdx, k, respTime);
                    }

                    // Utilization of a queueing place's embedded queue. Flush the
                    // busy-server integral to the final clock, then normalize: for a
                    // finite c-server queue U = busyTime/(simTime*c); for INF each token
                    // is always in service so U = busyTime/simTime = mean number in
                    // service (matching the LINE INF convention U == Q).
                    if (isQueueingPlace[placeListIdx]) {
                        accruePlaceBusy(placeListIdx);
                        int servers = placeNumServers[placeListIdx];
                        boolean infinite = (servers == Integer.MAX_VALUE);
                        for (int k = 0; k < numClasses; k++) {
                            double util = 0.0;
                            if (simTime > 0) {
                                double denom = infinite ? simTime : (simTime * servers);
                                util = totalPlaceBusyTime[placeListIdx][k] / denom;
                            }
                            UN.set(placeStationIdx, k, util);
                        }
                    }
                }
            }

            // Extract results for Join nodes
            Matrix dropRateJoin = null;
            if (!joinStations.isEmpty()) {
                dropRateJoin = new Matrix(numStations, numClasses);
                dropRateJoin.fill(0.0);
            }
            for (int joinListIdx = 0; joinListIdx < joinStations.size(); joinListIdx++) {
                int joinStationIdx = joinStations.get(joinListIdx);
                if (joinStationIdx >= 0 && joinStationIdx < numStations) {
                    double currentTime = ssjSim.time();
                    for (int k = 0; k < numClasses; k++) {
                        double elapsed = currentTime - lastJoinUpdateTime[joinListIdx][k];
                        if (elapsed > 0) {
                            totalJoinQueueTime[joinListIdx][k] +=
                                    currentJoinQueueLength[joinListIdx][k] * elapsed;
                        }
                        double qlen = (simTime > 0)
                                ? totalJoinQueueTime[joinListIdx][k] / simTime : 0.0;
                        QN.set(joinStationIdx, k, qlen);
                        double tput = (simTime > 0)
                                ? (double) joinCompletions[joinListIdx][k] / simTime : 0.0;
                        TN.set(joinStationIdx, k, tput);
                        double respTime = (joinResponseTimeTally[joinListIdx][k].numberObs() > 0)
                                ? joinResponseTimeTally[joinListIdx][k].average() : 0.0;
                        RN.set(joinStationIdx, k, respTime);
                        UN.set(joinStationIdx, k, 0.0);
                        double arvR = (simTime > 0)
                                ? (double) arrivedAtJoin[joinListIdx][k] / simTime : 0.0;
                        AN.set(joinStationIdx, k, arvR);
                        // Quorum/PARTIAL sibling-drop rate for getAvgLossTable.
                        if (dropRateJoin != null) {
                            dropRateJoin.set(joinStationIdx, k, (simTime > 0)
                                    ? (double) droppedByJoin[joinListIdx][k] / simTime : 0.0);
                        }
                    }
                }
            }

            // Extract results for each source
            for (int srcIdx = 0; srcIdx < sourceStations.size(); srcIdx++) {
                int sourceStation = sourceStations.get(srcIdx);
                for (int k = 0; k < numClasses; k++) {
                    if (lambdas[srcIdx][k] > 0) {
                        TN.set(sourceStation, k, lambdas[srcIdx][k]);
                    }
                }
            }

            // System metrics
            for (int k = 0; k < numClasses; k++) {
                if (isOpenClass[k]) {
                    double totalArrival = 0.0;
                    for (int srcIdx = 0; srcIdx < sourceStations.size(); srcIdx++) {
                        totalArrival += lambdas[srcIdx][k];
                    }
                    if (totalArrival > 0) {
                        XN.set(0, k, getSystemThroughput(k));
                        CN.set(0, k, getSystemResponseTime(k));
                        double sysTard = (systemTardinessTally[k].numberObs() > 0)
                                ? systemTardinessTally[k].average() : 0.0;
                        SysTardN.set(0, k, sysTard);
                    }
                } else {
                    boolean foundTput = false;
                    for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                        double tput = getThroughput(svcIdx, k);
                        if (tput > 0) {
                            XN.set(0, k, tput);
                            CN.set(0, k, getSystemResponseTime(k));
                            double sysTard = (systemTardinessTally[k].numberObs() > 0)
                                    ? systemTardinessTally[k].average() : 0.0;
                            SysTardN.set(0, k, sysTard);
                            foundTput = true;
                            break;
                        }
                    }
                    if (!foundTput && !placeNodes.isEmpty()) {
                        for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                            double tput;
                            if (mserEnabled && placeCompletionObservations != null
                                    && placeListIdx < placeCompletionObservations.length
                                    && k < placeCompletionObservations[placeListIdx].length
                                    && observationTimes != null) {
                                List<Integer> compObs = placeCompletionObservations[placeListIdx][k];
                                int truncationIdxLoc = mserTruncationBatch * effectiveMserBatchSize;
                                if (truncationIdxLoc < compObs.size()
                                        && truncationIdxLoc < observationTimes.size()) {
                                    int startCompletions = compObs.get(truncationIdxLoc);
                                    long endCompletions = placeCompletions[placeListIdx][k];
                                    double startTime = observationTimes.get(truncationIdxLoc);
                                    double endTime = ssjSim.time();
                                    double elapsed2 = endTime - startTime;
                                    tput = (elapsed2 > 0)
                                            ? (endCompletions - startCompletions) / elapsed2 : 0.0;
                                } else if (simTime > 0) {
                                    tput = (double) placeCompletions[placeListIdx][k] / simTime;
                                } else {
                                    tput = 0.0;
                                }
                            } else if (simTime > 0) {
                                tput = (double) placeCompletions[placeListIdx][k] / simTime;
                            } else {
                                tput = 0.0;
                            }
                            if (tput > 0) {
                                XN.set(0, k, tput);
                                break;
                            }
                        }
                    }
                }
            }

            // Apply control variates correction if enabled
            if (useControlVariates) {
                applyControlVariateCorrection(QN, UN, RN, TN);
            }

            // Compute arrival rates from throughputs via routing matrix
            Matrix computedAN = jline.api.sn.SnGetArvRFromTput.snGetArvRFromTput(
                    sn, TN, new jline.solvers.AvgHandle());

            LDESResult result = new LDESResult();
            result.QN = QN;
            result.UN = UN;
            result.RN = RN;
            result.TN = TN;
            result.AN = computedAN;
            result.TardN = TardN;
            result.SysTardN = SysTardN;
            result.CN = CN;
            result.XN = XN;
            result.DropRateJoin = dropRateJoin;
            result.sn = sn;

            // Compute FCR metrics if there are any finite capacity regions
            if (numRegions > 0) {
                double simTimeFcr = getActualSimTime();
                for (int f = 0; f < numRegions; f++) {
                    updateRegionTimeWeightedStats(f);
                }
                result.QNfcr = new Matrix(numRegions, numClasses);
                result.UNfcr = new Matrix(numRegions, numClasses);
                result.RNfcr = new Matrix(numRegions, numClasses);
                result.TNfcr = new Matrix(numRegions, numClasses);
                result.ANfcr = new Matrix(numRegions, numClasses);
                result.WNfcr = new Matrix(numRegions, numClasses);
                result.WeightNfcr = new Matrix(numRegions, numClasses);
                result.MemOccNfcr = new Matrix(numRegions, numClasses);
                result.DropRateNfcr = new Matrix(numRegions, numClasses);
                for (int f = 0; f < numRegions; f++) {
                    for (int k = 0; k < numClasses; k++) {
                        double qlen = (simTimeFcr > 0)
                                ? totalRegionJobTime[f][k] / simTimeFcr : 0.0;
                        result.QNfcr.set(f, k, qlen);
                        result.DropRateNfcr.set(f, k, (simTimeFcr > 0)
                                ? (double) droppedByRegion[f][k] / simTimeFcr : 0.0);
                        result.WeightNfcr.set(f, k, (simTimeFcr > 0)
                                ? totalRegionWeightTime[f][k] / simTimeFcr : 0.0);
                        result.MemOccNfcr.set(f, k, (simTimeFcr > 0)
                                ? totalRegionMemTime[f][k] / simTimeFcr : 0.0);
                        result.UNfcr.set(f, k, Double.NaN);
                        double tput = (simTimeFcr > 0)
                                ? (double) regionCompletions[f][k] / simTimeFcr : 0.0;
                        result.TNfcr.set(f, k, tput);
                        double respTime = (tput > 0) ? qlen / tput : 0.0;
                        result.RNfcr.set(f, k, respTime);
                        result.WNfcr.set(f, k, respTime);
                        double arrivalRate;
                        if (regionArrivalCount[f][k] > 1) {
                            double meanInterArrival = regionInterArrivalTimeSum[f][k]
                                    / (regionArrivalCount[f][k] - 1);
                            arrivalRate = (meanInterArrival > 0) ? 1.0 / meanInterArrival : 0.0;
                        } else {
                            arrivalRate = 0.0;
                        }
                        result.ANfcr.set(f, k, arrivalRate);
                    }
                }
            }

            // Variance reduction metadata
            result.varianceReductionMethod = (options.config != null && options.config.variates != null)
                    ? options.config.variates : "none";

            // Convergence data
            result.converged = hasConverged;
            result.stoppingReason = stoppingReason;
            result.convergenceBatches = (queueBatchMeans != null
                    && queueBatchMeans.length > 0
                    && queueBatchMeans[0].length > 0
                    && queueBatchMeans[0][0] != null)
                    ? queueBatchMeans[0][0].size() : 0;

            // Total simulated events
            result.totalSimulatedEvents = totalEventCount;

            // Per-metric sample counts: delegate to tail
            int truncationIdx = mserEnabled ? mserTruncationBatch * effectiveMserBatchSize : 0;
            getLDESResultTail(result, truncationIdx, simTime);

            exportCacheResults();

            populateRewardResult(result);

            return result;
        }

        /**
         * Copies the accumulated Markov reward metrics into the result: the
         * steady-state expectation E[r] = (1/T) integral r(X(t)) dt and the
         * single-run transient step trajectory r(X(t)).
         */
        private void populateRewardResult(LDESResult result) {
            if (exportStateHistogram && stateHistogram != null && !stateHistogram.isEmpty()) {
                int nstates = stateHistogram.size();
                Matrix space = new Matrix(nstates, rewardRowCols);
                Matrix time = new Matrix(nstates, 1);
                int si = 0;
                for (double[] hentry : stateHistogram.values()) {
                    time.set(si, 0, hentry[0]);
                    for (int c = 0; c < rewardRowCols; c++) {
                        space.set(si, c, hentry[c + 1]);
                    }
                    si++;
                }
                result.stateHistogramSpace = space;
                result.stateHistogramTime = time;

                // Integer-state trajectory for the transient reward, downsampled to a
                // bounded number of points so it stays serializable to JSON / a Matrix.
                if (stateTranSeries != null && !stateTranSeries.isEmpty()) {
                    final int MAX_TRAJ = 20000;
                    int total = stateTranSeries.size();
                    int stride = (total + MAX_TRAJ - 1) / MAX_TRAJ;
                    if (stride < 1) {
                        stride = 1;
                    }
                    int npts = (total + stride - 1) / stride;
                    Matrix trajSpace = new Matrix(npts, rewardRowCols);
                    Matrix trajTime = new Matrix(npts, 1);
                    int ti = 0;
                    for (int p = 0; p < total && ti < npts; p += stride) {
                        double[] e = stateTranSeries.get(p);
                        trajTime.set(ti, 0, e[0]);
                        for (int c = 0; c < rewardRowCols; c++) {
                            trajSpace.set(ti, c, e[c + 1]);
                        }
                        ti++;
                    }
                    result.stateTrajectorySpace = trajSpace;
                    result.stateTrajectoryTime = trajTime;
                }
            }
            if (!hasReward || rewardTotalTime <= 0) {
                return;
            }
            result.rewardNames = new ArrayList<String>(rewardNames);
            Map<String, Double> avg = new java.util.LinkedHashMap<String, Double>();
            for (int ri = 0; ri < rewardNames.size(); ri++) {
                avg.put(rewardNames.get(ri), rewardArea[ri] / rewardTotalTime);
            }
            result.avgReward = avg;

            int npts = rewardTranSeries.size();
            Matrix time = new Matrix(npts, 1);
            Map<String, double[]> tran = new java.util.LinkedHashMap<String, double[]>();
            for (int ri = 0; ri < rewardNames.size(); ri++) {
                tran.put(rewardNames.get(ri), new double[npts]);
            }
            for (int p = 0; p < npts; p++) {
                double[] entry = rewardTranSeries.get(p);
                time.set(p, 0, entry[0]);
                for (int ri = 0; ri < rewardNames.size(); ri++) {
                    tran.get(rewardNames.get(ri))[p] = entry[ri + 1];
                }
            }
            result.rewardTime = time;
            result.tranReward = tran;
        }

        /**
         * Export per-cache hit/miss ratios and, for retrieval (delayed-hit) systems, the
         * measured retrieval latency onto the model's Cache nodes. Ratios are indexed by the
         * accessing job class. With a retrieval system the miss ratio counts only
         * fetch-triggering requests (pi_{i,0}); the delayed-hit ratio is the remaining
         * 1 - hit - miss (phi_i), and the expected latency is measured directly from the
         * fetch sojourns and delayed-hit waits (CTMC/SSA report this as NaN).
         */
        private void exportCacheResults() {
            if (cacheStates == null) return;
            for (int cacheNodeIdx : cacheNodes) {
                CacheStateInfo cs = this.cacheStates[cacheNodeIdx];
                if (cs == null) continue;
                jline.lang.nodes.Node node = sn.nodes.get(cacheNodeIdx);
                if (!(node instanceof Cache)) continue;
                Cache cache = (Cache) node;

                Matrix hitProb = new Matrix(1, numClasses);
                Matrix delayedHitProb = new Matrix(1, numClasses);
                Matrix missProb = new Matrix(1, numClasses);
                Matrix latency = new Matrix(1, numClasses);
                hitProb.fill(Double.NaN);
                delayedHitProb.fill(Double.NaN);
                missProb.fill(Double.NaN);
                latency.fill(Double.NaN);

                double meanLatency = Double.NaN;
                if (cs.hasRetrieval && cs.totalDelayedHits != null) {
                    long totalReleased = 0;
                    for (int k = 0; k < numClasses; k++) totalReleased += cs.totalDelayedHits[k];
                    long denom = cs.completedFetches + totalReleased;
                    if (denom > 0) {
                        meanLatency = (cs.totalFetchTime + cs.totalDelayedHitWait) / denom;
                    }
                }

                for (int k = 0; k < numClasses; k++) {
                    long h = cs.totalHits[k];
                    long m = cs.totalMisses[k];
                    long dh = (cs.totalDelayedHits != null) ? cs.totalDelayedHits[k] : 0;
                    long tot = h + m + dh;
                    if (tot > 0) {
                        hitProb.set(0, k, (double) h / tot);
                        delayedHitProb.set(0, k, (double) dh / tot);
                        missProb.set(0, k, (double) m / tot);
                        if (cs.hasRetrieval) {
                            latency.set(0, k, meanLatency);
                        }
                    }
                }

                cache.setResultHitProb(hitProb);
                cache.setResultDelayedHitProb(delayedHitProb);
                cache.setResultMissProb(missProb);
                cache.setResultResidT(latency);

                // Per-list hit fractions [numClasses x h]: share of class-k requests that
                // hit in each list (rows of non-reading classes stay NaN).
                int h = cs.levels.length;
                Matrix hitProbList = new Matrix(numClasses, h);
                hitProbList.fill(Double.NaN);
                for (int k = 0; k < numClasses; k++) {
                    long dh = (cs.totalDelayedHits != null) ? cs.totalDelayedHits[k] : 0;
                    long tot = cs.totalHits[k] + cs.totalMisses[k] + dh;
                    if (tot > 0) {
                        for (int l = 0; l < h; l++) {
                            hitProbList.set(k, l, (double) cs.hitsPerList[k][l] / tot);
                        }
                    }
                }
                cache.setResultHitProbList(hitProbList);

                // Per-item steady-state probabilities [numItems x (h+1)]: column 0 is the
                // miss probability (item not cached), column 1+l the probability the item
                // is in list l, both as time averages over the measurement window.
                accumulateCacheOccupancy(cs);
                double window = cs.lastContentUpdateTime - cs.occupancyStartTime;
                Matrix itemProb = new Matrix(cs.numItems, h + 1);
                itemProb.fill(Double.NaN);
                if (window > 0) {
                    for (int it = 0; it < cs.numItems; it++) {
                        double inCache = 0.0;
                        for (int l = 0; l < h; l++) {
                            double p = cs.itemLevelTime[it][l] / window;
                            itemProb.set(it, l + 1, p);
                            inCache += p;
                        }
                        itemProb.set(it, 0, Math.max(0.0, 1.0 - inCache));
                    }
                }
                cache.setResultItemProb(itemProb);
            }
        }

        LDESResult getTransientLDESResult() {
            return getTransientLDESResultBody();
        }

        // =====================================================================
        // PART 2 TRANSLATION (Kotlin lines 1501-3000).
        //
        // Java port of the tail of the Kotlin {@code init} block (node
        // classification, derived arrays, FC regions) plus the
        // place/transition/cache/fork-join initialization helpers, the
        // distribution generator factories, the cache access helpers and the
        // body of {@link #simulate(double)}.
        //
        // Public surface methods ({@code simulate}, {@code getLDESResult},
        // {@code simulateTransient}, {@code getTransientLDESResult}) remain as
        // stubs above; subsequent chunks rewire the stubs to call into these
        // helpers as more of the source becomes available.
        // =====================================================================

        /**
         * Continuation of the Kotlin {@code init} block — performs node
         * classification, populates derived arrays and finite-capacity-region
         * tables.
         *
         * <p>Mirrors Kotlin lines 1471-1742 (the part not already inlined in
         * the Java constructor).</p>
         */
        private void initFromNetworkStructPart2() {
            // Identify nodes by type
            for (int i = 0; i < numNodes; i++) {
                NodeType nt = sn.nodetype.get(i);
                if (nt == NodeType.Source) {
                    sourceNodes.add(i);
                    sourceStations.add((int) sn.nodeToStation.get(i));
                } else if (nt == NodeType.Queue) {
                    serviceNodes.add(i);
                    serviceStations.add((int) sn.nodeToStation.get(i));
                    isDelayNode.add(Boolean.FALSE);
                } else if (nt == NodeType.Delay) {
                    serviceNodes.add(i);
                    serviceStations.add((int) sn.nodeToStation.get(i));
                    isDelayNode.add(Boolean.TRUE);
                } else if (nt == NodeType.Sink) {
                    sinkNodes.add(i);
                } else if (nt == NodeType.Logger) {
                    loggerNodes.add(i);
                } else if (nt == NodeType.Router) {
                    routerNodes.add(i);
                } else if (nt == NodeType.ClassSwitch) {
                    classSwitchNodes.add(i);
                } else if (nt == NodeType.Fork) {
                    forkNodes.add(i);
                } else if (nt == NodeType.Join) {
                    joinNodes.add(i);
                    joinStations.add((int) sn.nodeToStation.get(i));
                } else if (nt == NodeType.Place) {
                    placeNodes.add(i);
                } else if (nt == NodeType.Transition) {
                    transitionNodes.add(i);
                } else if (nt == NodeType.Cache) {
                    cacheNodes.add(i);
                }
                // ignore other node types
            }

            // Classify classes as open or closed based on sn.njobs.
            // Open: njobs == Inf, Closed: njobs == finite population.
            this.isOpenClass = new boolean[numClasses];
            this.isClosedClass = new boolean[numClasses];
            this.closedClassPopulation = new int[numClasses];
            this.referenceStation = new int[numClasses];
            for (int k = 0; k < numClasses; k++) {
                double n = sn.njobs.get(k);
                this.isOpenClass[k] = Double.isInfinite(n);
                this.isClosedClass[k] = !Double.isInfinite(n);
                this.closedClassPopulation[k] = this.isClosedClass[k] ? (int) n : 0;
                this.referenceStation[k] = (int) sn.refstat.get(k);
            }

            // Signal class detection
            this.isSignalClass = new boolean[numClasses];
            this.isNegativeSignal = new boolean[numClasses];
            for (int k = 0; k < numClasses; k++) {
                this.isSignalClass[k] = sn.issignal != null && sn.issignal.get(k, 0) > 0;
                this.isNegativeSignal[k] = this.isSignalClass[k]
                        && sn.signaltype != null
                        && sn.signaltype.get(k) == SignalType.NEGATIVE;
            }
            this.hasNegativeSignals = false;
            for (int k = 0; k < numClasses; k++) {
                if (this.isNegativeSignal[k]) { this.hasNegativeSignals = true; break; }
            }

            // Batch removal configuration
            this.signalRemovalDist = new DiscreteDistribution[numClasses];
            this.signalRemovalPolicy = new RemovalPolicy[numClasses];
            this.isCatastropheSignal = new boolean[numClasses];
            for (int k = 0; k < numClasses; k++) {
                if (sn.signalremdist != null && k < sn.signalremdist.size()) {
                    this.signalRemovalDist[k] = sn.signalremdist.get(k);
                } else {
                    this.signalRemovalDist[k] = null;
                }
                if (sn.signalrempolicy != null && k < sn.signalrempolicy.size()) {
                    this.signalRemovalPolicy[k] = sn.signalrempolicy.get(k);
                } else {
                    this.signalRemovalPolicy[k] = null;
                }
                this.isCatastropheSignal[k] = sn.iscatastrophe != null
                        && sn.iscatastrophe.get(k, 0) > 0;
            }
            this.hasCatastropheSignals = false;
            for (int k = 0; k < numClasses; k++) {
                if (this.isCatastropheSignal[k]) { this.hasCatastropheSignals = true; break; }
            }
            this.hasRemovalSignals = this.hasNegativeSignals || this.hasCatastropheSignals;

            // REPLY signal detection
            this.isReplySignal = new boolean[numClasses];
            for (int k = 0; k < numClasses; k++) {
                this.isReplySignal[k] = this.isSignalClass[k]
                        && sn.signaltype != null
                        && sn.signaltype.get(k) == SignalType.REPLY;
            }
            this.hasReplySignals = false;
            for (int k = 0; k < numClasses; k++) {
                if (this.isReplySignal[k]) { this.hasReplySignals = true; break; }
            }

            // Synchronous call reply class mapping
            this.synchCallReplyClass = new int[numClasses];
            if (sn.syncreply != null) {
                for (int k = 0; k < numClasses; k++) {
                    this.synchCallReplyClass[k] = (int) sn.syncreply.get(k, 0);
                }
            } else {
                for (int k = 0; k < numClasses; k++) {
                    this.synchCallReplyClass[k] = -1;
                }
            }

            // Spawn-on-completion mapping (LQN phase-2 continuations)
            this.spawnClassOf = new int[numClasses];
            for (int k = 0; k < numClasses; k++) {
                this.spawnClassOf[k] = (sn.classspawn != null && k < sn.classspawn.getNumRows())
                        ? (int) sn.classspawn.get(k, 0) : -1;
            }

            boolean hasOpenClasses = false;
            for (int k = 0; k < numClasses; k++) {
                if (this.isOpenClass[k]) { hasOpenClasses = true; break; }
            }
            boolean hasClosedClasses = false;
            for (int k = 0; k < numClasses; k++) {
                if (this.isClosedClass[k]) { hasClosedClasses = true; break; }
            }
            boolean hasPetriNet = !placeNodes.isEmpty() || !transitionNodes.isEmpty();
            boolean hasCache = !cacheNodes.isEmpty();

            // Validate network structure based on class types
            if (serviceNodes.isEmpty() && !hasPetriNet && !hasCache) {
                throw new RuntimeException(
                        "solver_ssj: Network must have Queue, Delay, Place, Transition, or Cache nodes");
            }
            if (hasOpenClasses && (sourceNodes.isEmpty() || sinkNodes.isEmpty())) {
                throw new RuntimeException(
                        "solver_ssj: Network with open classes must have Source and Sink nodes");
            }

            this.numServiceNodes = serviceNodes.size();
            this.numSources = sourceNodes.size();

            // Extract arrival rates
            this.lambdas = new double[numSources][numClasses];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    double rate = sn.rates.get(sourceStations.get(srcIdx), k);
                    if (rate <= 0 || Double.isNaN(rate) || Double.isInfinite(rate)) {
                        this.lambdas[srcIdx][k] = 0.0;
                    } else {
                        this.lambdas[srcIdx][k] = rate;
                    }
                }
            }

            // Extract service rates (with fallback to heterogeneous rates)
            this.mus = new double[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int stationIdx = serviceStations.get(svcIdx);
                jline.lang.nodes.Station station = sn.stations.get(stationIdx);
                for (int k = 0; k < numClasses; k++) {
                    double rate = sn.rates.get(stationIdx, k);
                    jline.lang.nodeparam.ServiceNodeParam snpRate = sn.getServiceParam(station);
                    if ((rate <= 0 || Double.isNaN(rate) || Double.isInfinite(rate))
                            && snpRate != null && snpRate.heterorates != null) {
                        Map<Integer, Map<Integer, Double>> heteroRates = snpRate.heterorates;
                        if (heteroRates != null) {
                            for (Map.Entry<Integer, Map<Integer, Double>> e : heteroRates.entrySet()) {
                                Map<Integer, Double> classRates = e.getValue();
                                Double heteroRate = classRates.get(k);
                                if (heteroRate != null && heteroRate > 0
                                        && heteroRate < Double.MAX_VALUE) {
                                    rate = heteroRate;
                                    break;
                                }
                            }
                        }
                    }
                    if (rate <= 0 || Double.isNaN(rate) || Double.isInfinite(rate)) {
                        this.mus[svcIdx][k] = Double.MAX_VALUE;
                    } else {
                        this.mus[svcIdx][k] = rate;
                    }
                }
            }

            // Extract number of servers (with support for heterogeneous server types)
            this.numServers = new int[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                if (this.isDelayNode.get(svcIdx)) {
                    this.numServers[svcIdx] = Integer.MAX_VALUE;
                } else {
                    int stationIdx = serviceStations.get(svcIdx);
                    int nservers = (int) sn.nservers.get(stationIdx);
                    jline.lang.nodes.Station station = sn.stations.get(stationIdx);
                    jline.lang.nodeparam.ServiceNodeParam snpServers = sn.getServiceParam(station);
                    if (nservers <= 1 && snpServers != null && snpServers.serverspertype != null) {
                        Matrix serversMatrix = snpServers.serverspertype;
                        if (serversMatrix != null && serversMatrix.length() > 0) {
                            int total = 0;
                            for (int t = 0; t < serversMatrix.length(); t++) {
                                total += (int) serversMatrix.get(t);
                            }
                            if (total > 0) {
                                nservers = total;
                            }
                        }
                    }
                    this.numServers[svcIdx] = nservers;
                }
            }

            // Extract buffer capacities
            this.bufferCapacities = new int[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int capValue = (int) sn.cap.get(serviceStations.get(svcIdx));
                if (capValue <= 0 || capValue == Integer.MAX_VALUE) {
                    this.bufferCapacities[svcIdx] = Integer.MAX_VALUE;
                } else {
                    this.bufferCapacities[svcIdx] = capValue;
                }
            }

            // Extract class capacities
            this.classCapacities = new int[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    double capValue = sn.classcap.get(serviceStations.get(svcIdx), k);
                    if (Double.isNaN(capValue) || capValue <= 0
                            || capValue >= Integer.MAX_VALUE) {
                        this.classCapacities[svcIdx][k] = Integer.MAX_VALUE;
                    } else {
                        this.classCapacities[svcIdx][k] = (int) capValue;
                    }
                }
            }

            // Extract scheduling strategies for each service node
            this.schedStrategies = new SchedStrategy[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                jline.lang.nodes.Station station = sn.stations.get(serviceStations.get(svcIdx));
                SchedStrategy s = (sn.sched != null) ? sn.sched.get(station) : null;
                this.schedStrategies[svcIdx] = (s != null) ? s : SchedStrategy.FCFS;
            }

            // Extract station-level drop rules (for BAS blocking).
            this.stationDropRule = new int[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int stationIdx = serviceStations.get(svcIdx);
                for (int k = 0; k < numClasses; k++) {
                    if (sn.droprule != null) {
                        jline.lang.nodes.Station station = sn.stations.get(stationIdx);
                        jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                        Map<jline.lang.JobClass, DropStrategy> stationMap = sn.droprule.get(station);
                        DropStrategy dropStrategy = (stationMap != null)
                                ? stationMap.get(jobClass) : null;
                        int dropStrategyID = (dropStrategy != null)
                                ? dropStrategy.getID() : DropStrategy.Drop.getID();
                        this.stationDropRule[svcIdx][k] = dropStrategyID;
                    } else {
                        this.stationDropRule[svcIdx][k] = DropStrategy.Drop.getID();
                    }
                }
            }

            // Initialize FC Regions
            this.fcRegionClassMax.fill((double) Integer.MAX_VALUE);
            for (int f = 0; f < numRegions; f++) {
                Matrix regionMatrix = (sn.region != null) ? sn.region.get(f) : null;
                if (regionMatrix == null) continue;

                // Read per-class drop rules
                if (sn.regionrule != null && f < sn.regionrule.getNumRows()) {
                    for (int r = 0; r < numClasses; r++) {
                        if (r < sn.regionrule.getNumCols()) {
                            int dropRuleId = (int) sn.regionrule.get(f, r);
                            this.fcRegionDropRule[f][r] = (dropRuleId == DropStrategy.Drop.getID());
                        }
                    }
                }

                Matrix regionMemMatrix = (sn.regionmaxmem != null) ? sn.regionmaxmem.get(f) : null;
                Matrix regionMemberMask = (sn.regionmembers != null && f < sn.regionmembers.size())
                        ? sn.regionmembers.get(f) : null;
                for (int i = 0; i < numStations; i++) {
                    int globalMax = (int) regionMatrix.get(i, numClasses);
                    double globalMaxMem = (regionMemMatrix != null) ? regionMemMatrix.get(i, 0) : -1.0;
                    // A station is a member of region f if the explicit member
                    // mask says so. The mask exists because membership cannot be
                    // recovered from the caps: -1 means "unbounded", which is
                    // indistinguishable from "not a member", so a region
                    // constrained only by linear admission constraints would
                    // read as empty and be silently ignored. Without a mask
                    // (older structs), fall back to cap-derived membership: a
                    // global job cap, a global memory budget, or any per-class
                    // cap on the row.
                    boolean isMember;
                    if (regionMemberMask != null && i < regionMemberMask.getNumRows()) {
                        isMember = regionMemberMask.get(i, 0) != 0;
                    } else {
                        boolean hasClassCap = false;
                        for (int r = 0; r < numClasses; r++) {
                            if (regionMatrix.get(i, r) >= 0) {
                                hasClassCap = true;
                                break;
                            }
                        }
                        isMember = globalMax >= 0 || globalMaxMem >= 0 || hasClassCap;
                    }
                    if (isMember) {
                        if (this.fcRegionIndices.get(i) < 0) {
                            this.fcRegionIndices.set(i, f);
                        } else if (this.fcRegionIndices.get(i) != f) {
                            // Overlapping finite capacity regions are not
                            // supported (each station binds to one region,
                            // first wins); warn instead of silently ignoring
                            // the later region's constraints.
                            System.err.println("LDES warning: station " + i
                                    + " belongs to finite capacity regions "
                                    + this.fcRegionIndices.get(i) + " and " + f
                                    + "; overlapping regions are unsupported and region "
                                    + f + "'s constraints are ignored at this station.");
                        }
                        if (globalMax >= 0 && globalMax < this.fcRegionGlobalMax.get(f)) {
                            this.fcRegionGlobalMax.set(f, globalMax);
                        }
                        if (globalMaxMem >= 0
                                && (this.fcRegionGlobalMaxMem[f] < 0
                                    || globalMaxMem < this.fcRegionGlobalMaxMem[f])) {
                            this.fcRegionGlobalMaxMem[f] = globalMaxMem;
                        }
                        for (int r = 0; r < numClasses; r++) {
                            int classMax = (int) regionMatrix.get(i, r);
                            if (classMax >= 0
                                    && classMax < this.fcRegionClassMax.get(f, r)) {
                                this.fcRegionClassMax.set(f, r, (double) classMax);
                            }
                        }
                    }
                }
            }

            // Initialize linear admission constraints for FC Regions
            if (sn.regionlincon != null) {
                for (int f = 0; f < numRegions; f++) {
                    MatrixCell linCon = sn.regionlincon.get(f);
                    if (linCon == null) continue;
                    Matrix matA = linCon.get(0);
                    Matrix matB = linCon.get(1);
                    if (matA == null || matB == null) continue;
                    int numConstraints = matA.getNumRows();
                    int numCols = matA.getNumCols();
                    // Pad to numClasses columns: classes appended to the model
                    // after the region was created (e.g. by the fork-join
                    // transformation) carry no coefficient and stay
                    // unconstrained, and the admission check indexes by class.
                    int padCols = Math.max(numCols, numClasses);
                    double[][] aRows = new double[numConstraints][padCols];
                    for (int c = 0; c < numConstraints; c++) {
                        for (int r = 0; r < numCols; r++) {
                            aRows[c][r] = matA.get(c, r);
                        }
                    }
                    double[] bVec = new double[numConstraints];
                    for (int c = 0; c < numConstraints; c++) {
                        bVec[c] = matB.get(c, 0);
                    }
                    this.fcRegionLinConA[f] = aRows;
                    this.fcRegionLinConb[f] = bVec;
                }
            }

            // Initialize Fork/Join, Place/Transition, Cache parameters
            initializeForkJoinParams();
            initializePlaceTransitionParams();
            initializeCacheParams();
        }

        /**
         * Initialize Place and Transition node parameters from network structure.
         *
         * <p>Mirrors Kotlin lines 1747-1941.</p>
         */
        @SuppressWarnings("unchecked")
        private void initializePlaceTransitionParams() {
            int numPlaces = placeNodes.size();
            int numTransitions = transitionNodes.size();

            if (numPlaces == 0 && numTransitions == 0) {
                this.placeTokens = new int[0][];
                this.isQueueingPlace = new boolean[0];
                this.placeNumServers = new int[0];
                this.placeDepository = new int[0][];
                this.placeBusy = new int[0];
                this.currentPlaceInService = new int[0][];
                this.totalPlaceBusyTime = new double[0][];
                this.lastPlaceBusyUpdateTime = new double[0];
                this.placeWaiting = newDequeArray(0);
                this.placeSvcType = new ProcessType[0][];
                this.placeSvcGen = new RandomVariateGen[0][];
                this.placeSvcProc = new MatrixCell[0][];
                this.placeSvcRng = new java.util.Random[0][];
                this.placeSvcMeSampler = new Me_sample.MeSampler[0][];
                this.transitionModes = (List<TransitionModeInfo>[]) new List<?>[0];
                this.transitionFiringGens = new RandomVariateGen[0][];
                this.transitionInService = new int[0][];
                this.totalPlaceTokenTime = new double[0][];
                this.placeCompletions = new int[0][];
                this.lastPlaceUpdateTime = new double[0];
                this.placeTransitTokenTime = new double[0][];
                this.placeTokensInTransit = new int[0][];
                this.lastPlaceTransitUpdateTime = new double[0];
                this.tokensInTransit = new int[0][];
                this.totalTransitTokenTime = new double[0][];
                this.lastTransitUpdateTime = new double[0];
                return;
            }

            // Initialize place tokens from initial state
            if (options.verbose == VerboseLevel.DEBUG) {
                System.out.println("LDES: Initializing place tokens, numPlaces=" + numPlaces
                        + ", numClasses=" + numClasses);
                System.out.println("LDES: sn.state is null: " + (sn.state == null));
                if (sn.state != null) {
                    System.out.println("LDES: sn.state size: " + sn.state.size());
                    for (Map.Entry<jline.lang.nodes.StatefulNode, Matrix> entry : sn.state.entrySet()) {
                        jline.lang.nodes.StatefulNode key = entry.getKey();
                        Matrix value = entry.getValue();
                        System.out.println("LDES:   Key: " + key.getName()
                                + " (" + System.identityHashCode(key) + ")"
                                + ", Value: " + value
                                + ", length=" + value.length());
                    }
                }
                System.out.println("LDES: placeNodes: " + placeNodes);
                for (int i = 0; i < placeNodes.size(); i++) {
                    int nodeIdx = placeNodes.get(i);
                    jline.lang.nodes.Node node = sn.nodes.get(nodeIdx);
                    System.out.println("LDES:   placeNodes[" + i + "] = nodeIdx=" + nodeIdx
                            + ", node=" + node.getName()
                            + " (" + System.identityHashCode(node) + ")");
                }
            }

            this.placeTokens = new int[numPlaces][numClasses];
            for (int placeListIdx = 0; placeListIdx < numPlaces; placeListIdx++) {
                int placeNodeIdx = placeNodes.get(placeListIdx);
                jline.lang.nodes.Node node = sn.nodes.get(placeNodeIdx);
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    Map<jline.lang.nodes.StatefulNode, Matrix> stateMap = sn.state;
                    if (stateMap != null && node instanceof jline.lang.nodes.StatefulNode) {
                        Matrix state = stateMap.get(node);
                        if (options.verbose == VerboseLevel.DEBUG) {
                            System.out.println("LDES: Looking up state for node "
                                    + node.getName()
                                    + " (" + System.identityHashCode(node) + "): found="
                                    + (state != null));
                        }
                        if (state != null && state.length() > classIdx) {
                            int tokens = (int) state.get(classIdx);
                            if (options.verbose == VerboseLevel.DEBUG) {
                                System.out.println("LDES:   class " + classIdx + ": "
                                        + tokens + " tokens");
                            }
                            this.placeTokens[placeListIdx][classIdx] = tokens;
                        } else {
                            if (options.verbose == VerboseLevel.DEBUG) {
                                System.out.println("LDES:   class " + classIdx
                                        + ": state null or too short, returning 0");
                            }
                            this.placeTokens[placeListIdx][classIdx] = 0;
                        }
                    } else {
                        if (options.verbose == VerboseLevel.DEBUG) {
                            System.out.println("LDES:   stateMap null or node not StatefulNode, returning 0");
                        }
                        this.placeTokens[placeListIdx][classIdx] = 0;
                    }
                }
            }

            // Initialize queueing-place (QPN embedded queue) state
            this.isQueueingPlace = new boolean[numPlaces];
            this.placeNumServers = new int[numPlaces];
            this.placeDepository = new int[numPlaces][numClasses];
            this.placeBusy = new int[numPlaces];
            this.currentPlaceInService = new int[numPlaces][numClasses];
            this.totalPlaceBusyTime = new double[numPlaces][numClasses];
            this.lastPlaceBusyUpdateTime = new double[numPlaces];
            this.placeWaiting = newDequeArray(numPlaces);
            this.placeSvcType = new ProcessType[numPlaces][numClasses];
            this.placeSvcGen = new RandomVariateGen[numPlaces][numClasses];
            this.placeSvcProc = new MatrixCell[numPlaces][numClasses];
            this.placeSvcRng = new java.util.Random[numPlaces][numClasses];
            this.placeSvcMeSampler = new Me_sample.MeSampler[numPlaces][numClasses];
            for (int placeListIdx = 0; placeListIdx < numPlaces; placeListIdx++) {
                this.placeWaiting[placeListIdx] = new java.util.ArrayDeque<Integer>();
                int placeNodeIdx = placeNodes.get(placeListIdx);
                int placeStationIdx = (int) sn.nodeToStation.get(placeNodeIdx);
                jline.lang.nodes.Node node = sn.nodes.get(placeNodeIdx);
                boolean queueing = (node instanceof jline.lang.nodes.Place)
                        && ((jline.lang.nodes.Place) node).isQueueing();
                this.isQueueingPlace[placeListIdx] = queueing;
                if (!queueing) {
                    continue;
                }
                jline.lang.nodes.Station station = sn.stations.get(placeStationIdx);
                // Number of servers of the embedded queue (INF -> Integer.MAX_VALUE)
                double ns = (sn.nservers != null) ? sn.nservers.get(placeStationIdx, 0) : 1.0;
                this.placeNumServers[placeListIdx] =
                        (Double.isInfinite(ns) || ns >= Integer.MAX_VALUE) ? Integer.MAX_VALUE : (int) ns;
                // Only FCFS-family and INF embedded scheduling are supported by the algorithm.
                SchedStrategy pss = (sn.sched != null) ? sn.sched.get(station) : SchedStrategy.FCFS;
                if (pss != null && pss != SchedStrategy.FCFS && pss != SchedStrategy.INF
                        && pss != SchedStrategy.SIRO && pss != SchedStrategy.LCFS) {
                    throw new RuntimeException(
                            "LDES: queueing place '" + node.getName() + "' uses scheduling strategy "
                            + pss + ", which the embedded-queue algorithm does not yet support "
                            + "(supported: FCFS, LCFS, SIRO, INF).");
                }
                // Build per-class embedded-queue service samplers from the network struct.
                for (int k = 0; k < numClasses; k++) {
                    jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                    Map<jline.lang.JobClass, ProcessType> pidmap =
                            (sn.procid != null) ? sn.procid.get(station) : null;
                    ProcessType pt = (pidmap != null) ? pidmap.get(jobClass) : null;
                    this.placeSvcType[placeListIdx][k] = (pt != null) ? pt : ProcessType.DISABLED;
                    Map<jline.lang.JobClass, MatrixCell> procmap =
                            (sn.proc != null) ? sn.proc.get(station) : null;
                    MatrixCell proc = (procmap != null) ? procmap.get(jobClass) : null;
                    this.placeSvcProc[placeListIdx][k] = proc;
                    MRG32k3a stream = new MRG32k3a();
                    if (seed > 0) {
                        long offset = ((long) (numSources + numServiceNodes + placeListIdx) * numClasses + k) * 10 + 5000;
                        stream.setSeed(new long[] { seed + offset, seed + offset + 1,
                                seed + offset + 2, seed + offset + 3,
                                seed + offset + 4, seed + offset + 5 });
                        this.placeSvcRng[placeListIdx][k] = new java.util.Random(seed + offset);
                    } else {
                        this.placeSvcRng[placeListIdx][k] = new java.util.Random();
                    }
                    double rate = (sn.rates != null) ? sn.rates.get(placeStationIdx, k) : Double.NaN;
                    if (pt == ProcessType.EXP && rate > 0 && rate < Double.MAX_VALUE) {
                        this.placeSvcGen[placeListIdx][k] =
                                new umontreal.ssj.randvar.ExponentialGen(stream, rate);
                    } else if (pt == ProcessType.DET && rate > 0 && rate < Double.MAX_VALUE) {
                        this.placeSvcGen[placeListIdx][k] = new ConstantGen(stream, 1.0 / rate);
                    } else if (pt == ProcessType.IMMEDIATE) {
                        this.placeSvcGen[placeListIdx][k] = new ConstantGen(stream, 1.0 / GlobalConstants.Immediate);
                    } else {
                        this.placeSvcGen[placeListIdx][k] = null; // renewal PH sampled from proc
                    }
                }
                // Enqueue any initial tokens for service (they must be served before
                // becoming available in the depository). placeTokens already holds the marking.
                for (int k = 0; k < numClasses; k++) {
                    for (int t = 0; t < placeTokens[placeListIdx][k]; t++) {
                        this.placeWaiting[placeListIdx].addLast(Integer.valueOf(k));
                    }
                }
            }

            // Initialize Place statistics tracking
            this.totalPlaceTokenTime = new double[numPlaces][numClasses];
            this.placeCompletions = new int[numPlaces][numClasses];
            this.lastPlaceUpdateTime = new double[numPlaces];
            this.placeTransitTokenTime = new double[numPlaces][numClasses];
            this.placeTokensInTransit = new int[numPlaces][numClasses];
            this.lastPlaceTransitUpdateTime = new double[numPlaces];

            // Initialize transition parameters
            this.transitionModes = (List<TransitionModeInfo>[]) new List<?>[numTransitions];
            for (int transListIdx = 0; transListIdx < numTransitions; transListIdx++) {
                int transNodeIdx = transitionNodes.get(transListIdx);
                jline.lang.nodes.Node node = sn.nodes.get(transNodeIdx);
                Object paramObj = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;

                if (paramObj instanceof TransitionNodeParam) {
                    TransitionNodeParam param = (TransitionNodeParam) paramObj;
                    int nmodes = param.nmodes;
                    List<TransitionModeInfo> modes = new ArrayList<TransitionModeInfo>(nmodes);
                    for (int modeIdx = 0; modeIdx < nmodes; modeIdx++) {
                        String modeName;
                        if (param.modenames != null && modeIdx < param.modenames.size()) {
                            modeName = param.modenames.get(modeIdx);
                        } else {
                            modeName = "Mode" + modeIdx;
                        }
                        TimingStrategy timing;
                        if (param.timing != null && modeIdx < param.timing.size()) {
                            timing = param.timing.get(modeIdx);
                            if (timing == null) timing = TimingStrategy.TIMED;
                        } else {
                            timing = TimingStrategy.TIMED;
                        }
                        int priority = (param.firingprio != null)
                                ? (int) param.firingprio.get(modeIdx) : 0;
                        double weight = (param.fireweight != null)
                                ? param.fireweight.get(modeIdx) : 1.0;
                        int servers = (param.nmodeservers != null)
                                ? (int) param.nmodeservers.get(modeIdx) : Integer.MAX_VALUE;

                        // Enabling conditions for this mode
                        int[][] enabling = new int[numPlaces][numClasses];
                        Matrix enabMatrix = (param.enabling != null
                                && modeIdx < param.enabling.size())
                                ? param.enabling.get(modeIdx) : null;
                        for (int pIdx = 0; pIdx < numPlaces; pIdx++) {
                            int placeNodeIdx = placeNodes.get(pIdx);
                            for (int cIdx = 0; cIdx < numClasses; cIdx++) {
                                if (enabMatrix != null
                                        && placeNodeIdx < enabMatrix.getNumRows()
                                        && cIdx < enabMatrix.getNumCols()) {
                                    enabling[pIdx][cIdx] = (int) enabMatrix.get(placeNodeIdx, cIdx);
                                } else {
                                    enabling[pIdx][cIdx] = 0;
                                }
                            }
                        }

                        // Inhibiting conditions for this mode
                        int[][] inhibiting = new int[numPlaces][numClasses];
                        Matrix inhibMatrix = (param.inhibiting != null
                                && modeIdx < param.inhibiting.size())
                                ? param.inhibiting.get(modeIdx) : null;
                        for (int pIdx = 0; pIdx < numPlaces; pIdx++) {
                            int placeNodeIdx = placeNodes.get(pIdx);
                            for (int cIdx = 0; cIdx < numClasses; cIdx++) {
                                if (inhibMatrix != null
                                        && placeNodeIdx < inhibMatrix.getNumRows()
                                        && cIdx < inhibMatrix.getNumCols()) {
                                    double value = inhibMatrix.get(placeNodeIdx, cIdx);
                                    if (Double.isInfinite(value) || Double.isNaN(value)
                                            || value <= 0) {
                                        inhibiting[pIdx][cIdx] = Integer.MAX_VALUE;
                                    } else {
                                        inhibiting[pIdx][cIdx] = (int) value;
                                    }
                                } else {
                                    inhibiting[pIdx][cIdx] = Integer.MAX_VALUE;
                                }
                            }
                        }

                        // Firing outcomes for this mode
                        int[][] firing = new int[numNodes][numClasses];
                        Matrix fireMatrix = (param.firing != null
                                && modeIdx < param.firing.size())
                                ? param.firing.get(modeIdx) : null;
                        for (int nIdx = 0; nIdx < numNodes; nIdx++) {
                            for (int cIdx = 0; cIdx < numClasses; cIdx++) {
                                if (fireMatrix != null
                                        && nIdx < fireMatrix.getNumRows()
                                        && cIdx < fireMatrix.getNumCols()) {
                                    firing[nIdx][cIdx] = (int) fireMatrix.get(nIdx, cIdx);
                                } else {
                                    firing[nIdx][cIdx] = 0;
                                }
                            }
                        }

                        SerializableFunction<Matrix, Double> firingDep =
                                (param.firingdep != null && modeIdx < param.firingdep.size())
                                        ? param.firingdep.get(modeIdx) : null;
                        modes.add(new TransitionModeInfo(modeIdx, modeName, timing,
                                priority, weight, servers, enabling, inhibiting, firing, firingDep));
                    }
                    this.transitionModes[transListIdx] = modes;
                } else {
                    this.transitionModes[transListIdx] = new ArrayList<TransitionModeInfo>();
                }
            }

            // Initialize firing distribution generators
            this.transitionFiringGens = new RandomVariateGen[numTransitions][];
            for (int transListIdx = 0; transListIdx < numTransitions; transListIdx++) {
                List<TransitionModeInfo> modes = this.transitionModes[transListIdx];
                this.transitionFiringGens[transListIdx] = new RandomVariateGen[modes.size()];
                int transNodeIdx = transitionNodes.get(transListIdx);
                jline.lang.nodes.Node node = sn.nodes.get(transNodeIdx);
                Object paramObj = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
                TransitionNodeParam param = (paramObj instanceof TransitionNodeParam)
                        ? (TransitionNodeParam) paramObj : null;

                for (int modeIdx = 0; modeIdx < modes.size(); modeIdx++) {
                    if (param != null
                            && modes.get(modeIdx).timingStrategy == TimingStrategy.TIMED) {
                        // Find Mode key with index == modeIdx + 1 (1-based).
                        jline.lang.Mode modeKey = null;
                        if (param.firingprocid != null) {
                            for (jline.lang.Mode m : param.firingprocid.keySet()) {
                                if (m.getIndex() == modeIdx + 1) {
                                    modeKey = m;
                                    break;
                                }
                            }
                        }
                        ProcessType procType = (modeKey != null && param.firingprocid != null)
                                ? param.firingprocid.get(modeKey) : null;
                        MatrixCell procCell = (modeKey != null && param.firingproc != null)
                                ? param.firingproc.get(modeKey) : null;

                        if (procType != null && procCell != null) {
                            this.transitionFiringGens[transListIdx][modeIdx] =
                                    createDistributionGenerator(procType, procCell,
                                            transListIdx, modeIdx);
                        } else {
                            // Default: exponential with rate 1
                            this.transitionFiringGens[transListIdx][modeIdx] =
                                    new umontreal.ssj.randvar.ExponentialGen(new MRG32k3a(), 1.0);
                        }
                    } else {
                        this.transitionFiringGens[transListIdx][modeIdx] = null; // Immediate
                    }
                }
            }

            // Initialize in-service counters
            this.transitionInService = new int[numTransitions][];
            for (int transListIdx = 0; transListIdx < numTransitions; transListIdx++) {
                this.transitionInService[transListIdx] =
                        new int[this.transitionModes[transListIdx].size()];
            }

            // Initialize tokens in transit tracking
            this.tokensInTransit = new int[numTransitions][numClasses];
            this.totalTransitTokenTime = new double[numTransitions][numClasses];
            this.lastTransitUpdateTime = new double[numTransitions];
        }

        /**
         * Create a distribution generator for transition firing time.
         *
         * <p>Mirrors Kotlin lines 1946-1997.</p>
         */
        private RandomVariateGen createDistributionGenerator(ProcessType procType,
                                                             MatrixCell procCell,
                                                             int transListIdx,
                                                             int modeIdx) {
            MRG32k3a stream = new MRG32k3a();

            switch (procType) {
                case EXP: {
                    Matrix d1 = procCell.get(1);
                    double rate = (d1 != null && d1.length() > 0) ? d1.get(0, 0) : 1.0;
                    return new umontreal.ssj.randvar.ExponentialGen(stream, rate);
                }
                case ERLANG: {
                    Matrix d0 = procCell.get(0);
                    int k = (d0 != null && d0.getNumRows() > 0) ? d0.getNumRows() : 1;
                    double phaseRate = (d0 != null && d0.length() > 0)
                            ? FastMath.abs(d0.get(0, 0)) : 1.0;
                    return new umontreal.ssj.randvar.ErlangGen(stream, k, phaseRate);
                }
                case HYPEREXP: {
                    Matrix d0 = procCell.get(0);
                    Matrix d1 = procCell.get(1);
                    double lambda1 = (d0 != null && d0.length() > 0) ? -d0.get(0, 0) : 1.0;
                    double lambda2 = (d0 != null && d0.getNumRows() > 1)
                            ? -d0.get(1, 1) : 1.0;
                    double p = (d1 != null && d1.length() > 0 && lambda1 > 0)
                            ? d1.get(0, 0) / lambda1 : 0.5;
                    return new HyperExponentialDistGen(stream,
                            new double[]{p, 1.0 - p},
                            new double[]{lambda1, lambda2});
                }
                case DMAP: {
                    return new DmapSampleGen(stream, procCell,
                            new Random((long) stream.nextInt(0, Integer.MAX_VALUE)));
                }
                case PH:
                case APH:
                case COXIAN:
                case COX2:
                case MAP:
                case MMPP2: {
                    return new MapSampleGen(stream, procCell,
                            new Random((long) stream.nextInt(0, Integer.MAX_VALUE)));
                }
                case ME: {
                    return new MeSampleGen(stream, procCell,
                            new Random((long) stream.nextInt(0, Integer.MAX_VALUE)));
                }
                case RAP: {
                    return new RapSampleGen(stream, procCell,
                            new Random((long) stream.nextInt(0, Integer.MAX_VALUE)));
                }
                case PARETO: {
                    // firingproc stores Pareto's native {shape, scale} (see
                    // Network.refreshStruct); mean = shape*scale/(shape-1).
                    if (procCell != null && procCell.size() >= 2) {
                        double shape = procCell.get(0).get(0, 0);
                        double scale = procCell.get(1).get(0, 0);
                        if (shape > 1.0 && scale > 0) {
                            return new umontreal.ssj.randvar.ParetoGen(stream, shape, scale);
                        }
                    }
                    return new umontreal.ssj.randvar.ExponentialGen(stream, 1.0);
                }
                case DET:
                case UNIFORM:
                case GAMMA:
                case WEIBULL:
                case LOGNORMAL: {
                    // These firing distributions are stored as {mean, scv} because a
                    // transition is not a station and has no sn.rates/sn.scv entry.
                    double mean = (procCell != null && procCell.size() >= 1)
                            ? procCell.get(0).get(0, 0) : 1.0;
                    double scv = (procCell != null && procCell.size() >= 2)
                            ? procCell.get(1).get(0, 0) : 1.0;
                    RandomVariateGen gen = firingGenFromMeanScv(procType, mean, scv, stream);
                    if (gen != null) return gen;
                    double rate = (mean > 0) ? 1.0 / mean : 1.0;
                    return new umontreal.ssj.randvar.ExponentialGen(stream, rate);
                }
                case IMMEDIATE:
                    return null;
                default:
                    return new umontreal.ssj.randvar.ExponentialGen(stream, 1.0);
            }
        }

        /**
         * Build a firing-time generator for a non-Markovian distribution specified by
         * its mean and squared coefficient of variation (the representation stored in
         * {@code firingproc} for DET/UNIFORM/GAMMA/WEIBULL/LOGNORMAL). The moment-to-
         * parameter mappings mirror the service-time generator so that a transition
         * firing time matches the identically distributed station service time.
         */
        private RandomVariateGen firingGenFromMeanScv(ProcessType procType, double mean,
                                                      double scv, MRG32k3a stream) {
            if (procType == ProcessType.DET) {
                return new ConstantGen(stream, mean);
            }
            if (!(mean > 0) || !(scv > 0)) {
                return null;
            }
            if (procType == ProcessType.GAMMA) {
                double shape = 1.0 / scv;
                double lambda = 1.0 / (mean * scv);
                return new umontreal.ssj.randvar.GammaGen(stream, shape, lambda);
            } else if (procType == ProcessType.WEIBULL) {
                double c = FastMath.sqrt(scv);
                double r = FastMath.pow(c, -1.086);
                double alpha = mean
                        / org.apache.commons.math3.special.Gamma.gamma(1.0 + 1.0 / r);
                if (r > 0 && alpha > 0) {
                    return new umontreal.ssj.randvar.WeibullGen(stream, r, 1.0 / alpha, 0.0);
                }
            } else if (procType == ProcessType.LOGNORMAL) {
                double c2plus1 = scv + 1.0;
                double mu = FastMath.log(mean / FastMath.sqrt(c2plus1));
                double sigma = FastMath.sqrt(FastMath.log(c2plus1));
                if (sigma > 0) {
                    return new umontreal.ssj.randvar.LognormalGen(stream, mu, sigma);
                }
            } else if (procType == ProcessType.UNIFORM) {
                // U[a,b] with mean (a+b)/2 and variance (b-a)^2/12.
                double halfWidth = mean * FastMath.sqrt(3.0 * scv);
                double minVal = Math.max(0.0, mean - halfWidth);
                double maxVal = mean + halfWidth;
                if (maxVal > minVal) {
                    return new umontreal.ssj.randvar.UniformGen(stream, minVal, maxVal);
                }
            }
            return null;
        }

        /**
         * Custom generator for hyper-exponential distribution (mixture of exponentials).
         *
         * <p>Mirrors Kotlin {@code HyperExponentialDistGen} (lines 2002-2018).</p>
         */
        private static final class HyperExponentialDistGen extends RandomVariateGen {
            private final RandomStream stream;
            private final double[] probs;
            private final double[] rates;

            HyperExponentialDistGen(RandomStream stream, double[] probs, double[] rates) {
                super(stream, null);
                this.stream = stream;
                this.probs = probs;
                this.rates = rates;
            }

            @Override
            public double nextDouble() {
                double u = stream.nextDouble();
                double cumProb = 0.0;
                for (int i = 0; i < probs.length; i++) {
                    cumProb += probs[i];
                    if (u <= cumProb) {
                        return -FastMath.log(stream.nextDouble()) / rates[i];
                    }
                }
                return -FastMath.log(stream.nextDouble()) / rates[rates.length - 1];
            }
        }

        /**
         * Custom generator wrapper for PH/MAP distributions using {@code map_sample}.
         *
         * <p>Mirrors Kotlin {@code MapSampleGen} (lines 2023-2037). Inner (non-static)
         * because the original was an inner class capturing the enclosing instance.</p>
         */
        private final class MapSampleGen extends RandomVariateGen {
            private final MatrixCell proc;
            private final Random rng;

            MapSampleGen(RandomStream stream, MatrixCell proc, Random rng) {
                super(stream, null);
                this.proc = proc;
                this.rng = rng;
            }

            @Override
            public double nextDouble() {
                if (proc.size() >= 2) {
                    Matrix d0 = proc.get(0);
                    Matrix d1 = proc.get(1);
                    double[] samples = jline.api.mam.Map_sample.map_sample(d0, d1, 1L, rng);
                    return samples[0];
                }
                return 0.0;
            }
        }

        /**
         * Recovers the ME initial vector alpha from a process cell {A, -A*e*alpha}.
         *
         * <p>An ME process always stores a rank-one D1, so a failure to recover
         * alpha means the cell is not an ME representation and must not be
         * sampled as one.</p>
         */
        private Matrix meAlphaOf(MatrixCell proc, String role) {
            Matrix alpha = Me_sample.me_alpha(proc);
            if (alpha == null) {
                throw new RuntimeException(
                        "LDES: ME " + role + " process has a D1 that is not the rank-one "
                        + "matrix (-A*e)*alpha implied by an ME representation, so its "
                        + "initial vector cannot be recovered. Declare the process as a "
                        + "RAP if it is a general rational arrival process.");
            }
            return alpha;
        }

        /**
         * Generator wrapper for ME distributions using inverse-CDF sampling.
         *
         * <p>Built on a cached {@link Me_sample.MeSampler} so that the inversion
         * table is constructed once per generator, not once per variate. The
         * uniform draws come from the SSJ stream so that seeding and
         * reproducibility match the other generators.</p>
         */
        private final class MeSampleGen extends RandomVariateGen {
            private final Me_sample.MeSampler sampler;
            private final Random rng;

            MeSampleGen(RandomStream stream, MatrixCell proc, Random rng) {
                super(stream, null);
                this.sampler = new Me_sample.MeSampler(meAlphaOf(proc, "service"), proc.get(0));
                this.rng = rng;
            }

            @Override
            public double nextDouble() {
                return sampler.next(rng);
            }
        }

        /**
         * Generator wrapper for RAP distributions using conditional-vector sampling.
         *
         * <p>The sampler is cached so that the conditional phase vector persists
         * across variates, which is what preserves the RAP autocorrelation (the
         * role {@link Map_sample.MapSampler} plays for a MAP).</p>
         */
        private final class RapSampleGen extends RandomVariateGen {
            private final Rap_sample.RapSampler sampler;
            private final Random rng;

            RapSampleGen(RandomStream stream, MatrixCell proc, Random rng) {
                super(stream, null);
                this.sampler = new Rap_sample.RapSampler(proc.get(0), proc.get(1));
                this.rng = rng;
            }

            @Override
            public double nextDouble() {
                return sampler.next(rng);
            }
        }

        /**
         * Custom generator for DMAP distributions using {@code dmap_sample}.
         *
         * <p>Mirrors Kotlin {@code DmapSampleGen} (lines 2042-2056).</p>
         */
        private final class DmapSampleGen extends RandomVariateGen {
            private final MatrixCell proc;
            private final Random rng;

            DmapSampleGen(RandomStream stream, MatrixCell proc, Random rng) {
                super(stream, null);
                this.proc = proc;
                this.rng = rng;
            }

            @Override
            public double nextDouble() {
                if (proc.size() >= 2) {
                    Matrix d0 = proc.get(0);
                    Matrix d1 = proc.get(1);
                    double[] samples = jline.api.mam.Dmap_sample.dmap_sample(d0, d1, 1, rng);
                    return samples[0];
                }
                return 0.0;
            }
        }

        /**
         * Create a random variate generator for setup/delayoff distributions.
         *
         * <p>Supports Exponential, Erlang, HyperExp, and PH/APH/Coxian distributions.
         * For PH-type distributions, uses the process representation with map_sample.
         * Mirrors Kotlin lines 2063-2110.</p>
         */
        private RandomVariateGen createSetupDelayoffGen(jline.lang.processes.Distribution dist,
                                                        MRG32k3a stream) {
            if (dist instanceof jline.lang.processes.Erlang) {
                MatrixCell proc = ((jline.lang.processes.Markovian) dist).getProcess();
                Matrix d0 = proc.get(0);
                int phases = d0.getNumRows();
                double phaseRate = FastMath.abs(d0.get(0, 0));
                return new umontreal.ssj.randvar.ErlangGen(stream, phases, phaseRate);
            }
            if (dist instanceof jline.lang.processes.HyperExp) {
                MatrixCell proc = ((jline.lang.processes.Markovian) dist).getProcess();
                Matrix d0 = proc.get(0);
                Matrix d1 = proc.get(1);
                int nPhases = d0.getNumRows();
                double[] rates = new double[nPhases];
                for (int i = 0; i < nPhases; i++) {
                    rates[i] = -d0.get(i, i);
                }
                double[] probs = new double[nPhases];
                if (nPhases == 2 && d1 != null && rates[0] > 0) {
                    probs[0] = d1.get(0, 0) / rates[0];
                    probs[1] = 1.0 - probs[0];
                } else {
                    double sum = 0.0;
                    for (int i = 0; i < nPhases; i++) {
                        probs[i] = (rates[i] > 0) ? d1.get(0, i) / rates[i] : 0.0;
                        sum += probs[i];
                    }
                    if (sum > 0) {
                        for (int i = 0; i < nPhases; i++) probs[i] /= sum;
                    }
                }
                return new HyperExponentialDistGen(stream, probs, rates);
            }
            if (dist instanceof jline.lang.processes.Markovian) {
                // PH, APH, Coxian, Cox2: use map_sample via MapSampleGen. ME and RAP
                // are not phase-type, so they get their own samplers.
                MatrixCell proc = ((jline.lang.processes.Markovian) dist).getProcess();
                if (proc != null && proc.size() >= 2 && !proc.get(0).hasNaN()) {
                    Random rng = new Random();
                    if (seed > 0) {
                        rng.setSeed(seed + (long) stream.hashCode());
                    }
                    if (dist instanceof jline.lang.processes.ME) {
                        return new MeSampleGen(stream, proc, rng);
                    }
                    if (dist instanceof jline.lang.processes.RAP) {
                        return new RapSampleGen(stream, proc, rng);
                    }
                    return new MapSampleGen(stream, proc, rng);
                }
            }
            // Default: Exponential with rate = 1/mean
            return new umontreal.ssj.randvar.ExponentialGen(stream, 1.0 / dist.getMean());
        }

        /**
         * Initialize Fork and Join node parameters from network structure.
         *
         * <p>Mirrors Kotlin lines 2115-2188.</p>
         */
        private void initializeForkJoinParams() {
            int numForks = forkNodes.size();
            int numJoins = joinNodes.size();

            // Initialize Fork parameters
            this.forkFanOut = new int[numForks];
            for (int i = 0; i < numForks; i++) this.forkFanOut[i] = 1;

            // Extract fanOut from ForkNodeParam for each fork node
            for (int forkListIdx = 0; forkListIdx < numForks; forkListIdx++) {
                int forkNodeIdx = forkNodes.get(forkListIdx);
                jline.lang.nodes.Node node = sn.nodes.get(forkNodeIdx);
                Object paramObj = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
                if (paramObj instanceof jline.lang.nodeparam.ForkNodeParam) {
                    jline.lang.nodeparam.ForkNodeParam param =
                            (jline.lang.nodeparam.ForkNodeParam) paramObj;
                    double fanOut = param.fanOut;
                    this.forkFanOut[forkListIdx] = (Double.isNaN(fanOut) || fanOut <= 0)
                            ? 1 : (int) fanOut;
                }
            }

            // Initialize Join parameters
            this.joinToForkMap = new int[numJoins];
            for (int i = 0; i < numJoins; i++) this.joinToForkMap[i] = -1;
            this.forkToJoinMap = new int[numForks];
            for (int i = 0; i < numForks; i++) this.forkToJoinMap[i] = -1;

            this.joinStrategies = new JoinStrategy[numJoins][numClasses];
            for (int j = 0; j < numJoins; j++) {
                for (int k = 0; k < numClasses; k++) {
                    this.joinStrategies[j][k] = JoinStrategy.STD;
                }
            }
            this.joinRequired = new int[numJoins][numClasses];
            for (int j = 0; j < numJoins; j++) {
                for (int k = 0; k < numClasses; k++) {
                    this.joinRequired[j][k] = -1; // -1 means all tasks
                }
            }

            // Build Fork-Join pairing from sn.fj matrix
            Matrix fjMatrix = sn.fj;
            if (fjMatrix != null) {
                for (int joinListIdx = 0; joinListIdx < numJoins; joinListIdx++) {
                    int joinNodeIdx = joinNodes.get(joinListIdx);
                    for (int forkListIdx = 0; forkListIdx < numForks; forkListIdx++) {
                        int forkNodeIdx = forkNodes.get(forkListIdx);
                        if (fjMatrix.get(forkNodeIdx, joinNodeIdx) > 0) {
                            this.joinToForkMap[joinListIdx] = forkListIdx;
                            this.forkToJoinMap[forkListIdx] = joinListIdx;
                            break;
                        }
                    }
                }
            }

            // Extract join strategy and required counts from JoinNodeParam
            for (int joinListIdx = 0; joinListIdx < numJoins; joinListIdx++) {
                int joinNodeIdx = joinNodes.get(joinListIdx);
                jline.lang.nodes.Node node = sn.nodes.get(joinNodeIdx);
                Object paramObj = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
                if (paramObj instanceof jline.lang.nodeparam.JoinNodeParam) {
                    jline.lang.nodeparam.JoinNodeParam param =
                            (jline.lang.nodeparam.JoinNodeParam) paramObj;
                    if (param.joinStrategy != null) {
                        for (Map.Entry<jline.lang.JobClass, JoinStrategy> e
                                : param.joinStrategy.entrySet()) {
                            int classIdx = sn.jobclasses.indexOf(e.getKey());
                            if (classIdx >= 0 && classIdx < numClasses) {
                                this.joinStrategies[joinListIdx][classIdx] = e.getValue();
                            }
                        }
                    }
                    if (param.fanIn != null) {
                        for (Map.Entry<jline.lang.JobClass, Double> e : param.fanIn.entrySet()) {
                            int classIdx = sn.jobclasses.indexOf(e.getKey());
                            if (classIdx >= 0 && classIdx < numClasses) {
                                double required = e.getValue();
                                this.joinRequired[joinListIdx][classIdx] =
                                        (required < 0) ? -1 : (int) required;
                            }
                        }
                    }
                }
            }

            // Initialize Join statistics arrays
            this.totalJoinQueueTime = new double[numJoins][numClasses];
            this.joinCompletions = new int[numJoins][numClasses];
            this.lastJoinUpdateTime = new double[numJoins][numClasses];
            this.currentJoinQueueLength = new int[numJoins][numClasses];
            this.joinResponseTimeTally = new Tally[numJoins][numClasses];
            for (int joinListIdx = 0; joinListIdx < numJoins; joinListIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    this.joinResponseTimeTally[joinListIdx][k] =
                            new Tally("Join response time J" + joinListIdx + " C" + k);
                }
            }
            this.arrivedAtJoin = new int[numJoins][numClasses];
            this.droppedByJoin = new int[numJoins][numClasses];
        }

        /**
         * Initialize Cache node parameters from network structure.
         *
         * <p>Sets up cache levels and replacement strategies. Access samplers are
         * deferred until runtime. Mirrors Kotlin lines 2194-2262.</p>
         */
        @SuppressWarnings("unchecked")
        private void initializeCacheParams() {
            int numCaches = cacheNodes.size();
            if (numCaches == 0) {
                this.cacheStates = new CacheStateInfo[numNodes];
                return;
            }

            this.cacheStates = new CacheStateInfo[numNodes];

            for (int cacheListIdx = 0; cacheListIdx < numCaches; cacheListIdx++) {
                int cacheNodeIdx = cacheNodes.get(cacheListIdx);
                jline.lang.nodes.Node node = sn.nodes.get(cacheNodeIdx);

                if (node instanceof Cache) {
                    Cache cache = (Cache) node;
                    int numItems = cache.getNumberOfItems();
                    // The list geometry and the replacement policy must be read from the
                    // same source as accost below. Network.refreshLocalVars rewrites CLIMB
                    // into its equivalent FIFO unit-capacity-list form, remapping itemcap,
                    // replacestrat and accost together; taking the capacity or the policy
                    // off the raw node while accost carries the rewritten geometry would
                    // mix two geometries and simulate a different cache than the model.
                    CacheNodeParam cnpGeom = null;
                    Object npGeomObj = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
                    if (npGeomObj instanceof CacheNodeParam) {
                        cnpGeom = (CacheNodeParam) npGeomObj;
                    }
                    Matrix itemLevelCap = (cnpGeom != null && cnpGeom.itemcap != null)
                            ? cnpGeom.itemcap : cache.getItemLevelCap();
                    // nLevels is the number of non-empty lists, as in Cache's constructor.
                    int numLevels = (cnpGeom != null && cnpGeom.itemcap != null)
                            ? itemLevelCap.getNonZeroLength() : cache.getnLevels();

                    int[] levelCapacities = new int[numLevels];
                    for (int level = 0; level < numLevels; level++) {
                        if (itemLevelCap != null && level < itemLevelCap.length()) {
                            levelCapacities[level] = (int) itemLevelCap.get(level);
                        } else {
                            levelCapacities[level] = 1;
                        }
                    }

                    ReplacementStrategy replacementStrategy = (cnpGeom != null && cnpGeom.replacestrat != null)
                            ? cnpGeom.replacestrat : cache.getReplacementStrategy();
                    if (replacementStrategy == null) {
                        replacementStrategy = ReplacementStrategy.LRU;
                    }

                    LinkedList<Integer>[] levels = (LinkedList<Integer>[]) new LinkedList<?>[numLevels];
                    for (int level = 0; level < numLevels; level++) {
                        levels[level] = new LinkedList<Integer>();
                    }

                    int[] hitClassArr = new int[numClasses];
                    int[] missClassArr = new int[numClasses];
                    for (int k = 0; k < numClasses; k++) {
                        hitClassArr[k] = k;
                        missClassArr[k] = k;
                    }

                    jline.lang.sections.CacheClassSwitcher cacheServer = cache.getCacheServer();
                    if (cacheServer != null) {
                        Matrix hitClassMatrix = cacheServer.hitClass;
                        Matrix missClassMatrix = cacheServer.missClass;
                        for (int k = 0; k < numClasses; k++) {
                            if (hitClassMatrix != null && k < hitClassMatrix.getNumCols()) {
                                int hitIdx = (int) hitClassMatrix.get(0, k);
                                if (hitIdx >= 0) hitClassArr[k] = hitIdx;
                            }
                            if (missClassMatrix != null && k < missClassMatrix.getNumCols()) {
                                int missIdx = (int) missClassMatrix.get(0, k);
                                if (missIdx >= 0) missClassArr[k] = missIdx;
                            }
                        }
                    }

                    long[] totalHits = new long[numClasses];
                    long[] totalMisses = new long[numClasses];
                    this.cacheStates[cacheNodeIdx] = new CacheStateInfo(
                            cacheNodeIdx,
                            numItems,
                            levelCapacities,
                            replacementStrategy,
                            levels,
                            hitClassArr,
                            missClassArr,
                            totalHits,
                            totalMisses);

                    // --- retrieval system (delayed hits): bind per-item retrieval classes ---
                    CacheStateInfo cs = this.cacheStates[cacheNodeIdx];
                    cs.qlru = cache.getAdmissionProb();
                    Object npObj = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
                    if (npObj instanceof CacheNodeParam) {
                        CacheNodeParam cnp = (CacheNodeParam) npObj;
                        if (cnp.retrievalSystemCapacity > 0 && cnp.retrievalClasses != null) {
                            Matrix rc = cnp.retrievalClasses;   // [items x classes] -> retrieval class index, -1 if none
                            cs.hasRetrieval = true;
                            cs.retrievalClass = new int[numItems][numClasses];
                            cs.retrievalClassToItem = new int[numClasses];
                            for (int c = 0; c < numClasses; c++) cs.retrievalClassToItem[c] = -1;
                            for (int it = 0; it < numItems; it++) {
                                for (int c = 0; c < numClasses; c++) {
                                    int rci = -1;
                                    if (it < rc.getNumRows() && c < rc.getNumCols()) {
                                        rci = (int) rc.get(it, c);
                                    }
                                    cs.retrievalClass[it][c] = rci;
                                    if (rci >= 0 && rci < numClasses) {
                                        cs.retrievalClassToItem[rci] = it;
                                    }
                                }
                            }
                            cs.inFlight = new boolean[numItems];
                            cs.heldDelayedHits = (List<HeldRequest>[]) new List<?>[numItems];
                            cs.totalDelayedHits = new long[numClasses];
                            cs.fetchStartTime = new double[numItems];
                        }
                    }

                    // Access-cost (accost) matrices [class][item] -> (h+1)x(h+1): row 0 drives
                    // miss insertion (column 0 = do-not-cache, column 1+l = insert into list l);
                    // row 1+i drives hits in list i (target list sampled over columns 1+i..1+h).
                    // Network.sanitize populates Cache.accessProb (default: super-diagonal =
                    // promote one list); accostFor falls back to that default when absent.
                    if (npObj instanceof CacheNodeParam && ((CacheNodeParam) npObj).accost != null) {
                        cs.accost = ((CacheNodeParam) npObj).accost;
                    } else {
                        cs.accost = cache.accessProb;
                    }

                    // Initial cache contents. sn.state carries [class counts | cache contents
                    // (totalCacheCapacity slots, 1-based item per position, 0 = empty) | per-item
                    // retrieval bitmap]; positions are list-major with the head of each list first.
                    // Without a usable state, fall back to the MATLAB initDefault convention of a
                    // warm cache pre-loaded with items 1..totalCacheCapacity.
                    int totalCap = 0;
                    for (int level = 0; level < numLevels; level++) totalCap += levelCapacities[level];
                    boolean seeded = false;
                    Matrix st = (sn.state != null && node instanceof jline.lang.nodes.StatefulNode)
                            ? sn.state.get((jline.lang.nodes.StatefulNode) node) : null;
                    if (st != null && !st.isEmpty()) {
                        int bitmapWidth = cs.hasRetrieval ? numItems : 0;
                        int offset = st.getNumCols() - totalCap - bitmapWidth;
                        if (offset >= 0) {
                            boolean valid = true;
                            boolean[] seen = new boolean[numItems + 1];
                            for (int p = 0; p < totalCap && valid; p++) {
                                int v = (int) st.get(0, offset + p);
                                if (v < 0 || v > numItems || (v > 0 && seen[v])) valid = false;
                                if (v > 0) seen[v] = true;
                            }
                            if (valid) {
                                for (int level = 0, p = 0; level < numLevels; level++) {
                                    for (int j = 0; j < levelCapacities[level]; j++, p++) {
                                        int v = (int) st.get(0, offset + p);
                                        if (v > 0) levels[level].addLast(v - 1);
                                    }
                                }
                                seeded = true;
                            }
                        }
                    }
                    if (!seeded) {
                        int itm = 0;
                        int fill = Math.min(totalCap, numItems);
                        for (int level = 0; level < numLevels && itm < fill; level++) {
                            for (int j = 0; j < levelCapacities[level] && itm < fill; j++) {
                                levels[level].addLast(itm);
                                itm++;
                            }
                        }
                    }
                }
            }
        }

        /**
         * Initialize cache access samplers. Must be called after routingRng is initialized.
         *
         * <p>Mirrors Kotlin lines 2267-2279.</p>
         */
        private void initializeCacheSamplers() {
            for (int cacheNodeIdx : cacheNodes) {
                CacheStateInfo cacheState = this.cacheStates[cacheNodeIdx];
                if (cacheState == null) continue;
                if (cacheState.accessSamplers == null) {
                    jline.lang.nodes.Node node = sn.nodes.get(cacheNodeIdx);
                    if (node instanceof Cache) {
                        Cache cache = (Cache) node;
                        cacheState.accessSamplers = createAccessSamplers(cache, cacheState.numItems);
                    }
                }
            }
        }

        /**
         * Create per-class item samplers from the per-class popularity distributions set
         * on the cache via {@code setRead()} (mirroring MATLAB {@code pread}). A class
         * without a discrete popularity gets a null sampler; if such a class actually
         * reads the cache, {@link #processCacheAccess} raises an error, matching the
         * MATLAB requirement that every read class defines a discrete popularity.
         */
        private DiscreteAccessSampler[] createAccessSamplers(Cache cache, int numItems) {
            DiscreteAccessSampler[] samplers = new DiscreteAccessSampler[numClasses];
            CacheNodeParam cnp = null;
            Object npObj = (sn.nodeparam != null) ? sn.nodeparam.get(cache) : null;
            if (npObj instanceof CacheNodeParam) {
                cnp = (CacheNodeParam) npObj;
            }
            for (int r = 0; r < numClasses; r++) {
                double[] probs = null;
                // Preferred source: sn nodeparam pread (per-class PMF over items 1..n)
                if (cnp != null && cnp.pread != null && cnp.pread.get(r) != null) {
                    List<Double> pr = cnp.pread.get(r);
                    if (pr.size() >= numItems) {
                        probs = new double[numItems];
                        double sum = 0.0;
                        for (int i = 0; i < numItems; i++) {
                            double v = pr.get(i);
                            if (Double.isNaN(v) || v < 0) {
                                probs = null;
                                break;
                            }
                            probs[i] = v;
                            sum += v;
                        }
                        if (probs != null && sum <= 0) probs = null;
                    }
                }
                // Fallback: evaluate the class's popularity distribution directly
                if (probs == null) {
                    jline.lang.processes.Distribution popularityDist = cache.popularityGet(r);
                    if (popularityDist instanceof DiscreteDistribution && !popularityDist.isDisabled()) {
                        List<Double> itemIndices = new ArrayList<Double>(numItems);
                        for (int i = 1; i <= numItems; i++) {
                            itemIndices.add((double) i);
                        }
                        Matrix pmfMatrix = ((DiscreteDistribution) popularityDist).evalPMF(itemIndices);
                        if (pmfMatrix.length() >= numItems) {
                            probs = new double[numItems];
                            for (int i = 0; i < numItems; i++) {
                                probs[i] = pmfMatrix.get(i);
                            }
                        }
                    }
                }
                samplers[r] = (probs != null) ? new DiscreteAccessSampler(routingRng, probs) : null;
            }
            return samplers;
        }

        /**
         * Discrete sampler for cache item accesses.
         *
         * <p>Mirrors Kotlin {@code DiscreteAccessSampler} (lines 2333-2364).</p>
         */
        private static final class DiscreteAccessSampler extends RandomVariateGen {
            private final MRG32k3a rng;
            private final double[] cdf;

            DiscreteAccessSampler(MRG32k3a rng, double[] probs) {
                super(rng, null);
                this.rng = rng;
                this.cdf = new double[probs.length];
                double sum = 0.0;
                for (int i = 0; i < probs.length; i++) {
                    sum += probs[i];
                    this.cdf[i] = sum;
                }
                if (sum > 0) {
                    for (int i = 0; i < this.cdf.length; i++) {
                        this.cdf[i] /= sum;
                    }
                }
            }

            @Override
            public double nextDouble() {
                double u = rng.nextDouble();
                for (int i = 0; i < cdf.length; i++) {
                    if (u <= cdf[i]) return (double) i;
                }
                return (double) (cdf.length - 1);
            }

            int nextItem() {
                return (int) nextDouble();
            }
        }

        /**
         * Process a cache access for a job, returning the new class after hit/miss
         * determination.
         *
         * <p>Mirrors Kotlin lines 2372-2399.</p>
         *
         * @param cacheNodeIdx the cache node index
         * @param currentClass the current job class
         * @return the new class after cache access (hit class or miss class)
         */
        private int processCacheAccess(int cacheNodeIdx, int currentClass, long jobId) {
            CacheStateInfo cacheState = this.cacheStates[cacheNodeIdx];
            if (cacheState == null) return currentClass;

            // Retrieval system (delayed hits): a returning retrieval-class job completes a fetch.
            if (cacheState.hasRetrieval) {
                int fetchedItem = cacheState.retrievalClassToItem[currentClass];
                if (fetchedItem >= 0) {
                    // Fetch complete: the item enters the cache and the in-flight flag clears.
                    // The miss was already counted when the fetch was triggered, so do not
                    // count again here; route the returning request on to its miss class.
                    handleCacheMiss(cacheState, currentClass, fetchedItem);
                    cacheState.inFlight[fetchedItem] = false;
                    cacheState.totalFetchTime += (ssjSim.time() - cacheState.fetchStartTime[fetchedItem]);
                    cacheState.completedFetches++;
                    // Release every request parked while the item was in flight: each reads the
                    // now-cached item and completes instantaneously, counted as a delayed hit.
                    releaseDelayedHits(cacheState, fetchedItem);
                    return cacheState.missClass[currentClass];
                }
            }

            DiscreteAccessSampler sampler = (cacheState.accessSamplers != null
                    && currentClass < cacheState.accessSamplers.length)
                    ? cacheState.accessSamplers[currentClass] : null;
            if (sampler == null) {
                throw new RuntimeException("solver_ssj: class "
                        + sn.jobclasses.get(currentClass).getName() + " reads cache "
                        + sn.nodes.get(cacheNodeIdx).getName()
                        + " but has no discrete popularity distribution (Cache.setRead)");
            }
            int itemIdx = sampler.nextItem();

            int hitLevel = -1;
            for (int level = 0; level < cacheState.levels.length; level++) {
                if (cacheState.levels[level].contains(Integer.valueOf(itemIdx))) {
                    hitLevel = level;
                    break;
                }
            }

            if (hitLevel >= 0) {
                cacheState.totalHits[currentClass]++;
                cacheState.hitsPerList[currentClass][hitLevel]++;
                handleCacheHit(cacheState, currentClass, itemIdx, hitLevel);
                return cacheState.hitClass[currentClass];
            }

            // Miss. With a retrieval system the miss triggers (or joins) a fetch.
            if (cacheState.hasRetrieval) {
                int rClass = (currentClass < cacheState.retrievalClass[itemIdx].length)
                        ? cacheState.retrievalClass[itemIdx][currentClass] : -1;
                if (rClass >= 0) {
                    if (cacheState.inFlight[itemIdx]) {
                        // Item already being fetched: park this request as a delayed hit.
                        if (cacheState.heldDelayedHits[itemIdx] == null) {
                            cacheState.heldDelayedHits[itemIdx] = new ArrayList<HeldRequest>();
                        }
                        cacheState.heldDelayedHits[itemIdx].add(
                                new HeldRequest(currentClass, ssjSim.time(), jobId));
                        return CACHE_HELD;
                    }
                    // Trigger a fetch: count the miss now and switch to the per-item retrieval
                    // class; the existing routing (rtnodes injected by setRetrievalSystem) carries
                    // it through the retrieval queues and back to the cache.
                    cacheState.totalMisses[currentClass]++;
                    cacheState.inFlight[itemIdx] = true;
                    cacheState.fetchStartTime[itemIdx] = ssjSim.time();
                    return rClass;
                }
            }

            // Plain miss (no retrieval system, or no retrieval class for this class).
            cacheState.totalMisses[currentClass]++;
            handleCacheMiss(cacheState, currentClass, itemIdx);
            return cacheState.missClass[currentClass];
        }

        /**
         * Release all requests parked at {@code cacheState} while item {@code itemIdx} was
         * being fetched. The fetch has just completed and the item is now cached, so each
         * held request reads it and completes instantaneously (a delayed hit), routing on
         * via its hit class.
         */
        private void releaseDelayedHits(CacheStateInfo cacheState, int itemIdx) {
            List<HeldRequest> held = cacheState.heldDelayedHits[itemIdx];
            if (held == null || held.isEmpty()) return;
            double now = ssjSim.time();
            for (HeldRequest hr : held) {
                cacheState.totalDelayedHits[hr.accessClass]++;
                cacheState.totalDelayedHitWait += (now - hr.holdTime);
                int hitCls = cacheState.hitClass[hr.accessClass];
                deliverCacheCompletion(cacheState.nodeIdx, hitCls, hr.holdTime, hr.jobId);
            }
            held.clear();
        }

        /**
         * Route a request that has just completed at a cache (a released delayed hit) onward
         * from the cache node with the given (hit) class, preserving its original arrival time
         * for response-time accounting. Mirrors the dispatch in routeSignalToNextDestination.
         */
        private void deliverCacheCompletion(int cacheNode, int classId, double systemArrivalTime, long jobId) {
            RoutingResult rr = selectDestination(cacheNode, classId);
            int destNode = rr.destNode;
            int destClassId = rr.destClassId;
            if (destNode < 0) {
                return;
            }
            if (sinkNodes.contains(Integer.valueOf(destNode))) {
                if (warmupDone) {
                    systemResponseTimeTally[classId].add(ssjSim.time() - systemArrivalTime);
                    systemCompletedCustomers[classId]++;
                }
                checkEventCountStop();
                return;
            }
            if (forkNodes.contains(Integer.valueOf(destNode))) {
                long parentJobId = nextJobId++;
                handleForkArrival(destNode, parentJobId, destClassId, systemArrivalTime);
                return;
            }
            if (placeNodes.contains(Integer.valueOf(destNode))) {
                handlePlaceArrival(destNode, destClassId, systemArrivalTime);
                return;
            }
            int nextQueueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
            if (nextQueueIdx >= 0) {
                long jid = (synchCallReplyClass[destClassId] >= 0) ? nextJobId++ : -1L;
                Customer customer = new Customer(
                        destClassId, classPrio[destClassId],
                        systemArrivalTime, ssjSim.time(),
                        siroRng.nextDouble(),
                        -1.0, jid,
                        ssjSim.time() + classDeadline[destClassId],
                        -1, null);
                arriveAtQueue(nextQueueIdx, customer);
            }
        }

        /**
         * Handle a cache hit under the segmented (multi-list) cache model.
         *
         * <p>Mirrors {@code afterEventCache}: the target list {@code inew} is sampled
         * from row {@code 1+i} of the class/item access-cost matrix over columns
         * {@code 1+i..1+h} (no demotion; default super-diagonal = promote one list).
         * If the item stays ({@code inew == i}), LRU moves it to the head and
         * RR/FIFO/SFIFO leave it in place. On promotion the displaced item of the
         * target list (its tail for LRU/FIFO/SFIFO, a uniformly random occupant for
         * RR) returns to the hit list: to the head for LRU/SFIFO, to the freed slot
         * for RR/FIFO.</p>
         */
        private void handleCacheHit(CacheStateInfo cacheState, int accessClass, int itemIdx, int hitLevel) {
            accumulateCacheOccupancy(cacheState);
            ReplacementStrategy rs = cacheState.replacementStrategy;
            int inew = sampleCacheHitTarget(cacheState, accessClass, itemIdx, hitLevel);
            if (inew <= hitLevel) {
                // Item stays in its list: LRU/HLRU/QLRU refresh recency (move to head),
                // others leave in place.
                if (rs == ReplacementStrategy.LRU || rs == ReplacementStrategy.HLRU
                        || rs == ReplacementStrategy.QLRU) {
                    cacheState.levels[hitLevel].remove(Integer.valueOf(itemIdx));
                    cacheState.levels[hitLevel].addFirst(itemIdx);
                }
                return;
            }
            // Promote from list hitLevel to list inew.
            LinkedList<Integer> lo = cacheState.levels[hitLevel];
            LinkedList<Integer> hi = cacheState.levels[inew];
            int kpos = lo.indexOf(Integer.valueOf(itemIdx));
            lo.remove(Integer.valueOf(itemIdx));
            if (hi.size() < cacheState.levelCapacities[inew]) {
                // Higher list not yet full (warm-up): move up without demotion.
                hi.addFirst(itemIdx);
                return;
            }
            // Higher list full: promote the item and demote one item back down.
            int demoted;
            if (rs == ReplacementStrategy.RR) {
                int r = (int) (routingRng.nextDouble() * hi.size());
                if (r < 0) r = 0;
                if (r > hi.size() - 1) r = hi.size() - 1;
                demoted = hi.remove(r);
                hi.add(r, Integer.valueOf(itemIdx));
            } else {
                // LRU/HLRU/QLRU/FIFO/SFIFO: promote to the head, demote the tail.
                demoted = hi.removeLast();
                hi.addFirst(itemIdx);
            }
            if (rs == ReplacementStrategy.LRU || rs == ReplacementStrategy.SFIFO
                    || rs == ReplacementStrategy.HLRU || rs == ReplacementStrategy.QLRU) {
                lo.addFirst(Integer.valueOf(demoted));   // demoted becomes most-recent in list i
            } else {
                int at = (kpos < 0 || kpos > lo.size()) ? lo.size() : kpos;
                lo.add(at, Integer.valueOf(demoted));    // RR/FIFO: demoted takes the freed slot
            }
        }

        /**
         * Handle a cache miss under the segmented (multi-list) cache model.
         *
         * <p>Mirrors {@code afterEventCache}: the insertion target is sampled from the
         * miss row (row 0) of the class/item access-cost matrix — column 0 rejects the
         * item (pass through without caching), column {@code 1+l} inserts it into list
         * {@code l} (default: the entry list). RR replaces a uniformly random occupant;
         * LRU/FIFO/SFIFO insert at the head and evict the tail when the list is full.</p>
         */
        private void handleCacheMiss(CacheStateInfo cacheState, int accessClass, int itemIdx) {
            accumulateCacheOccupancy(cacheState);
            int target = sampleCacheMissTarget(cacheState, accessClass, itemIdx);
            if (target < 0) {
                return;   // cache reject: item is not cached
            }
            ReplacementStrategy rs = cacheState.replacementStrategy;
            LinkedList<Integer> list = cacheState.levels[target];
            int cap = cacheState.levelCapacities[target];
            if (rs == ReplacementStrategy.QLRU && routingRng.nextDouble() > cacheState.qlru) {
                return;   // q-LRU: not admitted on this miss, pass through uncached
            }
            if (rs == ReplacementStrategy.RR) {
                if (list.size() >= cap && !list.isEmpty()) {
                    int r = (int) (routingRng.nextDouble() * list.size());
                    if (r < 0) r = 0;
                    if (r > list.size() - 1) r = list.size() - 1;
                    list.set(r, Integer.valueOf(itemIdx));   // replace a random occupant
                } else {
                    list.addFirst(itemIdx);
                }
            } else {
                // LRU/HLRU/QLRU/FIFO/SFIFO: insert at head, evict tail when full.
                list.addFirst(itemIdx);
                if (list.size() > cap) {
                    list.removeLast();
                }
            }
        }

        /**
         * Access-cost matrix for the given class and item, or null when unavailable
         * (callers then apply the default super-diagonal behaviour).
         */
        private Matrix accostFor(CacheStateInfo cacheState, int accessClass, int itemIdx) {
            Matrix[][] ac = cacheState.accost;
            if (ac == null || accessClass < 0 || accessClass >= ac.length) return null;
            Matrix[] row = ac[accessClass];
            if (row == null || itemIdx < 0 || itemIdx >= row.length) return null;
            Matrix m = row[itemIdx];
            if (m == null || m.getNumRows() < cacheState.levels.length + 1
                    || m.getNumCols() < cacheState.levels.length + 1) return null;
            return m;
        }

        /**
         * Sample the target list for a hit in list {@code hitLevel}: row {@code 1+hitLevel}
         * of the access-cost matrix, columns {@code 1+hitLevel..1+h} (normalized as in
         * MATLAB afterEventCache). Returns a 0-based list index in {@code [hitLevel, h-1]};
         * defaults to promote-one-list when no matrix is available.
         */
        private int sampleCacheHitTarget(CacheStateInfo cacheState, int accessClass, int itemIdx, int hitLevel) {
            int h = cacheState.levels.length;
            Matrix ac = accostFor(cacheState, accessClass, itemIdx);
            if (ac == null) {
                return Math.min(hitLevel + 1, h - 1);
            }
            double sum = 0.0;
            for (int c = hitLevel + 1; c <= h; c++) {
                sum += ac.get(hitLevel + 1, c);
            }
            if (sum <= 0.0) {
                return Math.min(hitLevel + 1, h - 1);
            }
            double u = routingRng.nextDouble() * sum;
            double acc = 0.0;
            for (int c = hitLevel + 1; c <= h; c++) {
                acc += ac.get(hitLevel + 1, c);
                if (u <= acc) {
                    return c - 1;
                }
            }
            return h - 1;
        }

        /**
         * Sample the insertion list for a missed item from the miss row (row 0) of the
         * access-cost matrix. Returns the 0-based target list, or -1 for cache reject
         * (column 0). Defaults to the entry list (list 0) when no matrix is available.
         */
        private int sampleCacheMissTarget(CacheStateInfo cacheState, int accessClass, int itemIdx) {
            int h = cacheState.levels.length;
            Matrix ac = accostFor(cacheState, accessClass, itemIdx);
            if (ac == null) {
                return 0;
            }
            double sum = 0.0;
            for (int c = 0; c <= h; c++) {
                sum += ac.get(0, c);
            }
            if (sum <= 0.0) {
                return 0;
            }
            double u = routingRng.nextDouble() * sum;
            double acc = 0.0;
            for (int c = 0; c <= h; c++) {
                acc += ac.get(0, c);
                if (u <= acc) {
                    return c - 1;
                }
            }
            return h - 1;
        }

        /**
         * Accumulate the time-integrated per-item, per-list presence up to the current
         * simulation time. Called before every cache content mutation and at export.
         */
        private void accumulateCacheOccupancy(CacheStateInfo cacheState) {
            double now = ssjSim.time();
            double dt = now - cacheState.lastContentUpdateTime;
            if (dt > 0) {
                for (int level = 0; level < cacheState.levels.length; level++) {
                    for (Integer item : cacheState.levels[level]) {
                        cacheState.itemLevelTime[item][level] += dt;
                    }
                }
            }
            cacheState.lastContentUpdateTime = now;
        }

        /**
         * Body of the steady-state {@code simulate} entry point — runs the
         * simulation until the specified number of service-completion events.
         * Uses MSER-5 for automatic transient detection or falls back to a
         * fixed warmup if disabled.
         *
         * <p>Mirrors Kotlin lines 2519-2600. Named {@code simulateBody} (rather
         * than overriding {@link #simulate(double)}) so that the public
         * {@code simulate} stub can be wired in by a later translation chunk.</p>
         *
         * @param maxEventCount maximum number of service-completion events
         */
        private void simulateBody(double maxEventCount) {
            ssjSim.init();

            // Initialize transient detection and CI configuration from LDESOptions
            initializeTransientAndCIConfig();

            // Initialize variance reduction mode from options
            String variatesMode = (options.config != null && options.config.variates != null)
                    ? options.config.variates : "none";
            this.useAntitheticVariates = variatesMode.equals("antithetic")
                    || variatesMode.equals("both");
            this.useControlVariates = variatesMode.equals("control")
                    || variatesMode.equals("both");

            if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                if (this.useAntitheticVariates) {
                    System.out.println("LDES: Antithetic variates enabled");
                }
                if (this.useControlVariates) {
                    System.out.println("LDES: Control variates enabled");
                }
            }

            // Initialize event-count based stopping
            this.maxEvents = (long) maxEventCount;
            this.totalEventCount = 0L;
            this.totalSimEvents = 0L;
            if (options instanceof LDESOptions) {
                this.maxSimEventsLimit = (long) ((LDESOptions) options).maxSimEvents;
                this.maxTimeLimit = ((LDESOptions) options).maxTime;
            }
            // Fall back to the generic wall-clock budget on the base SolverOptions
            // (the wrapper builds a base SolverOptions, not an LDESOptions).
            if (Double.isInfinite(this.maxTimeLimit) && options != null
                    && Double.isFinite(options.timeout) && options.timeout > 0) {
                this.maxTimeLimit = options.timeout;
            }
            this.simStartNanos = System.nanoTime();
            this.warmupEventThreshold = (long) (this.maxEvents * this.effectiveWarmupFraction);

            // Initialize trace logging if DEBUG verbose level
            initializeTracing();

            // Initialize Logger nodes for CSV output
            initializeLoggers();

            // Initialize routing strategies
            initializeRouting();

            // Initialize random generators
            initializeGenerators();

            // Initialize cache samplers (requires routingRng to be initialized)
            initializeCacheSamplers();

            // Initialize queues and statistics
            initializeState();

            // Validate closed class routing (no routes to Sink)
            validateClosedClassRouting();

            // Initialize closed class populations - use init_sol if provided (e.g. derived
            // from Network.initFromMarginal), otherwise place jobs at reference stations
            this.initializing = true;
            if (initSol != null && !initSol.isEmpty()) {
                initFromInitSol();
            } else {
                initClosedClassPopulations();
            }
            this.initializing = false;

            if (this.mserEnabled) {
                // MSER-5 mode: collect observations throughout simulation,
                // determine truncation at the end
                initializeMSEREventBased(this.maxEvents);
                this.warmupDone = true;
                this.warmupEndTime = 0.0;
            } else {
                // Legacy fixed warmup mode — warmup triggered by event count
                this.warmupDone = false;
            }

            // Initialize convergence checking (runs alongside MSER)
            initializeConvergence(this.maxEvents);

            // Schedule initial arrivals from each source for OPEN classes only.
            // Closed classes have fixed population already initialized at reference stations.
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    if (this.isOpenClass[k]) {
                        if (isMarkedNonCarrier(srcIdx, k)) {
                            continue; // driven by the marked group's carrier class
                        }
                        double interarrivalTime = generateInterarrivalTime(srcIdx, k);
                        if (interarrivalTime > 0) {
                            new ExternalArrival(srcIdx, k).schedule(interarrivalTime);
                        }
                    }
                }
            }

            // Start service for any initial tokens sitting in queueing places' embedded queues.
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                if (isQueueingPlace[placeListIdx]) {
                    tryStartPlaceService(placeListIdx);
                }
            }

            // Kick off transition firings from the initial marking. Closed Petri
            // nets are triggered via initClosedClassPopulations, but open nets with
            // initial place tokens (and no closed-class injection) need an explicit
            // start. Idempotent: in-flight counts prevent double-scheduling.
            if (!placeNodes.isEmpty() && !transitionNodes.isEmpty()) {
                checkAndFireTransitions();
            }

            // Run simulation (stops when maxEvents reached via checkEventCountStop())
            ssjSim.start();
        }

        // =====================================================================
        // Forward-reference stubs for methods translated in PART 3+ (Kotlin
        // lines 3000+). The {@link #simulateBody(double)} above calls into
        // these; later chunks will replace each stub with the real
        // implementation.
        // =====================================================================

        // initializeTransientAndCIConfig() — implementation in PART 3 below.

        // initializeRouting() — translated in PART 6 (Kotlin lines 6044-6140).

        // validateClosedClassRouting(), initClosedClassPopulations(),
        // initializeMSEREventBased(long), initializeConvergence(long) —
        // implementations in PART 3 below.

        /**
         * True when class k belongs to a marked (MMAP) group at this source but
         * is not its carrier: such classes receive jobs through the carrier's
         * ExternalArrival events (mark->class routing) and must not schedule
         * their own arrival stream.
         */
        private boolean isMarkedNonCarrier(int sourceIdx, int k) {
            int[] group = markedGroupClass[sourceIdx];
            if (group == null || k == markedCarrierClass[sourceIdx]) {
                return false;
            }
            for (int m = 1; m < group.length; m++) {
                if (group[m] == k) {
                    return true;
                }
            }
            return false;
        }

        /** Generates an interarrival time for source/class. Kotlin line 12410. */
        /**
         * Validates that a sampled interval lies on the slot lattice and returns
         * it with the floating-point drift removed.
         *
         * <p>In continuous mode the value passes through untouched. In slotted
         * mode a value that is not a lattice point is an error: rounding it would
         * silently substitute a different distribution for the modelled one, and
         * the resulting metrics would look plausible while describing a system the
         * user never specified. A zero interval is also rejected, because the
         * engine reads zero as "no further event" and on a lattice the smallest
         * meaningful interval is one slot.
         *
         * @param value the sampled interval
         * @param what  description of the sample, used in the error message
         * @return the value snapped to the lattice
         */
        /**
         * Resolves the zero atom of a counting distribution.
         *
         * <p>Bernoulli, Binomial and Poisson are supported on {0,1,...}, so used
         * as an interval they produce 0 with positive probability. The engine
         * reads a zero interarrival time as "no further event" (see the
         * {@code interarrivalTime > 0} guards at the scheduling sites), so
         * leaving it alone would silently retire the source rather than model a
         * zero-length interval.
         *
         * <p>In continuous time a zero-length interval is exactly what
         * {@code Immediate} denotes, and the engine already represents that as
         * {@code 1/GlobalConstants.Immediate}; reusing the same constant keeps
         * the two consistent. In slotted mode there is no such constant: the
         * smallest interval on the lattice is one slot, so a zero sample is
         * rejected instead, matching {@link #slotSnap}'s refusal to invent a
         * lattice point.
         *
         * @param value    the sampled interval
         * @param procType the process type it came from
         * @param what     description of the sample, used in the error message
         * @return the interval to use
         */
        private double resolveZeroAtom(double value, ProcessType procType, String what) {
            if (value != 0.0) {
                return value;
            }
            if (procType != ProcessType.BERNOULLI && procType != ProcessType.BINOMIAL
                    && procType != ProcessType.POISSON) {
                return value;
            }
            if (slotted) {
                throw new RuntimeException("LDES slotted mode: sampled " + what
                        + " is 0, which is not a point of the slot lattice. " + procType
                        + " is supported on {0,1,...}; a slotted model needs an interval "
                        + "distribution supported on {1,2,...}, such as Geometric.");
            }
            return 1.0 / GlobalConstants.Immediate;
        }

        private double slotSnap(double value, String what) {
            if (!slotted) {
                return value;
            }
            if (value == 0.0) {
                return 0.0;
            }
            double slots = value / slotLength;
            double rounded = FastMath.rint(slots);
            if (rounded < 1.0 || FastMath.abs(slots - rounded) > 1e-9 * FastMath.max(1.0, slots)) {
                throw new RuntimeException("LDES slotted mode: sampled " + what + " is "
                        + value + ", which is not a positive multiple of the slot length "
                        + slotLength + ". Discrete-time simulation requires every "
                        + "interarrival and service time to be lattice valued, e.g. "
                        + "Geometric or Det with an integral number of slots.");
            }
            return rounded * slotLength;
        }

        /**
         * Rejects models whose semantics are not defined on a discrete time scale.
         *
         * <p>Two families are excluded. Rate-schedule processes (NHPP) integrate a
         * cumulative intensity over continuous time, so their inverse transform
         * does not land on the lattice. Processor-sharing disciplines divide the
         * server continuously among the jobs in service, which has no counterpart
         * in a model where the server completes at most one job per slot.
         */
        private void assertSlottedModelIsSupported(NetworkStruct sn) {
            if (sn.procid != null) {
                for (Map<JobClass, ProcessType> byClass : sn.procid.values()) {
                    if (byClass == null) continue;
                    for (ProcessType pt : byClass.values()) {
                        if (pt == ProcessType.NHPP) {
                            throw new RuntimeException("LDES slotted mode does not support NHPP: "
                                    + "a rate schedule is integrated over continuous time and "
                                    + "does not yield lattice-valued intervals.");
                        }
                    }
                }
            }
            if (sn.sched != null) {
                for (Map.Entry<Station, SchedStrategy> e : sn.sched.entrySet()) {
                    SchedStrategy s = e.getValue();
                    if (s == SchedStrategy.PS || s == SchedStrategy.DPS || s == SchedStrategy.GPS
                            || s == SchedStrategy.PSPRIO || s == SchedStrategy.DPSPRIO
                            || s == SchedStrategy.GPSPRIO || s == SchedStrategy.LPS) {
                        throw new RuntimeException("LDES slotted mode does not support "
                                + s + " at station '" + e.getKey().getName() + "': processor "
                                + "sharing divides the server continuously, which has no "
                                + "discrete-time counterpart.");
                    }
                }
            }
        }

        private double generateInterarrivalTime(int sourceIdx, int classId) {
            double interarrivalTime = 0.0;

            // Check for Replayer/trace distribution first
            ProcessType procType = arrivalProcessType[sourceIdx][classId];
            if (procType == ProcessType.REPLAYER) {
                TraceSampler traceSampler = arrivalTraceSamplers[sourceIdx][classId];
                if (traceSampler != null) {
                    interarrivalTime = traceSampler.nextSample();
                }
            } else {
                RandomVariateGen gen = arrivalGens[sourceIdx][classId];
                if (gen != null) {
                    interarrivalTime = gen.nextDouble();
                } else if (procType == ProcessType.BMAP) {
                    // BMAP: sample inter-arrival time and batch size, carrying the
                    // modulating phase across arrivals to preserve autocorrelation.
                    MatrixCell proc = arrivalProc[sourceIdx][classId];
                    Random rng = arrivalRng[sourceIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 3) {
                        Map_sample.BmapSampler s = arrivalBmapSampler[sourceIdx][classId];
                        if (s == null) {
                            s = new Map_sample.BmapSampler(proc);
                            arrivalBmapSampler[sourceIdx][classId] = s;
                        }
                        jline.api.mam.BmapSample sample = s.next(rng);
                        interarrivalTime = sample.getInterarrivalTime();
                        arrivalBatchSize[sourceIdx][classId] = sample.getBatchSize();
                    }
                } else if (procType == ProcessType.MMAP) {
                    // MMAP: stateful marked sampling that carries the modulating
                    // phase across arrivals (autocorrelation preserved) and
                    // records the 1-based mark of the scheduled arrival, consumed
                    // by ExternalArrival to route the job to the class bound to
                    // that mark (Source.setMarkedArrival / sn.markidx).
                    MatrixCell proc = arrivalProc[sourceIdx][classId];
                    Random rng = arrivalRng[sourceIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 2) {
                        jline.api.mam.Mmap_sample.MmapSampler s = arrivalMmapSampler[sourceIdx][classId];
                        if (s == null) {
                            s = new jline.api.mam.Mmap_sample.MmapSampler(proc);
                            arrivalMmapSampler[sourceIdx][classId] = s;
                        }
                        jline.api.mam.Mmap_sample.MarkedSample ms = s.next(rng);
                        interarrivalTime = ms.interarrivalTime;
                        arrivalPendingMark[sourceIdx][classId] = ms.mark;
                    }
                } else if (procType == ProcessType.MAP || procType == ProcessType.MMPP2) {
                    // Correlated arrival process: carry the modulating phase across
                    // arrivals so the inter-arrival autocorrelation is preserved.
                    MatrixCell proc = arrivalProc[sourceIdx][classId];
                    Random rng = arrivalRng[sourceIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 2) {
                        Map_sample.MapSampler s = arrivalMapSampler[sourceIdx][classId];
                        if (s == null) {
                            s = new Map_sample.MapSampler(proc.get(0), proc.get(1));
                            arrivalMapSampler[sourceIdx][classId] = s;
                        }
                        interarrivalTime = s.next(rng);
                    }
                } else if (procType == ProcessType.RAP) {
                    // Rational arrival process: the phase is a conditional vector
                    // rather than a CTMC state, so the autocorrelation is carried by
                    // the RAP vector update instead of a MAP phase walk.
                    MatrixCell proc = arrivalProc[sourceIdx][classId];
                    Random rng = arrivalRng[sourceIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 2) {
                        Rap_sample.RapSampler s = arrivalRapSampler[sourceIdx][classId];
                        if (s == null) {
                            s = new Rap_sample.RapSampler(proc.get(0), proc.get(1));
                            arrivalRapSampler[sourceIdx][classId] = s;
                        }
                        interarrivalTime = s.next(rng);
                    }
                } else if (procType == ProcessType.ME) {
                    // Matrix-exponential renewal arrivals: inverse-CDF sampling from a
                    // cached table (map_sample's CTMC walk is invalid for a genuine ME).
                    MatrixCell proc = arrivalProc[sourceIdx][classId];
                    Random rng = arrivalRng[sourceIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 2) {
                        Me_sample.MeSampler s = arrivalMeSampler[sourceIdx][classId];
                        if (s == null) {
                            s = new Me_sample.MeSampler(meAlphaOf(proc, "arrival"), proc.get(0));
                            arrivalMeSampler[sourceIdx][classId] = s;
                        }
                        interarrivalTime = s.next(rng);
                    }
                } else if (procType == ProcessType.PH || procType == ProcessType.APH
                        || procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN
                        || procType == ProcessType.COX2) {
                    // Renewal phase-type: per-arrival resampling from the marginal is exact.
                    MatrixCell proc = arrivalProc[sourceIdx][classId];
                    Random rng = arrivalRng[sourceIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 2) {
                        Matrix D0 = proc.get(0);
                        Matrix D1 = proc.get(1);
                        double[] samples =
                                jline.api.mam.Map_sample.map_sample(D0, D1, 1L, rng);
                        interarrivalTime = samples[0];
                    }
                } else if (procType == ProcessType.NHPP) {
                    // Rate schedule: inverse transform on the cumulative intensity,
                    // integrating the schedule forward from the current clock.
                    double[][] sched = arrivalSchedule[sourceIdx][classId];
                    if (sched != null) {
                        interarrivalTime = sampleScheduleInterarrival(
                                sched[0], sched[1], sched[2][0] != 0.0,
                                arrivalRng[sourceIdx][classId]);
                    }
                }
            }

            interarrivalTime = resolveZeroAtom(interarrivalTime, procType,
                    "interarrival time at source " + sourceIdx + " class " + classId);
            interarrivalTime = slotSnap(interarrivalTime, "interarrival time at source "
                    + sourceIdx + " class " + classId);

            // Explicit batch-size law. Sampled alongside the interarrival time
            // and left in arrivalBatchSize for ExternalArrival to consume, the
            // same handshake the BMAP branch above uses. A BMAP already carries
            // its own batch sizes, so the two are mutually exclusive by
            // construction: sn.arrivalbatch is only populated from
            // Source.setArrivalBatch.
            DiscreteDistribution batchDist = arrivalBatchDist[sourceIdx][classId];
            if (batchDist != null) {
                double sampled = batchDist.sample(1, arrivalBatchRng[sourceIdx][classId])[0];
                long batch = FastMath.round(sampled);
                if (batch < 1 || FastMath.abs(sampled - batch) > 1e-9) {
                    throw new RuntimeException("LDES: arrival batch size at source " + sourceIdx
                            + " class " + classId + " sampled " + sampled
                            + ", which is not a positive integer. A batch-size law must be "
                            + "supported on {1,2,...}.");
                }
                arrivalBatchSize[sourceIdx][classId] = (int) batch;
            }

            // Track sample for control variates if enabled and after warmup
            if (useControlVariates && warmupDone && interarrivalTime > 0) {
                arrivalSampleSum[sourceIdx][classId] += interarrivalTime;
                arrivalSampleCount[sourceIdx][classId]++;
            }

            return interarrivalTime;
        }

        /**
         * Samples the time to the next event of a process carrying a
         * piecewise-constant rate schedule, by inverse transform on the
         * cumulative intensity.
         *
         * <p>Draws E ~ Exp(1) and returns the T solving
         * int_{now}^{now+T} lambda(u) du = E, walking the schedule forward from
         * the current clock position and consuming E segment by segment. This is
         * exact for a non-homogeneous Poisson process: conditional on no event
         * since the last one, the residual is governed by lambda from the current
         * instant onward, so a holding time drawn under a rate that has since
         * changed is not a sample from this process. Integrating at sampling time
         * makes segment boundaries exact without any rate-change event, since the
         * schedule is a deterministic function of the clock.
         *
         * <p>The NHPP schedule reaches this in breakpoint form, cyclic or not,
         * so there is a single sampler.
         *
         * <p>Returns 0 when no further event can occur: a degenerate schedule, or
         * a non-cyclic horizon that is exhausted. Callers read 0 as "never fires"
         * and schedule no event.
         *
         * @param breakpoints segment boundaries, length n+1
         * @param rates       rate on each segment, length n
         * @param cyclic      whether the schedule repeats with the horizon as period
         */
        private double sampleScheduleInterarrival(double[] breakpoints, double[] rates,
                boolean cyclic, Random rng) {
            if (breakpoints == null || rates == null || rng == null || rates.length == 0
                    || breakpoints.length != rates.length + 1) {
                return 0.0;
            }
            double horizon = breakpoints[rates.length] - breakpoints[0];
            double totalMass = 0.0;
            for (int i = 0; i < rates.length; i++) {
                // A NaN schedule means the rates never reached sn.proc; that used
                // to fall through as a zero holding time, i.e. an infinitely fast
                // station reported as a valid result.
                if (Double.isNaN(rates[i]) || Double.isNaN(breakpoints[i])) {
                    throw new RuntimeException(
                            "LDES: rate schedule contains NaN at segment " + i
                            + ". The schedule did not reach sn.proc; this is a"
                            + " model-compilation defect, not a modelling error.");
                }
                if (rates[i] < 0.0) {
                    throw new RuntimeException(
                            "LDES: rate schedule has a negative rate at segment " + i + ".");
                }
                totalMass += rates[i] * (breakpoints[i + 1] - breakpoints[i]);
            }
            if (!(horizon > 0.0) || !(totalMass > 0.0)) {
                return 0.0;
            }
            // Exp(1) budget of cumulative intensity to consume.
            double residual = -Math.log(1.0 - rng.nextDouble());
            double offset = ssjSim.time() - breakpoints[0];
            if (cyclic) {
                offset = offset % horizon;
                if (offset < 0.0) {
                    offset += horizon;
                }
            } else if (offset >= horizon) {
                return 0.0; // horizon exhausted: the intensity is zero from here on
            } else if (offset < 0.0) {
                offset = 0.0;
            }
            double pos = breakpoints[0] + offset;
            int idx = 0;
            while (idx < rates.length - 1 && pos >= breakpoints[idx + 1]) {
                idx++;
            }
            double elapsed = 0.0;
            // A cyclic schedule consumes totalMass > 0 per cycle, so this
            // terminates; a non-cyclic one exits at the horizon.
            while (true) {
                double remainingInSegment = breakpoints[idx + 1] - pos;
                double massInSegment = rates[idx] * remainingInSegment;
                // The rate guard also keeps a zero-rate segment from dividing 0/0
                // on the measure-zero draw residual == 0.
                if (rates[idx] > 0.0 && massInSegment >= residual) {
                    return elapsed + residual / rates[idx];
                }
                residual -= massInSegment;
                elapsed += remainingInSegment;
                idx++;
                if (idx >= rates.length) {
                    if (!cyclic) {
                        return 0.0;
                    }
                    idx = 0;
                }
                pos = breakpoints[idx];
            }
        }

        /**
         * Unpacks a rate schedule from sn.proc into breakpoint form.
         *
         * <p>NHPP stores {breakpoints, rates, cyclic}, so the engine sees one
         * schedule representation and uses one sampler.
         *
         * @return {breakpoints, rates, {cyclic ? 1 : 0}}, or null if absent
         */
        private double[][] unpackRateSchedule(MatrixCell proc) {
            if (proc == null || proc.get(0) == null || proc.get(1) == null) {
                return null;
            }
            Matrix bpm = proc.get(0);
            Matrix rm = proc.get(1);
            int n = rm.getNumCols();
            double[] breakpoints = new double[n + 1];
            double[] rates = new double[n];
            for (int i = 0; i < n; i++) {
                breakpoints[i] = bpm.get(0, i);
                rates[i] = rm.get(0, i);
            }
            breakpoints[n] = bpm.get(0, n);
            boolean cyclic = proc.get(2) == null || proc.get(2).get(0, 0) != 0.0;
            double[] cyclicFlag = new double[]{cyclic ? 1.0 : 0.0};
            return new double[][]{breakpoints, rates, cyclicFlag};
        }

        // ==================== Operational-time (work-clock) transform ====================
        //
        // A time-varying service rate mu(t) is handled by measuring service in
        // units of cumulative intensity tau = int mu(u) du (operational time)
        // rather than wall clock. In tau-units the process is unit-rate, so the
        // engine's residual/PS/preemption bookkeeping -- which all assumes work
        // is delivered at unit rate -- stays correct once residuals are carried
        // in tau and the two wall<->tau primitives below bridge to the event
        // clock. For a station whose service is NOT an NHPP schedule the rate is
        // constant unity, so both primitives reduce to the identity and every
        // caller is numerically unchanged.

        /** True if (queueIdx, classId) carries a time-varying NHPP service rate. */
        private boolean hasServiceSchedule(int queueIdx, int classId) {
            return serviceProcessType[queueIdx][classId] == ProcessType.NHPP
                    && serviceSchedule != null
                    && serviceSchedule[queueIdx][classId] != null;
        }

        /**
         * Cumulative service intensity int_{b0}^{t} mu(u) du at a station whose
         * piecewise-constant rate schedule is (breakpoints, rates), cyclic with
         * the horizon as period or, when not cyclic, zero outside the horizon.
         * Monotone nondecreasing in t since rates are nonnegative.
         */
        private double scheduleCumulativeIntensity(double[] breakpoints, double[] rates,
                boolean cyclic, double t) {
            int n = rates.length;
            double b0 = breakpoints[0];
            double horizon = breakpoints[n] - b0;
            double totalMass = 0.0;
            for (int i = 0; i < n; i++) {
                totalMass += rates[i] * (breakpoints[i + 1] - breakpoints[i]);
            }
            double x = t - b0;
            double mass;
            double rem;
            if (cyclic) {
                double periods = Math.floor(x / horizon);
                rem = x - periods * horizon; // in [0, horizon)
                mass = periods * totalMass;
            } else if (x <= 0.0) {
                return 0.0;
            } else if (x >= horizon) {
                return totalMass;
            } else {
                rem = x;
                mass = 0.0;
            }
            double segStart = 0.0;
            for (int i = 0; i < n; i++) {
                double segLen = breakpoints[i + 1] - breakpoints[i];
                if (rem <= segStart + segLen) {
                    mass += rates[i] * (rem - segStart);
                    break;
                }
                mass += rates[i] * segLen;
                segStart += segLen;
            }
            return mass;
        }

        /**
         * Operational-time work int_{t0}^{t1} mu(u) du delivered to a job at full
         * server share between two wall-clock instants at (queueIdx, classId).
         * Reduces to t1 - t0 when the station carries no time-varying rate.
         */
        private double serviceWorkBetween(int queueIdx, int classId, double t0, double t1) {
            if (t1 <= t0) {
                return 0.0;
            }
            if (!hasServiceSchedule(queueIdx, classId)) {
                return t1 - t0;
            }
            double[][] sched = serviceSchedule[queueIdx][classId];
            boolean cyclic = sched[2][0] != 0.0;
            return scheduleCumulativeIntensity(sched[0], sched[1], cyclic, t1)
                    - scheduleCumulativeIntensity(sched[0], sched[1], cyclic, t0);
        }

        /**
         * Inverse of {@link #serviceWorkBetween} from a fixed start: the wall
         * clock instant at which `work` units of operational time have been
         * delivered at full share since t0 at (queueIdx, classId). Reduces to
         * t0 + work when the station carries no time-varying rate. Returns
         * {@link Double#POSITIVE_INFINITY} when the (non-cyclic) intensity is
         * exhausted before the work is met, which callers read as "never fires".
         */
        private double serviceWallAfterWork(int queueIdx, int classId, double t0, double work) {
            if (work <= 0.0) {
                return t0;
            }
            if (!hasServiceSchedule(queueIdx, classId)) {
                return t0 + work;
            }
            double[][] sched = serviceSchedule[queueIdx][classId];
            double[] breakpoints = sched[0];
            double[] rates = sched[1];
            boolean cyclic = sched[2][0] != 0.0;
            int n = rates.length;
            double b0 = breakpoints[0];
            double horizon = breakpoints[n] - b0;
            double offset = t0 - b0;
            if (cyclic) {
                offset = offset % horizon;
                if (offset < 0.0) {
                    offset += horizon;
                }
            } else if (offset >= horizon) {
                return Double.POSITIVE_INFINITY;
            } else if (offset < 0.0) {
                offset = 0.0;
            }
            double pos = b0 + offset;
            int idx = 0;
            while (idx < n - 1 && pos >= breakpoints[idx + 1]) {
                idx++;
            }
            double residual = work;
            double elapsed = 0.0;
            while (true) {
                double remainingInSegment = breakpoints[idx + 1] - pos;
                double massInSegment = rates[idx] * remainingInSegment;
                if (rates[idx] > 0.0 && massInSegment >= residual) {
                    return t0 + elapsed + residual / rates[idx];
                }
                residual -= massInSegment;
                elapsed += remainingInSegment;
                idx++;
                if (idx >= n) {
                    if (!cyclic) {
                        return Double.POSITIVE_INFINITY;
                    }
                    idx = 0;
                }
                pos = breakpoints[idx];
            }
        }

        /**
         * Wall-clock delay from now until `work` operational-time units complete
         * at (queueIdx, classId). Identity (returns work) off a time-varying rate.
         */
        private double serviceWallDelay(int queueIdx, int classId, double work) {
            double now = ssjSim.time();
            double t1 = serviceWallAfterWork(queueIdx, classId, now, work);
            return t1 - now;
        }

        /**
         * External arrival from a source node (Kotlin lines 7487-7546).
         * For BMAP (Batch Markovian Arrival Process), multiple customers
         * may arrive simultaneously.
         */
        private final class ExternalArrival extends SimEvent {
            @Override
            protected double slotPhase() {
                return SLOT_PHASE_ARRIVAL;
            }

            final int srcIdx;
            final int classIdx;

            ExternalArrival(int srcIdx, int classIdx) {
                this.srcIdx = srcIdx;
                this.classIdx = classIdx;
            }

            @Override
            public void actions() {
                trackEvent();
                int sourceIdx = this.srcIdx;
                int classId = this.classIdx;
                // Get batch size before generating next arrival (BMAP sets this
                // during generateInterarrivalTime)
                int batchSize = arrivalBatchSize[sourceIdx][classId];
                // Reset batch size to 1 (will be updated by next
                // generateInterarrivalTime if BMAP)
                arrivalBatchSize[sourceIdx][classId] = 1;
                // Get the mark of THIS arrival before generating the next one
                // (MMAP sets the pending mark during generateInterarrivalTime,
                // same pattern as the BMAP batch size)
                int pendingMark = arrivalPendingMark[sourceIdx][classId];
                arrivalPendingMark[sourceIdx][classId] = 0;

                // Schedule next arrival using appropriate generator
                // (EXP or PH/Erlang/HyperExp/BMAP/MMAP)
                double interarrivalTime = generateInterarrivalTime(sourceIdx, classId);
                if (interarrivalTime > 0) {
                    new ExternalArrival(sourceIdx, classId).schedule(interarrivalTime);
                }

                // Marked (MMAP) source: the sampled mark selects the class of
                // the arriving job; the event itself stays on the carrier class.
                int arrivalClassId = classId;
                int[] markGroup = markedGroupClass[sourceIdx];
                if (pendingMark > 0 && markGroup != null
                        && pendingMark < markGroup.length && markGroup[pendingMark] >= 0) {
                    arrivalClassId = markGroup[pendingMark];
                }
                final int classIdArv = arrivalClassId;

                // Process batch: route each customer to destination
                for (int b = 0; b < batchSize; b++) {
                    // Route customer to first destination
                    int sourceNode = sourceNodes.get(sourceIdx).intValue();
                    RoutingResult routingResult = selectDestination(sourceNode, classIdArv);
                    int destNode = routingResult.destNode;
                    int destClassId = routingResult.destClassId;

                    if (destNode >= 0 && sinkNodes.contains(Integer.valueOf(destNode))) {
                        // Customer leaves system immediately
                        // (e.g., Source -> Cache -> Sink)
                        // Track completion in original class for consistency
                        // with service node departures
                        systemResponseTimeTally[classIdArv].add(0.0);
                        systemCompletedCustomers[classIdArv]++;
                        checkEventCountStop();
                    } else if (destNode >= 0 && forkNodes.contains(Integer.valueOf(destNode))) {
                        long parentJobId = nextJobId++;
                        handleForkArrival(destNode, parentJobId, destClassId, ssjSim.time());
                    } else if (destNode >= 0 && placeNodes.contains(Integer.valueOf(destNode))) {
                        handlePlaceArrival(destNode, destClassId, ssjSim.time());
                    } else if (destNode >= 0 && transitionNodes.contains(Integer.valueOf(destNode))) {
                        checkAndFireTransitions();
                    } else if (destNode >= 0) {
                        int queueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
                        if (queueIdx >= 0) {
                            // Check if this is a negative signal or catastrophe
                            if ((hasNegativeSignals && isNegativeSignal[destClassId])
                                    || (hasCatastropheSignals && isCatastropheSignal[destClassId])) {
                                handleNegativeSignalArrival(queueIdx, destClassId, ssjSim.time());
                            } else {
                                // Use destClassId (possibly switched from original classId)
                                long jobId = (synchCallReplyClass[destClassId] >= 0) ? nextJobId++ : -1L;
                                double currentTime = ssjSim.time();
                                Customer customer = new Customer(
                                        destClassId, classPrio[destClassId],
                                        currentTime, currentTime,
                                        siroRng.nextDouble(),
                                        -1.0, jobId,
                                        currentTime + classDeadline[destClassId],
                                        -1, null);
                                arriveAtQueue(queueIdx, customer);
                            }
                        }
                    }
                }
            }
        }

        // =====================================================================
        // PART 3 TRANSLATION (Kotlin lines 3001-4500).
        //
        // This block contains:
        //   * Tail of {@code getLDESResult} (Kotlin lines 3001-3098).
        //   * Helpers {@code setRespTimeSamples},
        //     {@code computeImpatienceStatistics},
        //     {@code applyControlVariateCorrection}.
        //   * Transient-mode entry point {@code simulateTransient}, plus
        //     {@code validateClosedClassRouting}, {@code initClosedClassPopulations},
        //     {@code initFromInitSol}, {@code initializeTransient},
        //     inner {@code TransientSampleEvent}, and
        //     {@code getTransientLDESResult}.
        //   * MSER initialization ({@code initializeMSER},
        //     {@code initializeMSEREventBased}),
        //     {@code initializeTransientAndCIConfig},
        //     {@code initializeConvergence}, event-count tracking
        //     ({@code trackEvent}, {@code checkEventCountStop}),
        //     {@code pushStreamingMetrics}, {@code collectMSERSample},
        //     {@code collectConvergenceSample}, {@code finalizeBatch},
        //     {@code resetBatchAccumulators}, {@code recordResponseTimeForBatch},
        //     {@code finishSimulation}, {@code checkConvergence},
        //     {@code isMetricConverged}, {@code getTCriticalValue},
        //     {@code computeFinalCIMatrices}, {@code computeCIForMetric},
        //     inner {@code MSERSampleEvent}, and the start of
        //     {@code computeMSER5TruncationPoint} (continued in PART 4).
        // =====================================================================

        /**
         * Transient-mode flag: true while {@link #simulateTransient(double)} is
         * executing, false during steady-state simulation.
         */
        private boolean transientMode = false;

        /** Time stamps recorded by the transient sampler (Kotlin: {@code transientTimes}). */
        private List<Double> transientTimes;

        /** Per-(svcIdx,class) transient queue-length time series. */
        private List<Double>[][] transientQueueLengths;

        /** Per-(svcIdx,class) transient utilization time series. */
        private List<Double>[][] transientUtilizations;

        /** Per-(svcIdx,class) transient throughput time series. */
        private List<Double>[][] transientThroughputs;

        /** Per-(svcIdx,class) transient cumulative-completion time series. */
        private List<Integer>[][] transientCompletions;

        /** Sampling interval used by {@link TransientSampleEvent}. */
        private double transientSamplingInterval = 0.0;

        /** Wall clock time of the most recent transient sample. */
        private double lastTransientSampleTime = 0.0;

        /** Cumulative queue time as of the last transient sample, per (svcIdx,class). */
        private double[][] lastTransientQueueTime;

        /** Cumulative busy time as of the last transient sample, per (svcIdx,class). */
        private double[][] lastTransientBusyTime;

        /** Completion count as of the last transient sample, per (svcIdx,class). */
        private int[][] lastTransientCompletions;

        /** Event counter for progress reporting (Kotlin: {@code lastProgressEventCount}). */
        private long lastProgressEventCount = 0L;

        // ---------------------------------------------------------------------
        // getLDESResult tail (Kotlin lines 3001-3098).
        //
        // The full Kotlin {@code getLDESResult} extends from line 2602 to line
        // 3098.  Lines 2602-3000 (the head) are handled by an earlier chunk;
        // this method translates only the tail (lines 3001-3098).
        //
        // Callers in earlier chunks pass in the partially-built {@code result}
        // along with the post-warmup {@code truncationIdx} and the simulation
        // time at which the tail is being constructed.
        // ---------------------------------------------------------------------

        /**
         * Tail of {@code getLDESResult} (Kotlin lines 3001-3098).
         *
         * Populates QNSamples / UNSamples / RNSamples / TNSamples, the various
         * CI / relative-precision matrices, impatience statistics and
         * response-time samples.
         */
        private void getLDESResultTail(LDESResult result, int truncationIdx, double simTime) {
            Matrix QNSamples = new Matrix(numStations, numClasses);
            QNSamples.fill(0.0);
            Matrix UNSamples = new Matrix(numStations, numClasses);
            UNSamples.fill(0.0);
            Matrix RNSamples = new Matrix(numStations, numClasses);
            RNSamples.fill(0.0);
            Matrix TNSamples = new Matrix(numStations, numClasses);
            TNSamples.fill(0.0);
            for (int svcIdx = 0; svcIdx < serviceStations.size(); svcIdx++) {
                int serviceStation = serviceStations.get(svcIdx);
                for (int k = 0; k < numClasses; k++) {
                    // QN samples from queue length observations (post-warmup)
                    if (queueLengthObservations != null
                            && svcIdx < queueLengthObservations.length
                            && k < queueLengthObservations[svcIdx].length) {
                        int totalObs = queueLengthObservations[svcIdx][k].size();
                        QNSamples.set(serviceStation, k, (double) Math.max(0, totalObs - truncationIdx));
                    }
                    // TN/UN samples from throughput observations (post-warmup, minus 1 for differencing)
                    if (throughputObservations != null
                            && svcIdx < throughputObservations.length
                            && k < throughputObservations[svcIdx].length) {
                        int totalObs = throughputObservations[svcIdx][k].size();
                        int postWarmup = Math.max(0, totalObs - truncationIdx);
                        int rates = (postWarmup >= 2) ? postWarmup - 1 : 0;
                        TNSamples.set(serviceStation, k, (double) rates);
                        UNSamples.set(serviceStation, k, (double) rates);
                    }
                    // RN samples from response time tally
                    if (responseTimeTally != null
                            && svcIdx < responseTimeTally.length
                            && k < responseTimeTally[svcIdx].length
                            && responseTimeTally[svcIdx][k] != null) {
                        RNSamples.set(serviceStation, k,
                                (double) responseTimeTally[svcIdx][k].numberObs());
                    }
                }
            }
            result.QNSamples = QNSamples;
            result.UNSamples = UNSamples;
            result.RNSamples = RNSamples;
            result.TNSamples = TNSamples;

            // Transfer CI half-widths
            Matrix QNCI = new Matrix(numStations, numClasses);
            QNCI.fill(0.0);
            Matrix UNCI = new Matrix(numStations, numClasses);
            UNCI.fill(0.0);
            Matrix RNCI = new Matrix(numStations, numClasses);
            RNCI.fill(0.0);
            Matrix TNCI = new Matrix(numStations, numClasses);
            TNCI.fill(0.0);

            // Transfer relative precision
            Matrix QNRelPrec = new Matrix(numStations, numClasses);
            QNRelPrec.fill(0.0);
            Matrix UNRelPrec = new Matrix(numStations, numClasses);
            UNRelPrec.fill(0.0);
            Matrix RNRelPrec = new Matrix(numStations, numClasses);
            RNRelPrec.fill(0.0);
            Matrix TNRelPrec = new Matrix(numStations, numClasses);
            TNRelPrec.fill(0.0);

            if (convergenceEnabled && finalQNCI != null) {
                for (int svcIdx = 0; svcIdx < serviceStations.size(); svcIdx++) {
                    int serviceStation = serviceStations.get(svcIdx);
                    for (int k = 0; k < numClasses; k++) {
                        QNCI.set(serviceStation, k, finalQNCI[svcIdx][k]);
                        UNCI.set(serviceStation, k, finalUNCI[svcIdx][k]);
                        RNCI.set(serviceStation, k, finalRNCI[svcIdx][k]);
                        TNCI.set(serviceStation, k, finalTNCI[svcIdx][k]);

                        QNRelPrec.set(serviceStation, k, finalQNRelPrec[svcIdx][k]);
                        UNRelPrec.set(serviceStation, k, finalUNRelPrec[svcIdx][k]);
                        RNRelPrec.set(serviceStation, k, finalRNRelPrec[svcIdx][k]);
                        TNRelPrec.set(serviceStation, k, finalTNRelPrec[svcIdx][k]);
                    }
                }
            }

            // Compute OBM confidence intervals if enabled
            if (options.confint > 0 && mserEnabled) {
                CIResults ciResults = computeOBMConfidenceIntervals(options.confint);
                result.QNCI = ciResults.QNCI;
                result.UNCI = ciResults.UNCI;
                result.RNCI = ciResults.RNCI;
                result.TNCI = ciResults.TNCI;
                result.ANCI = ciResults.ANCI;
                result.WNCI = ciResults.WNCI;
            } else {
                result.QNCI = QNCI;
                result.UNCI = UNCI;
                result.RNCI = RNCI;
                result.TNCI = TNCI;
            }
            result.QNRelPrec = QNRelPrec;
            result.UNRelPrec = UNRelPrec;
            result.RNRelPrec = RNRelPrec;
            result.TNRelPrec = TNRelPrec;

            // Compute impatience statistics
            computeImpatienceStatistics(result, simTime);

            // Store response time samples for CDF computation
            setRespTimeSamples(result, numStations, numClasses);
        }

        /**
         * Helper to set response time samples array on LDESResult.
         */
        private void setRespTimeSamples(LDESResult result, int numStations, int numClasses) {
            result.initRespTimeSamples(numStations, numClasses);

            for (int ist = 0; ist < numStations; ist++) {
                for (int k = 0; k < numClasses; k++) {
                    int svcIdx = serviceStations.indexOf(ist);
                    if (svcIdx >= 0
                            && responseTimeSamples != null
                            && svcIdx < responseTimeSamples.length
                            && k < responseTimeSamples[svcIdx].length
                            && responseTimeSamples[svcIdx][k] != null) {
                        for (Double sample : responseTimeSamples[svcIdx][k]) {
                            result.respTimeSamples[ist][k].add(sample);
                        }
                    }
                }
            }
        }

        /**
         * Compute impatience (reneging, balking, retrial) statistics and populate LDESResult.
         */
        private void computeImpatienceStatistics(LDESResult result, double simTime) {
            Matrix renegedMatrix = new Matrix(numStations, numClasses);
            Matrix avgRenegingWaitMatrix = new Matrix(numStations, numClasses);
            Matrix renegingRateMatrix = new Matrix(numStations, numClasses);
            Matrix balkedMatrix = new Matrix(numStations, numClasses);
            Matrix balkingProbMatrix = new Matrix(numStations, numClasses);
            Matrix retriedMatrix = new Matrix(numStations, numClasses);
            Matrix retrialDroppedMatrix = new Matrix(numStations, numClasses);
            Matrix avgOrbitMatrix = new Matrix(numStations, numClasses);

            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int ist = serviceStations.get(svcIdx);
                for (int k = 0; k < numClasses; k++) {
                    // Reneging statistics
                    int reneged = renegedCustomers[svcIdx][k];
                    renegedMatrix.set(ist, k, (double) reneged);

                    double avgRenegingWait = (reneged > 0)
                            ? totalRenegingWaitTime[svcIdx][k] / reneged
                            : 0.0;
                    avgRenegingWaitMatrix.set(ist, k, avgRenegingWait);

                    int completed = completedCustomers[svcIdx][k];
                    int dropped = droppedCustomers[svcIdx][k];
                    int total = completed + reneged + dropped;
                    double renegingRate = (total > 0) ? ((double) reneged) / total : 0.0;
                    renegingRateMatrix.set(ist, k, renegingRate);

                    // Balking statistics
                    int balked = balkedCustomers[svcIdx][k];
                    balkedMatrix.set(ist, k, (double) balked);

                    int arrived = arrivedCustomers[svcIdx][k];
                    double balkProb = (arrived > 0) ? ((double) balked) / arrived : 0.0;
                    balkingProbMatrix.set(ist, k, balkProb);

                    // Retrial statistics
                    int retried = retriedCustomers[svcIdx][k];
                    retriedMatrix.set(ist, k, (double) retried);

                    int retrialDrop = maxRetriesExceeded[svcIdx][k];
                    retrialDroppedMatrix.set(ist, k, (double) retrialDrop);

                    // Average orbit size (time-weighted)
                    updateOrbitTimeStats(svcIdx, k);
                    double avgOrbit = (simTime > 0) ? totalOrbitTime[svcIdx][k] / simTime : 0.0;
                    avgOrbitMatrix.set(ist, k, avgOrbit);
                }
            }

            result.renegedCustomers = renegedMatrix;
            result.avgRenegingWaitTime = avgRenegingWaitMatrix;
            result.renegingRate = renegingRateMatrix;
            result.balkedCustomers = balkedMatrix;
            result.balkingProbability = balkingProbMatrix;
            result.retriedCustomers = retriedMatrix;
            result.retrialDropped = retrialDroppedMatrix;
            result.avgOrbitSize = avgOrbitMatrix;
        }

        /**
         * Apply control variates correction to performance metrics.
         * Uses the deviation of sampled service times from theoretical means to adjust estimates.
         */
        private void applyControlVariateCorrection(Matrix QN, Matrix UN, Matrix RN, Matrix TN) {
            for (int svcIdx = 0; svcIdx < serviceStations.size(); svcIdx++) {
                int serviceStation = serviceStations.get(svcIdx);
                for (int k = 0; k < numClasses; k++) {
                    long count = serviceSampleCount[svcIdx][k];
                    double expectedMean = serviceExpectedMean[svcIdx][k];

                    if (count > 0 && expectedMean > 0) {
                        double actualMean = serviceSampleSum[svcIdx][k] / (double) count;

                        if (actualMean > 0 && Double.isFinite(actualMean)) {
                            double correctionFactor = expectedMean / actualMean;

                            // Only apply correction if factor is reasonable (within 50%)
                            if (correctionFactor > 0.5 && correctionFactor < 2.0) {
                                // Correct utilization: U = λ * E[S], so if E[S] biased, U biased
                                double correctedU = UN.get(serviceStation, k) * correctionFactor;
                                if (correctedU < 0.0) {
                                    correctedU = 0.0;
                                } else if (correctedU > 1.0) {
                                    correctedU = 1.0;
                                }
                                UN.set(serviceStation, k, correctedU);

                                // Response time also scales with service time deviation
                                double originalR = RN.get(serviceStation, k);
                                if (originalR > 0) {
                                    // Estimate waiting time contribution
                                    double waitingTime = Math.max(0.0, originalR - actualMean);
                                    // Corrected R = W + E[S]_theoretical
                                    RN.set(serviceStation, k, waitingTime + expectedMean);
                                }
                            }
                        }
                    }
                }
            }
        }

        // ---------------------------------------------------------------------
        // Transient simulation entry point (Kotlin lines 3242-3320).
        // ---------------------------------------------------------------------

        /**
         * Run transient simulation for the specified time horizon.
         * No warmup period - collects metrics from time 0.
         * MSER-5 is disabled during transient analysis.
         */
        @SuppressWarnings("unchecked")
        void simulateTransientBody(double timeHorizon) {
            transientMode = true;
            mserEnabled = false;  // Disable MSER-5 for transient analysis
            ssjSim.init();

            // Warn that timespan is used instead of samples for transient analysis
            if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                System.out.println("LDES transient: using options.timespan=["
                        + options.timespan[0] + ", " + options.timespan[1]
                        + "] (options.samples ignored)");
            }

            // Disable event-count based stopping for transient mode
            maxEvents = Long.MAX_VALUE;
            totalEventCount = 0L;
            totalSimEvents = 0L;
            maxSimEventsLimit = -1L;  // No event limit for transient mode

            // Initialize trace logging if DEBUG verbose level
            initializeTracing();

            // Initialize Logger nodes for CSV output
            initializeLoggers();

            // Initialize routing strategies
            initializeRouting();

            // Initialize random generators
            initializeGenerators();

            // Initialize cache samplers (requires routingRng to be initialized)
            initializeCacheSamplers();

            // Initialize queues and statistics
            initializeState();

            // Validate closed class routing (no routes to Sink)
            validateClosedClassRouting();

            // Initialize queue populations - use init_sol if provided, otherwise reference stations
            if (initSol != null && !initSol.isEmpty()) {
                initFromInitSol();
            } else {
                initClosedClassPopulations();
            }

            // Initialize transient data collection
            initializeTransient(timeHorizon);

            // Initialize convergence checking structures (but disable for transient mode)
            initializeConvergence(Long.MAX_VALUE);
            convergenceEnabled = false;  // Force disable for transient mode

            // Mark warmup as done immediately (no warmup for transient)
            warmupDone = true;
            warmupEndTime = 0.0;

            // Schedule periodic sampling for transient metrics
            new TransientSampleEvent().schedule(transientSamplingInterval);

            // Schedule end of simulation
            new EndOfSimulation().schedule(timeHorizon);

            // Schedule initial arrivals from each source for OPEN classes only
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    if (isOpenClass[k]) {
                        if (isMarkedNonCarrier(srcIdx, k)) {
                            continue; // driven by the marked group's carrier class
                        }
                        double interarrivalTime = generateInterarrivalTime(srcIdx, k);
                        if (interarrivalTime > 0) {
                            new ExternalArrival(srcIdx, k).schedule(interarrivalTime);
                        }
                    }
                }
            }

            // Start service for any initial tokens sitting in queueing places' embedded queues.
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                if (isQueueingPlace[placeListIdx]) {
                    tryStartPlaceService(placeListIdx);
                }
            }

            // Run simulation
            ssjSim.start();
        }

        /**
         * Validate that closed classes do not route to Sink nodes.
         */
        private void validateClosedClassRouting() {
            for (int k = 0; k < numClasses; k++) {
                if (isClosedClass[k]) {
                    for (int fromNode = 0; fromNode < numNodes; fromNode++) {
                        for (Integer sinkNode : sinkNodes) {
                            int fromIdx = fromNode * numClasses + k;
                            int toIdx = sinkNode * numClasses + k;
                            if (fromIdx < sn.rtnodes.getNumRows()
                                    && toIdx < sn.rtnodes.getNumCols()) {
                                double prob = sn.rtnodes.get(fromIdx, toIdx);
                                if (prob > 0) {
                                    throw new RuntimeException(
                                            "solver_ssj: Closed class " + k
                                                    + " has routing to Sink "
                                                    + "(violates population conservation)");
                                }
                            }
                        }
                    }
                }
            }
        }

        /**
         * Initialize closed class populations by injecting jobs at their reference stations.
         */
        private void initClosedClassPopulations() {
            for (int k = 0; k < numClasses; k++) {
                if (isClosedClass[k]) {
                    int population = closedClassPopulation[k];
                    int refStationIdx = referenceStation[k];

                    // Find the service node index for this reference station
                    int refNodeIdx = (int) sn.stationToNode.get(refStationIdx);
                    NodeType nodeType = sn.nodetype.get(refNodeIdx);

                    // Check if reference station is a Queue/Delay or a Place node
                    if (nodeType == NodeType.Place) {
                        // For Petri nets: tokens already initialized; just check transitions
                        if (!placeNodes.isEmpty() && !transitionNodes.isEmpty()) {
                            checkAndFireTransitions();
                        }
                        continue;
                    }

                    int queueIdx = serviceNodes.indexOf(refNodeIdx);

                    if (queueIdx < 0) {
                        throw new RuntimeException(
                                "solver_ssj: Reference station " + refStationIdx
                                        + " for closed class " + k
                                        + " is not a service node (Queue, Delay, or Place)");
                    }

                    // Inject all jobs at reference station at time 0. If the reference
                    // station has finite capacity and cannot admit the full population,
                    // spill the remainder to subsequent service nodes: silently dropping
                    // a closed-class job would break population conservation and freeze
                    // the simulation at a lower effective population.
                    boolean needsJobId = synchCallReplyClass[k] >= 0;
                    for (int jobIdx = 0; jobIdx < population; jobIdx++) {
                        long jobId = needsJobId ? nextJobId++ : -1L;
                        Customer customer = new Customer(
                                k,
                                classPrio[k],
                                0.0,
                                0.0,
                                siroRng.nextDouble(),
                                -1.0,
                                jobId,
                                0.0 + classDeadline[k],
                                -1,
                                null);
                        boolean admitted = false;
                        for (int off = 0; off < numServiceNodes && !admitted; off++) {
                            int qIdx = (queueIdx + off) % numServiceNodes;
                            if (getTotalCustomersAtStation(qIdx) < bufferCapacities[qIdx]) {
                                admitted = arriveAtQueue(qIdx, customer);
                            }
                        }
                        if (!admitted) {
                            throw new RuntimeException(
                                    "solver_ssj: Cannot place closed class " + k
                                            + " population " + population
                                            + " - total finite capacity is insufficient");
                        }
                    }
                }
            }
        }

        /**
         * Initialize queue populations from init_sol matrix.
         */
        private void initFromInitSol() {
            if (initSol == null || initSol.isEmpty()) {
                // Fallback to standard initialization
                initClosedClassPopulations();
                return;
            }
            if (!placeNodes.isEmpty()) {
                // Petri-net models keep their tokens in Places, not in the
                // service-node queues that init_sol describes: the station-based
                // parser below would place nothing and the closed-population
                // conservation check would fail. Use the SPN marking-based
                // initialization instead.
                initClosedClassPopulations();
                return;
            }

            // Parse init_sol: row-major [station0_class0, station0_class1, ..., stationM-1_classK-1]
            int idx = 0;
            for (int stationIdx = 0; stationIdx < numStations; stationIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    int initialQueueLength;
                    if (idx < initSol.length()) {
                        initialQueueLength = (int) FastMath.max(0.0, initSol.get(idx));
                    } else {
                        initialQueueLength = 0;
                    }
                    idx++;

                    if (initialQueueLength <= 0) continue;

                    // Find the service node index for this station
                    int nodeIdx = (int) sn.stationToNode.get(stationIdx);
                    NodeType nodeType = sn.nodetype.get(nodeIdx);

                    // Skip non-service nodes (Source, Sink, Router, etc.)
                    if (nodeType != NodeType.Queue && nodeType != NodeType.Delay) {
                        continue;
                    }

                    int queueIdx = serviceNodes.indexOf(nodeIdx);
                    if (queueIdx < 0) {
                        // This station is not a service node in our list
                        continue;
                    }

                    // Inject jobs at this station at time 0
                    boolean needsJobId = synchCallReplyClass[k] >= 0;
                    for (int jobIdx = 0; jobIdx < initialQueueLength; jobIdx++) {
                        long jobId = needsJobId ? nextJobId++ : -1L;
                        Customer customer = new Customer(
                                k,
                                classPrio[k],
                                0.0,
                                0.0,
                                siroRng.nextDouble(),
                                -1.0,
                                jobId,
                                0.0 + classDeadline[k],
                                -1,
                                null);
                        arriveAtQueue(queueIdx, customer);
                    }
                }
            }

            // For closed classes, verify population conservation
            for (int k = 0; k < numClasses; k++) {
                if (isClosedClass[k]) {
                    int totalInSystem = 0;
                    for (int queueIdx = 0; queueIdx < numServiceNodes; queueIdx++) {
                        totalInSystem += currentQueueLength[queueIdx][k];
                    }
                    if (totalInSystem != closedClassPopulation[k]) {
                        throw new RuntimeException(
                                "init_sol population mismatch for closed class " + k + ": "
                                        + "expected " + closedClassPopulation[k]
                                        + ", got " + totalInSystem);
                    }
                }
            }
        }

        /**
         * Initialize transient data collection structures.
         */
        @SuppressWarnings("unchecked")
        private void initializeTransient(double timeHorizon) {
            // Target ~100 time points for smooth curves
            int targetObservations = 100;
            transientSamplingInterval = timeHorizon / targetObservations;

            transientTimes = new ArrayList<Double>();
            transientQueueLengths = new List[numServiceNodes][numClasses];
            transientUtilizations = new List[numServiceNodes][numClasses];
            transientThroughputs = new List[numServiceNodes][numClasses];
            transientCompletions = new List[numServiceNodes][numClasses];
            for (int q = 0; q < numServiceNodes; q++) {
                for (int k = 0; k < numClasses; k++) {
                    transientQueueLengths[q][k] = new ArrayList<Double>();
                    transientUtilizations[q][k] = new ArrayList<Double>();
                    transientThroughputs[q][k] = new ArrayList<Double>();
                    transientCompletions[q][k] = new ArrayList<Integer>();
                }
            }

            lastTransientSampleTime = 0.0;
            lastTransientQueueTime = new double[numServiceNodes][numClasses];
            lastTransientBusyTime = new double[numServiceNodes][numClasses];
            lastTransientCompletions = new int[numServiceNodes][numClasses];

            // Record initial state at time 0
            transientTimes.add(0.0);
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    transientQueueLengths[qIdx][k].add(
                            (double) currentQueueLength[qIdx][k]);
                    int busyServers = currentBusyServers[qIdx][k];
                    double initialUtil;
                    if (isDelayNode.get(qIdx)) {
                        initialUtil = (double) currentQueueLength[qIdx][k];
                    } else {
                        initialUtil = ((double) busyServers) / numServers[qIdx];
                    }
                    transientUtilizations[qIdx][k].add(initialUtil);
                    transientThroughputs[qIdx][k].add(0.0);  // No throughput at t=0
                    transientCompletions[qIdx][k].add(0);
                }
            }
        }

        /**
         * Transient sampling event - records metrics at regular intervals.
         */
        private final class TransientSampleEvent extends SimEvent {
            @Override
            protected double slotPhase() {
                return SLOT_PHASE_BOOKKEEPING;
            }

            @Override
            public void actions() {
                double currentTime = ssjSim.time();
                double intervalDuration = currentTime - lastTransientSampleTime;

                // Update queue stats before sampling
                for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                    for (int k = 0; k < numClasses; k++) {
                        updateQueueStats(qIdx, k);
                        updateBusyStats(qIdx, k);
                    }
                }

                transientTimes.add(currentTime);

                // Record time-weighted metrics over this interval
                for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                    for (int k = 0; k < numClasses; k++) {
                        if (intervalDuration > 0) {
                            double intervalQueueTime = totalQueueTime[qIdx][k]
                                    - lastTransientQueueTime[qIdx][k];
                            double avgQueueLength = intervalQueueTime / intervalDuration;
                            transientQueueLengths[qIdx][k].add(avgQueueLength);
                            lastTransientQueueTime[qIdx][k] = totalQueueTime[qIdx][k];

                            double intervalBusyTime = totalBusyTime[qIdx][k]
                                    - lastTransientBusyTime[qIdx][k];
                            double avgUtilization;
                            if (isDelayNode.get(qIdx)) {
                                avgUtilization = avgQueueLength;
                            } else {
                                avgUtilization = intervalBusyTime
                                        / (intervalDuration * numServers[qIdx]);
                            }
                            transientUtilizations[qIdx][k].add(avgUtilization);
                            lastTransientBusyTime[qIdx][k] = totalBusyTime[qIdx][k];

                            int intervalCompletions = completedCustomers[qIdx][k]
                                    - lastTransientCompletions[qIdx][k];
                            double throughput = ((double) intervalCompletions) / intervalDuration;
                            transientThroughputs[qIdx][k].add(throughput);
                            lastTransientCompletions[qIdx][k] = completedCustomers[qIdx][k];
                        } else {
                            transientQueueLengths[qIdx][k].add(
                                    (double) currentQueueLength[qIdx][k]);
                            transientUtilizations[qIdx][k].add(0.0);
                            transientThroughputs[qIdx][k].add(0.0);
                        }
                        transientCompletions[qIdx][k].add(completedCustomers[qIdx][k]);
                    }
                }
                lastTransientSampleTime = currentTime;

                // Schedule next sample if not at end
                if (currentTime + transientSamplingInterval < Double.MAX_VALUE) {
                    new TransientSampleEvent().schedule(transientSamplingInterval);
                }
            }
        }

        /**
         * Get transient analysis results.
         */
        LDESResult getTransientLDESResultBody() {
            LDESResult result = new LDESResult();
            result.sn = sn;

            // Create time points matrix
            int numTimePoints = transientTimes.size();
            result.t = new Matrix(numTimePoints, 1);
            for (int i = 0; i < numTimePoints; i++) {
                result.t.set(i, 0, transientTimes.get(i));
            }

            // Create transient matrices [stations x classes] at each time point
            result.QNt = new Matrix[numStations][numClasses];
            result.UNt = new Matrix[numStations][numClasses];
            result.TNt = new Matrix[numStations][numClasses];
            for (int i = 0; i < numStations; i++) {
                for (int k = 0; k < numClasses; k++) {
                    result.QNt[i][k] = new Matrix(0, 0);
                    result.UNt[i][k] = new Matrix(0, 0);
                    result.TNt[i][k] = new Matrix(0, 0);
                }
            }

            // Fill in results for service nodes
            for (int svcIdx = 0; svcIdx < serviceStations.size(); svcIdx++) {
                int serviceStation = serviceStations.get(svcIdx);
                for (int k = 0; k < numClasses; k++) {
                    if (svcIdx < transientQueueLengths.length
                            && k < transientQueueLengths[svcIdx].length) {
                        List<Double> qObs = transientQueueLengths[svcIdx][k];
                        List<Double> uObs = transientUtilizations[svcIdx][k];
                        List<Double> tObs = transientThroughputs[svcIdx][k];
                        int obsCount = Math.min(qObs.size(), transientTimes.size());

                        Matrix qMatrix = new Matrix(obsCount, 2);
                        Matrix uMatrix = new Matrix(obsCount, 2);
                        Matrix tMatrix = new Matrix(obsCount, 2);

                        for (int i = 0; i < obsCount; i++) {
                            qMatrix.set(i, 0, qObs.get(i));
                            qMatrix.set(i, 1, transientTimes.get(i));
                            uMatrix.set(i, 0, uObs.get(i));
                            uMatrix.set(i, 1, transientTimes.get(i));
                            tMatrix.set(i, 0, tObs.get(i));
                            tMatrix.set(i, 1, transientTimes.get(i));
                        }

                        result.QNt[serviceStation][k] = qMatrix;
                        result.UNt[serviceStation][k] = uMatrix;
                        result.TNt[serviceStation][k] = tMatrix;
                    }
                }
            }

            // Also compute final steady-state estimates (last interval averages)
            Matrix QN = new Matrix(numStations, numClasses);
            Matrix UN = new Matrix(numStations, numClasses);
            Matrix RN = new Matrix(numStations, numClasses);
            Matrix TN = new Matrix(numStations, numClasses);
            Matrix CN = new Matrix(1, numClasses);
            Matrix XN = new Matrix(1, numClasses);

            // Build set of classes that can receive jobs via class-switching
            Set<Integer> classSwitchClasses2 = new HashSet<Integer>();
            for (int c = 0; c < sn.nchains; c++) {
                List<Integer> classesInChain = new ArrayList<Integer>();
                for (int k = 0; k < numClasses; k++) {
                    if (sn.chains.get(c, k) > 0) {
                        classesInChain.add(k);
                    }
                }
                if (classesInChain.size() > 1) {
                    boolean chainHasJobs = false;
                    for (Integer k : classesInChain) {
                        if (sn.njobs.get(k) > 0) {
                            chainHasJobs = true;
                            break;
                        }
                    }
                    if (chainHasJobs) {
                        classSwitchClasses2.addAll(classesInChain);
                    }
                }
            }

            // Spawn targets are live even in zero-population chains (LQN
            // phase-2 continuations): unmask their chains.
            if (spawnClassOf != null) {
                for (int k = 0; k < spawnClassOf.length; k++) {
                    int sc = spawnClassOf[k];
                    if (sc < 0) continue;
                    for (int c = 0; c < sn.nchains; c++) {
                        if (sn.chains.get(c, sc) > 0) {
                            for (int k2 = 0; k2 < numClasses; k2++) {
                                if (sn.chains.get(c, k2) > 0) {
                                    classSwitchClasses2.add(k2);
                                }
                            }
                        }
                    }
                }
            }

            for (int svcIdx = 0; svcIdx < serviceStations.size(); svcIdx++) {
                int serviceStation = serviceStations.get(svcIdx);
                for (int k = 0; k < numClasses; k++) {
                    double njobs = sn.njobs.get(k);
                    boolean isClosedWithZeroPopulation = !isOpenClass[k]
                            && Double.isFinite(njobs) && njobs == 0.0;
                    boolean canReceiveViaClassSwitch2 = classSwitchClasses2.contains(k);
                    if ((!isClosedWithZeroPopulation || canReceiveViaClassSwitch2)
                            && mus[svcIdx][k] < Double.MAX_VALUE) {
                        QN.set(serviceStation, k, getAvgQueueLength(svcIdx, k));
                        UN.set(serviceStation, k, getUtilization(svcIdx, k));
                        RN.set(serviceStation, k, getAvgResponseTime(svcIdx, k));
                        TN.set(serviceStation, k, getThroughput(svcIdx, k));
                    }
                }
            }

            for (int srcIdx = 0; srcIdx < sourceStations.size(); srcIdx++) {
                int sourceStation = sourceStations.get(srcIdx);
                for (int k = 0; k < numClasses; k++) {
                    if (lambdas[srcIdx][k] > 0) {
                        TN.set(sourceStation, k, lambdas[srcIdx][k]);
                    }
                }
            }

            for (int k = 0; k < numClasses; k++) {
                double totalArrival = 0.0;
                for (int srcIdx = 0; srcIdx < sourceStations.size(); srcIdx++) {
                    totalArrival += lambdas[srcIdx][k];
                }
                if (totalArrival > 0) {
                    XN.set(0, k, getSystemThroughput(k));
                    CN.set(0, k, getSystemResponseTime(k));
                }
            }

            result.QN = QN;
            result.UN = UN;
            result.RN = RN;
            result.TN = TN;
            result.CN = CN;
            result.XN = XN;

            // Store response time samples for CDF computation
            setRespTimeSamples(result, numStations, numClasses);

            return result;
        }

        // ---------------------------------------------------------------------
        // MSER initialization (Kotlin lines 3722-3812).
        // ---------------------------------------------------------------------

        /**
         * Initialize MSER-5 data collection structures.
         */
        @SuppressWarnings("unchecked")
        private void initializeMSER(double timeHorizon) {
            // Target ~1000 observations for good MSER-5 estimation (200 batches of 5)
            int targetObservations = 1000;
            mserSamplingInterval = timeHorizon / targetObservations;

            queueLengthObservations = new List[numServiceNodes][numClasses];
            throughputObservations = new List[numServiceNodes][numClasses];
            busyTimeObservations = new List[numServiceNodes][numClasses];
            queueTimeObservations = new List[numServiceNodes][numClasses];
            blockingTimeObservations = new List[numServiceNodes][numClasses];
            placeCompletionObservations = new List[placeNodes.size()][numClasses];
            for (int q = 0; q < numServiceNodes; q++) {
                for (int k = 0; k < numClasses; k++) {
                    queueLengthObservations[q][k] = new ArrayList<Double>();
                    throughputObservations[q][k] = new ArrayList<Integer>();
                    busyTimeObservations[q][k] = new ArrayList<Double>();
                    queueTimeObservations[q][k] = new ArrayList<Double>();
                    blockingTimeObservations[q][k] = new ArrayList<Double>();
                }
            }
            for (int p = 0; p < placeNodes.size(); p++) {
                for (int k = 0; k < numClasses; k++) {
                    placeCompletionObservations[p][k] = new ArrayList<Integer>();
                }
            }
            observationTimes = new ArrayList<Double>();

            // Initialize MSER tracking arrays
            lastMserQueueTime = new double[numServiceNodes][numClasses];
            lastMserSampleTime = 0.0;
        }

        /**
         * Initialize MSER-5 data collection structures for event-based stopping.
         */
        @SuppressWarnings("unchecked")
        private void initializeMSEREventBased(long maxEventCount) {
            // Target ~1000 observations for good MSER-5 estimation (200 batches of 5)
            int targetObservations = 1000;
            mserEventInterval = maxEventCount / targetObservations;
            if (mserEventInterval < 1) mserEventInterval = 1;
            lastMserEventCount = 0L;

            queueLengthObservations = new List[numServiceNodes][numClasses];
            throughputObservations = new List[numServiceNodes][numClasses];
            busyTimeObservations = new List[numServiceNodes][numClasses];
            queueTimeObservations = new List[numServiceNodes][numClasses];
            blockingTimeObservations = new List[numServiceNodes][numClasses];
            placeCompletionObservations = new List[placeNodes.size()][numClasses];
            for (int q = 0; q < numServiceNodes; q++) {
                for (int k = 0; k < numClasses; k++) {
                    queueLengthObservations[q][k] = new ArrayList<Double>();
                    throughputObservations[q][k] = new ArrayList<Integer>();
                    busyTimeObservations[q][k] = new ArrayList<Double>();
                    queueTimeObservations[q][k] = new ArrayList<Double>();
                    blockingTimeObservations[q][k] = new ArrayList<Double>();
                }
            }
            for (int p = 0; p < placeNodes.size(); p++) {
                for (int k = 0; k < numClasses; k++) {
                    placeCompletionObservations[p][k] = new ArrayList<Integer>();
                }
            }
            observationTimes = new ArrayList<Double>();

            // Initialize MSER tracking arrays
            lastMserQueueTime = new double[numServiceNodes][numClasses];
            lastMserSampleTime = 0.0;
        }

        /**
         * Initialize transient detection and CI configuration from LDESOptions.
         */
        private void initializeTransientAndCIConfig() {
            LDESOptions ldesOptions = (options instanceof LDESOptions)
                    ? (LDESOptions) options : null;

            // Transient detection configuration
            String tranfilter = (ldesOptions != null) ? ldesOptions.tranfilter : "mser5";
            effectiveMserBatchSize = (ldesOptions != null)
                    ? ldesOptions.mserbatch : LDESOptions.DEFAULT_MSER_BATCH;
            effectiveWarmupFraction = (ldesOptions != null)
                    ? ldesOptions.warmupfrac : LDESOptions.DEFAULT_WARMUP_FRAC;

            // Derive mserEnabled from tranfilter
            mserEnabled = "mser5".equals(tranfilter);

            // CI configuration
            effectiveCiMethod = (ldesOptions != null) ? ldesOptions.cimethod : "obm";
            effectiveObmOverlap = (ldesOptions != null)
                    ? ldesOptions.obmoverlap : LDESOptions.DEFAULT_OBM_OVERLAP;
            effectiveCiMinBatch = (ldesOptions != null)
                    ? ldesOptions.ciminbatch : LDESOptions.DEFAULT_CI_MIN_BATCH;
            effectiveCiMinObs = (ldesOptions != null)
                    ? ldesOptions.ciminobs : LDESOptions.DEFAULT_CI_MIN_OBS;
            effectiveSpectralLowFreqFrac = (ldesOptions != null)
                    ? ldesOptions.spectralLowFreqFrac : LDESOptions.DEFAULT_SPECTRAL_LOW_FREQ_FRAC;

            // Transient-filter and CI-method notices are DEBUG-only: at STD the
            // LDES output mirrors SSA, which prints just the sample-count line
            // and the analysis summary with no filter/estimator preamble.
            if (options.verbose == VerboseLevel.DEBUG) {
                String suffix;
                if ("mser5".equals(tranfilter)) {
                    suffix = " (batch size = " + effectiveMserBatchSize + ")";
                } else if ("fixed".equals(tranfilter)) {
                    suffix = " (warmup fraction = " + effectiveWarmupFraction + ")";
                } else {
                    suffix = "";
                }
                System.out.println("LDES: Transient filter = " + tranfilter + suffix);
                String ciSuffix = "obm".equals(effectiveCiMethod)
                        ? " (overlap = " + effectiveObmOverlap + ")" : "";
                System.out.println("LDES: CI method = " + effectiveCiMethod + ciSuffix);
            }
        }

        /**
         * Initialize convergence checking data structures.
         */
        @SuppressWarnings("unchecked")
        private void initializeConvergence(long maxEventCount) {
            // Read convergence options (cast to LDESOptions if available)
            LDESOptions ldesOptions = (options instanceof LDESOptions)
                    ? (LDESOptions) options : null;
            convergenceEnabled = (ldesOptions != null) ? ldesOptions.cnvgon : false;
            convergenceTolerance = (ldesOptions != null) ? ldesOptions.cnvgtol : 0.05;
            convergenceMinBatches = (ldesOptions != null) ? ldesOptions.cnvgbatch : 20;

            // Auto-calculate check interval if not specified
            int configuredInterval = (ldesOptions != null) ? ldesOptions.cnvgchk : 0;
            if (configuredInterval > 0) {
                convergenceCheckInterval = (long) configuredInterval;
            } else {
                convergenceCheckInterval = Math.max(1L, maxEventCount / 50);
            }
            lastConvergenceCheckEventCount = 0L;

            // Initialize batch means arrays
            queueBatchMeans = new List[numServiceNodes][numClasses];
            utilBatchMeans = new List[numServiceNodes][numClasses];
            respTimeBatchMeans = new List[numServiceNodes][numClasses];
            throughputBatchMeans = new List[numServiceNodes][numClasses];
            for (int q = 0; q < numServiceNodes; q++) {
                for (int k = 0; k < numClasses; k++) {
                    queueBatchMeans[q][k] = new ArrayList<Double>();
                    utilBatchMeans[q][k] = new ArrayList<Double>();
                    respTimeBatchMeans[q][k] = new ArrayList<Double>();
                    throughputBatchMeans[q][k] = new ArrayList<Double>();
                }
            }

            // Initialize batch accumulators
            batchStartQueueTime = new double[numServiceNodes][numClasses];
            batchStartBusyTime = new double[numServiceNodes][numClasses];
            batchStartCompletions = new int[numServiceNodes][numClasses];
            currentBatchRespTimeSum = new double[numServiceNodes][numClasses];
            currentBatchRespTimeCount = new int[numServiceNodes][numClasses];

            currentBatchObservations = 0;
            batchStartTime = 0.0;
            hasConverged = false;
            stoppingReason = "max_events";

            // Initialize final CI arrays
            finalQNCI = new double[numServiceNodes][numClasses];
            finalUNCI = new double[numServiceNodes][numClasses];
            finalRNCI = new double[numServiceNodes][numClasses];
            finalTNCI = new double[numServiceNodes][numClasses];
            finalQNRelPrec = new double[numServiceNodes][numClasses];
            finalUNRelPrec = new double[numServiceNodes][numClasses];
            finalRNRelPrec = new double[numServiceNodes][numClasses];
            finalTNRelPrec = new double[numServiceNodes][numClasses];
        }

        /**
         * Track a simulation event and check if the total event limit has been reached.
         */
        private void trackEvent() {
            totalSimEvents++;
            if (maxSimEventsLimit > 0 && totalSimEvents >= maxSimEventsLimit) {
                stoppingReason = "max_sim_events";
                finishSimulation();
                return;
            }
            // Cooperative wall-clock time budget (LDESOptions.maxTime, seconds).
            if (!Double.isInfinite(maxTimeLimit) && maxTimeLimit > 0
                    && (System.nanoTime() - simStartNanos) / 1e9 > maxTimeLimit) {
                stoppingReason = "max_time";
                finishSimulation();
            }
        }

        /**
         * Check if event count thresholds have been reached and take appropriate action.
         */
        private void checkEventCountStop() {
            // Firings that settle the initial marking are model setup, not
            // simulation events: skip counting and all stop/sampling checks
            // (the MSER/convergence samplers are not yet initialized here).
            if (initializing) {
                return;
            }
            totalEventCount++;

            // Cooperative wall-clock time budget (LDESOptions.maxTime /
            // SolverOptions.timeout, seconds). Checked every event so the
            // simulation stops promptly once the budget is exceeded.
            if (!Double.isInfinite(maxTimeLimit) && maxTimeLimit > 0
                    && (System.nanoTime() - simStartNanos) / 1e9 > maxTimeLimit) {
                stoppingReason = "max_time";
                finishSimulation();
                return;
            }

            // Handle warmup completion (non-MSER mode)
            if (!mserEnabled && !warmupDone && totalEventCount >= warmupEventThreshold) {
                resetStatistics();
                warmupDone = true;
            }

            // Handle MSER sampling
            if (mserEnabled && (totalEventCount - lastMserEventCount) >= mserEventInterval) {
                collectMSERSample();
                lastMserEventCount = totalEventCount;
            }

            // Handle convergence sampling and checking at configured intervals
            if (convergenceEnabled
                    && (totalEventCount - lastConvergenceCheckEventCount) >= convergenceCheckInterval) {
                // Collect convergence sample
                collectConvergenceSample();
                lastConvergenceCheckEventCount = totalEventCount;

                // Check for convergence
                if (checkConvergence()) {
                    hasConverged = true;
                    stoppingReason = "convergence";
                    if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                        System.out.println();
                        System.out.println("LDES: Convergence detected at event " + totalEventCount);
                    }
                    finishSimulation();
                    return;
                }
            }

            // Handle progress reporting
            long progressInterval = maxEvents / 50;
            if (progressInterval > 0
                    && (totalEventCount - lastProgressEventCount) >= progressInterval) {
                if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                    if (lastProgressEventCount == 0L) {
                        System.out.printf("LDES events: %6d ", totalEventCount);
                        System.out.flush();
                    } else {
                        System.out.printf("\b\b\b\b\b\b\b %6d", totalEventCount);
                        System.out.flush();
                    }
                }
                lastProgressEventCount = totalEventCount;
            }

            // Handle streaming if collector is active
            if (stream != null) {
                jline.streaming.StreamingOptions streamOpts = stream.getOptions();
                double currentTime = ssjSim.time();
                boolean shouldStream;
                if (streamOpts.mode == jline.streaming.StreamingOptions.StreamMode.SAMPLED) {
                    shouldStream = (totalEventCount - lastStreamEventCount)
                            >= streamOpts.sampleFrequency;
                } else if (streamOpts.mode == jline.streaming.StreamingOptions.StreamMode.TIME_WINDOW) {
                    shouldStream = (currentTime - lastStreamTime)
                            >= streamOpts.timeWindowSeconds;
                } else {
                    shouldStream = false;
                }
                if (shouldStream) {
                    pushStreamingMetrics(currentTime);
                    lastStreamEventCount = totalEventCount;
                    lastStreamTime = currentTime;
                }
            }

            // Check if max events reached
            if (totalEventCount >= maxEvents) {
                stoppingReason = "max_events";
                finishSimulation();
            }
        }

        /**
         * Push current queue state metrics to the streaming collector.
         */
        private void pushStreamingMetrics(double currentTime) {
            if (stream == null) return;

            // Build queue length matrix from current state
            Matrix nir = new Matrix(numServiceNodes, numClasses);
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    nir.set(qIdx, k, (double) currentQueueLength[qIdx][k]);
                }
            }

            // Calculate time delta since last push (or from simulation start)
            double dt = (lastStreamTime > 0.0) ? currentTime - lastStreamTime : currentTime;

            // Record state with the streaming collector
            stream.recordState(currentTime, dt, nir, null, null);
        }

        /**
         * Collect an MSER sample based on current state.
         */
        private void collectMSERSample() {
            double currentTime = ssjSim.time();
            double intervalDuration = currentTime - lastMserSampleTime;

            // Update queue and busy stats before sampling
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                if (isPSScheduling(schedStrategies[qIdx])) {
                    updatePSBusyStats(qIdx);
                }
                for (int k = 0; k < numClasses; k++) {
                    updateQueueStats(qIdx, k);
                    if (!isPSScheduling(schedStrategies[qIdx])) {
                        updateBusyStats(qIdx, k);
                    }
                }
            }

            observationTimes.add(currentTime);

            // Record time-weighted average queue length over this interval
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    if (intervalDuration > 0) {
                        double intervalQueueTime = totalQueueTime[qIdx][k]
                                - lastMserQueueTime[qIdx][k];
                        double avgQueueLength = intervalQueueTime / intervalDuration;
                        queueLengthObservations[qIdx][k].add(avgQueueLength);
                        lastMserQueueTime[qIdx][k] = totalQueueTime[qIdx][k];
                    } else {
                        // Include blocked jobs (BAS, BBS) but NOT FCR-blocked trying to enter this queue
                        int effectiveQueueLength = currentQueueLength[qIdx][k]
                                + basBlockedAtDest[qIdx][k] + bbsBlockedAtDest[qIdx][k];
                        queueLengthObservations[qIdx][k].add((double) effectiveQueueLength);
                    }
                    throughputObservations[qIdx][k].add(completedCustomers[qIdx][k]);
                    busyTimeObservations[qIdx][k].add(totalBusyTime[qIdx][k]);
                    queueTimeObservations[qIdx][k].add(totalQueueTime[qIdx][k]);
                    blockingTimeObservations[qIdx][k].add(totalBlockingTime[qIdx][k]);
                }
            }

            // Record Place completion counts for MSER truncation
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    placeCompletionObservations[placeListIdx][k]
                            .add(placeCompletions[placeListIdx][k]);
                }
            }

            lastMserSampleTime = currentTime;
        }

        /**
         * Collect a convergence sample based on current state.
         */
        private void collectConvergenceSample() {
            if (!convergenceEnabled) return;

            double currentTime = ssjSim.time();

            // Update stats before sampling
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                if (isPSScheduling(schedStrategies[qIdx])) {
                    updatePSBusyStats(qIdx);
                }
                for (int k = 0; k < numClasses; k++) {
                    updateQueueStats(qIdx, k);
                    if (!isPSScheduling(schedStrategies[qIdx])) {
                        updateBusyStats(qIdx, k);
                    }
                }
            }

            currentBatchObservations++;

            // Check if batch is complete (use same batch size as MSER)
            if (currentBatchObservations >= effectiveMserBatchSize) {
                finalizeBatch(currentTime);
            }
        }

        /**
         * Finalize a batch - compute batch means for all metrics and add to lists.
         */
        private void finalizeBatch(double currentTime) {
            double batchDuration = currentTime - batchStartTime;
            if (batchDuration <= 0) {
                // Reset batch without recording
                resetBatchAccumulators(currentTime);
                return;
            }

            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    // Queue length batch mean (time-weighted)
                    double queueTimeDelta = totalQueueTime[qIdx][k]
                            - batchStartQueueTime[qIdx][k];
                    double avgQueueLength = queueTimeDelta / batchDuration;
                    queueBatchMeans[qIdx][k].add(avgQueueLength);

                    // Per-server Utilization batch mean (time-weighted)
                    double busyTimeDelta = totalBusyTime[qIdx][k]
                            - batchStartBusyTime[qIdx][k];
                    double avgUtil;
                    if (isDelayNode.get(qIdx)) {
                        avgUtil = avgQueueLength;
                    } else {
                        avgUtil = busyTimeDelta / (batchDuration * numServers[qIdx]);
                    }
                    utilBatchMeans[qIdx][k].add(avgUtil);

                    // Response time batch mean
                    int respTimeCount = currentBatchRespTimeCount[qIdx][k];
                    double avgRespTime;
                    if (respTimeCount > 0) {
                        avgRespTime = currentBatchRespTimeSum[qIdx][k] / respTimeCount;
                    } else {
                        // Fallback: use overall average if available
                        if (responseTimeTally[qIdx][k] != null
                                && responseTimeTally[qIdx][k].numberObs() > 0) {
                            avgRespTime = responseTimeTally[qIdx][k].average();
                        } else {
                            avgRespTime = 0.0;
                        }
                    }
                    respTimeBatchMeans[qIdx][k].add(avgRespTime);

                    // Throughput batch mean
                    int completionsDelta = completedCustomers[qIdx][k]
                            - batchStartCompletions[qIdx][k];
                    double avgThroughput = ((double) completionsDelta) / batchDuration;
                    throughputBatchMeans[qIdx][k].add(avgThroughput);
                }
            }

            // Reset batch accumulators for next batch
            resetBatchAccumulators(currentTime);
        }

        /**
         * Reset batch accumulators for the start of a new batch.
         */
        private void resetBatchAccumulators(double currentTime) {
            batchStartTime = currentTime;
            currentBatchObservations = 0;

            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    batchStartQueueTime[qIdx][k] = totalQueueTime[qIdx][k];
                    batchStartBusyTime[qIdx][k] = totalBusyTime[qIdx][k];
                    batchStartCompletions[qIdx][k] = completedCustomers[qIdx][k];
                    currentBatchRespTimeSum[qIdx][k] = 0.0;
                    currentBatchRespTimeCount[qIdx][k] = 0;
                }
            }
        }

        /**
         * Record a response time observation for the current batch.
         */
        private void recordResponseTimeForBatch(int qIdx, int k, double responseTime) {
            if (!convergenceEnabled) return;
            if (qIdx < numServiceNodes && k < numClasses) {
                currentBatchRespTimeSum[qIdx][k] += responseTime;
                currentBatchRespTimeCount[qIdx][k]++;
            }
        }

        /**
         * Finish the simulation - called when max events reached.
         */
        private void finishSimulation() {
            // Final update of queue and busy statistics
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                if (isPSScheduling(schedStrategies[qIdx])) {
                    updatePSBusyStats(qIdx);
                    for (int k = 0; k < numClasses; k++) {
                        updateQueueStats(qIdx, k);
                    }
                } else {
                    for (int k = 0; k < numClasses; k++) {
                        updateQueueStats(qIdx, k);
                        updateBusyStats(qIdx, k);
                    }
                }
            }

            // Apply MSER-5 truncation to determine warmup period
            if (mserEnabled) {
                applyMSER5Truncation();
            }

            // Compute final CI matrices for convergence results
            computeFinalCIMatrices();

            // Close trace writer
            closeTracing();

            // Close Logger file writers
            closeLoggers();

            if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                // Print final event count before newline
                System.out.printf("\b\b\b\b\b\b\b %6d", totalEventCount);
                System.out.println();
            }

            ssjSim.stop();
        }

        /**
         * Check if all metrics have converged.
         */
        private boolean checkConvergence() {
            if (!convergenceEnabled) return false;

            // Check if we have enough batches
            int numBatches;
            if (queueBatchMeans.length > 0 && queueBatchMeans[0].length > 0) {
                numBatches = queueBatchMeans[0][0].size();
            } else {
                numBatches = 0;
            }

            if (numBatches < convergenceMinBatches) return false;

            // Check all metrics for all active station/class combinations
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    // Skip inactive combinations (disabled service rate)
                    if (mus[qIdx][k] >= Double.MAX_VALUE) continue;

                    // Check queue length convergence
                    if (!isMetricConverged(queueBatchMeans[qIdx][k],
                            "Q[" + qIdx + "][" + k + "]")) return false;

                    // Check utilization convergence
                    if (!isMetricConverged(utilBatchMeans[qIdx][k],
                            "U[" + qIdx + "][" + k + "]")) return false;

                    // Check response time convergence (skip if no data)
                    List<Double> rtBatch = respTimeBatchMeans[qIdx][k];
                    if (!rtBatch.isEmpty()) {
                        boolean anyPositive = false;
                        for (Double v : rtBatch) {
                            if (v > 0) { anyPositive = true; break; }
                        }
                        if (anyPositive) {
                            if (!isMetricConverged(rtBatch,
                                    "R[" + qIdx + "][" + k + "]")) return false;
                        }
                    }

                    // Check throughput convergence
                    if (!isMetricConverged(throughputBatchMeans[qIdx][k],
                            "T[" + qIdx + "][" + k + "]")) return false;
                }
            }

            return true;
        }

        /**
         * Check if a single metric has converged based on batch means.
         */
        private boolean isMetricConverged(List<Double> batchMeans, String metricName) {
            if (batchMeans.size() < 2) return false;

            // Compute mean
            int n = batchMeans.size();
            double sum = 0.0;
            for (Double x : batchMeans) sum += x;
            double mean = sum / n;

            // Skip if mean is effectively zero (treat as converged)
            if (Math.abs(mean) < 1e-12) return true;

            // Compute sample variance
            double variance = 0.0;
            for (Double x : batchMeans) {
                double diff = x - mean;
                variance += diff * diff;
            }
            variance /= (n - 1);

            // Standard error
            double stdErr = Math.sqrt(variance / n);

            // Get confidence level (default to 0.95 if not configured)
            double confintLevel = (options.confint > 0) ? options.confint : 0.95;
            double alpha = 1.0 - confintLevel;

            // t-critical value
            double tCrit = getTCriticalValue(alpha, n - 1);

            // CI half-width
            double ciHalfWidth = tCrit * stdErr;

            // Relative precision
            double relPrec = ciHalfWidth / Math.abs(mean);

            if (options.verbose == VerboseLevel.DEBUG) {
                System.out.println("Convergence check " + metricName + ": mean=" + mean
                        + ", relPrec=" + relPrec + ", tol=" + convergenceTolerance);
            }

            return relPrec <= convergenceTolerance;
        }

        /**
         * Get t-distribution critical value for given alpha and degrees of freedom.
         */
        private double getTCriticalValue(double alpha, int df) {
            // Common t-critical values for two-tailed test
            double halfAlpha = alpha / 2.0;

            // For large df (>30), t approaches normal distribution
            if (df > 30) {
                if (halfAlpha <= 0.005) return 2.576;  // 99% CI
                if (halfAlpha <= 0.025) return 1.96;   // 95% CI
                if (halfAlpha <= 0.05) return 1.645;   // 90% CI
                return 1.28;                            // 80% CI
            }

            // Table values for smaller df (95% CI values)
            Map<Integer, Double> t95 = new HashMap<Integer, Double>();
            t95.put(1, 12.71); t95.put(2, 4.30); t95.put(3, 3.18); t95.put(4, 2.78); t95.put(5, 2.57);
            t95.put(6, 2.45); t95.put(7, 2.36); t95.put(8, 2.31); t95.put(9, 2.26); t95.put(10, 2.23);
            t95.put(11, 2.20); t95.put(12, 2.18); t95.put(13, 2.16); t95.put(14, 2.14); t95.put(15, 2.13);
            t95.put(16, 2.12); t95.put(17, 2.11); t95.put(18, 2.10); t95.put(19, 2.09); t95.put(20, 2.09);
            t95.put(25, 2.06); t95.put(30, 2.04);

            // For 95% CI (most common case)
            if (halfAlpha >= 0.02 && halfAlpha <= 0.03) {
                Double v = t95.get(df);
                if (v != null) return v;
                v = t95.get(Math.min(df, 30));
                return (v != null) ? v : 2.0;
            }

            // Rough scaling for other confidence levels
            Double baseTBoxed = t95.get(df);
            if (baseTBoxed == null) {
                baseTBoxed = t95.get(Math.min(df, 30));
            }
            double baseT = (baseTBoxed != null) ? baseTBoxed : 2.0;
            if (halfAlpha <= 0.005) return baseT * 1.32;  // 99% CI
            if (halfAlpha <= 0.05) return baseT * 0.84;   // 90% CI
            return baseT * 0.65;                           // 80% CI
        }

        /**
         * Compute final CI matrices for all metrics.
         */
        private void computeFinalCIMatrices() {
            if (!convergenceEnabled) return;

            double confintLevel = (options.confint > 0) ? options.confint : 0.95;
            double alpha = 1.0 - confintLevel;

            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    // Queue length CI
                    jline.util.Pair<Double, Double> qResult = computeCIForMetric(
                            queueBatchMeans[qIdx][k], alpha);
                    finalQNCI[qIdx][k] = qResult.getLeft();
                    finalQNRelPrec[qIdx][k] = qResult.getRight();

                    // Utilization CI
                    jline.util.Pair<Double, Double> uResult = computeCIForMetric(
                            utilBatchMeans[qIdx][k], alpha);
                    finalUNCI[qIdx][k] = uResult.getLeft();
                    finalUNRelPrec[qIdx][k] = uResult.getRight();

                    // Response time CI
                    jline.util.Pair<Double, Double> rResult = computeCIForMetric(
                            respTimeBatchMeans[qIdx][k], alpha);
                    finalRNCI[qIdx][k] = rResult.getLeft();
                    finalRNRelPrec[qIdx][k] = rResult.getRight();

                    // Throughput CI
                    jline.util.Pair<Double, Double> tResult = computeCIForMetric(
                            throughputBatchMeans[qIdx][k], alpha);
                    finalTNCI[qIdx][k] = tResult.getLeft();
                    finalTNRelPrec[qIdx][k] = tResult.getRight();
                }
            }
        }

        /**
         * Compute CI half-width and relative precision for a single metric.
         *
         * @return Pair of (CI half-width, relative precision)
         */
        private jline.util.Pair<Double, Double> computeCIForMetric(
                List<Double> batchMeans, double alpha) {
            if (batchMeans.size() < 2) {
                return new jline.util.Pair<Double, Double>(0.0, 0.0);
            }

            int n = batchMeans.size();
            double sum = 0.0;
            for (Double x : batchMeans) sum += x;
            double mean = sum / n;

            if (Math.abs(mean) < 1e-12) {
                return new jline.util.Pair<Double, Double>(0.0, 0.0);
            }

            double variance = 0.0;
            for (Double x : batchMeans) {
                double diff = x - mean;
                variance += diff * diff;
            }
            variance /= (n - 1);

            double stdErr = Math.sqrt(variance / n);
            double tCrit = getTCriticalValue(alpha, n - 1);
            double ciHalfWidth = tCrit * stdErr;
            double relPrec = ciHalfWidth / Math.abs(mean);

            return new jline.util.Pair<Double, Double>(ciHalfWidth, relPrec);
        }

        /** Cumulative queue time at last MSER sample. */
        private double[][] lastMserQueueTime;
        /** Time of last MSER sample. */
        private double lastMserSampleTime = 0.0;

        /**
         * MSER-5 sampling event - records time-weighted average queue length over each interval.
         */
        private final class MSERSampleEvent extends SimEvent {
            @Override
            protected double slotPhase() {
                return SLOT_PHASE_BOOKKEEPING;
            }

            @Override
            public void actions() {
                double currentTime = ssjSim.time();
                double intervalDuration = currentTime - lastMserSampleTime;

                // Update queue and busy stats before sampling
                for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                    if (isPSScheduling(schedStrategies[qIdx])) {
                        updatePSBusyStats(qIdx);
                    }
                    for (int k = 0; k < numClasses; k++) {
                        updateQueueStats(qIdx, k);
                        if (!isPSScheduling(schedStrategies[qIdx])) {
                            updateBusyStats(qIdx, k);
                        }
                    }
                }

                observationTimes.add(currentTime);

                // Record time-weighted average queue length over this interval
                for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                    for (int k = 0; k < numClasses; k++) {
                        if (intervalDuration > 0) {
                            double intervalQueueTime = totalQueueTime[qIdx][k]
                                    - lastMserQueueTime[qIdx][k];
                            double avgQueueLength = intervalQueueTime / intervalDuration;
                            queueLengthObservations[qIdx][k].add(avgQueueLength);
                            lastMserQueueTime[qIdx][k] = totalQueueTime[qIdx][k];
                        } else {
                            // Include blocked jobs (BAS, BBS) but NOT FCR-blocked trying to enter this queue
                            int effectiveQueueLength = currentQueueLength[qIdx][k]
                                    + basBlockedAtDest[qIdx][k] + bbsBlockedAtDest[qIdx][k];
                            queueLengthObservations[qIdx][k].add((double) effectiveQueueLength);
                        }
                        throughputObservations[qIdx][k].add(completedCustomers[qIdx][k]);
                        busyTimeObservations[qIdx][k].add(totalBusyTime[qIdx][k]);
                        blockingTimeObservations[qIdx][k].add(totalBlockingTime[qIdx][k]);
                    }
                }
                lastMserSampleTime = currentTime;

                // Schedule next sample
                if (currentTime + mserSamplingInterval < Double.MAX_VALUE) {
                    new MSERSampleEvent().schedule(mserSamplingInterval);
                }
            }
        }

        // ---------------------------------------------------------------------
        // Container for OBM confidence-interval results (Kotlin: CIResults).
        // The full implementation lives in PART 4+.
        // ---------------------------------------------------------------------

        private static final class CIResults {
            final Matrix QNCI;
            final Matrix UNCI;
            final Matrix RNCI;
            final Matrix TNCI;
            final Matrix ANCI;
            final Matrix WNCI;

            CIResults(Matrix QNCI, Matrix UNCI, Matrix RNCI, Matrix TNCI,
                      Matrix ANCI, Matrix WNCI) {
                this.QNCI = QNCI;
                this.UNCI = UNCI;
                this.RNCI = RNCI;
                this.TNCI = TNCI;
                this.ANCI = ANCI;
                this.WNCI = WNCI;
            }
        }

        // ---------------------------------------------------------------------
        // Forward-reference stubs for helpers translated in PART 4+ (Kotlin
        // lines 4500+).  The PART 3 code above calls into these.
        // ---------------------------------------------------------------------

        /**
         * Compute confidence intervals for all metrics using configured CI method.
         * Uses MSER truncation point and applies batch means to post-warmup observations.
         * @param confintLevel Confidence level (e.g., 0.95)
         * @return CIResults containing CI half-widths for all metrics
         */
        private CIResults computeOBMConfidenceIntervals(double confintLevel) {
            Matrix QNCI = new Matrix(numStations, numClasses);
            QNCI.fill(0.0);
            Matrix UNCI = new Matrix(numStations, numClasses);
            UNCI.fill(0.0);
            Matrix RNCI = new Matrix(numStations, numClasses);
            RNCI.fill(0.0);
            Matrix TNCI = new Matrix(numStations, numClasses);
            TNCI.fill(0.0);
            Matrix ANCI = new Matrix(numStations, numClasses);
            ANCI.fill(0.0);
            Matrix WNCI = new Matrix(numStations, numClasses);
            WNCI.fill(0.0);

            // Skip CI computation if disabled
            if ("none".equals(effectiveCiMethod)) {
                return new CIResults(QNCI, UNCI, RNCI, TNCI, ANCI, WNCI);
            }

            // Get truncation index from MSER
            int truncationIdx = mserTruncationBatch * effectiveMserBatchSize;

            // Compute CI for each service station and class
            for (int svcIdx = 0; svcIdx < serviceStations.size(); svcIdx++) {
                int serviceStation = serviceStations.get(svcIdx);
                for (int k = 0; k < numClasses; k++) {
                    // Queue length CI from observations
                    if (queueLengthObservations != null
                            && svcIdx < queueLengthObservations.length
                            && k < queueLengthObservations[svcIdx].length) {
                        List<Double> allObs = queueLengthObservations[svcIdx][k];
                        if (allObs != null && truncationIdx < allObs.size()) {
                            List<Double> postWarmupObs = allObs.subList(truncationIdx, allObs.size());
                            if (postWarmupObs.size() >= effectiveCiMinObs) {
                                int batchSize = Math.max(effectiveCiMinBatch,
                                        (int) Math.sqrt((double) postWarmupObs.size()));
                                StatTriple stats = computeCIStatistics(postWarmupObs, batchSize);
                                if (stats != null) {
                                    double tCrit = getTCriticalValueInternal(confintLevel, stats.df);
                                    QNCI.set(serviceStation, k, tCrit * stats.stdError);
                                }
                            }
                        }
                    }

                    // Utilization CI - derive from throughput and service rate
                    if (throughputObservations != null
                            && svcIdx < throughputObservations.length
                            && k < throughputObservations[svcIdx].length) {
                        List<Integer> compObs = throughputObservations[svcIdx][k];
                        if (compObs != null && truncationIdx < compObs.size()
                                && observationTimes != null && truncationIdx < observationTimes.size()) {
                            // Convert cumulative completions to per-interval rates
                            List<Integer> postWarmupComp = compObs.subList(truncationIdx, compObs.size());
                            List<Double> postWarmupTimes = observationTimes.subList(truncationIdx, observationTimes.size());
                            if (postWarmupComp.size() >= 2) {
                                ArrayList<Double> rates = new ArrayList<Double>();
                                for (int i = 1; i < postWarmupComp.size(); i++) {
                                    int deltaComp = postWarmupComp.get(i) - postWarmupComp.get(i - 1);
                                    double deltaTime = postWarmupTimes.get(i) - postWarmupTimes.get(i - 1);
                                    if (deltaTime > 0) {
                                        rates.add(((double) deltaComp) / deltaTime);
                                    }
                                }
                                if (rates.size() >= effectiveCiMinObs) {
                                    int batchSize = Math.max(effectiveCiMinBatch,
                                            (int) Math.sqrt((double) rates.size()));
                                    StatTriple stats = computeCIStatistics(rates, batchSize);
                                    if (stats != null) {
                                        double tCrit = getTCriticalValueInternal(confintLevel, stats.df);
                                        TNCI.set(serviceStation, k, tCrit * stats.stdError);
                                        // Utilization CI: convert throughput CI using service rate
                                        double mu = mus[svcIdx][k];
                                        if (mu > 0 && mu < Double.MAX_VALUE) {
                                            int nServers = (int) sn.nservers.get(serviceStations.get(svcIdx));
                                            UNCI.set(serviceStation, k,
                                                    (tCrit * stats.stdError) / (mu * nServers));
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Response time CI from tally variance
                    if (responseTimeTally != null
                            && svcIdx < responseTimeTally.length
                            && k < responseTimeTally[svcIdx].length
                            && responseTimeTally[svcIdx][k] != null) {
                        Tally tally = responseTimeTally[svcIdx][k];
                        if (tally.numberObs() > 30) {
                            double stdErr = tally.standardDeviation() / Math.sqrt((double) tally.numberObs());
                            double tCrit = getTCriticalValueInternal(confintLevel, tally.numberObs() - 1);
                            RNCI.set(serviceStation, k, tCrit * stdErr);
                        }
                    }
                }
            }

            // WNCI = RNCI (residence time same as response time for single-visit)
            for (int i = 0; i < numStations; i++) {
                for (int k = 0; k < numClasses; k++) {
                    WNCI.set(i, k, RNCI.get(i, k));
                }
            }

            return new CIResults(QNCI, UNCI, RNCI, TNCI, ANCI, WNCI);
        }

        /**
         * Get t-distribution critical value for given confidence level and degrees of freedom.
         * Internal version using confintLevel directly (not alpha).
         * @param confintLevel Confidence level (e.g., 0.95 for 95% CI)
         * @param df Degrees of freedom
         * @return Critical value from t-distribution
         */
        private double getTCriticalValueInternal(double confintLevel, int df) {
            // Lookup based on confidence level directly for better precision
            double[] tTable95 = new double[]{
                12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
                2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
                2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042
            };
            double[] tTable90 = new double[]{
                6.314, 2.920, 2.353, 2.132, 2.015, 1.943, 1.895, 1.860, 1.833, 1.812,
                1.796, 1.782, 1.771, 1.761, 1.753, 1.746, 1.740, 1.734, 1.729, 1.725,
                1.721, 1.717, 1.714, 1.711, 1.708, 1.706, 1.703, 1.701, 1.699, 1.697
            };
            double[] tTable99 = new double[]{
                63.657, 9.925, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.250, 3.169,
                3.106, 3.055, 3.012, 2.977, 2.947, 2.921, 2.898, 2.878, 2.861, 2.845,
                2.831, 2.819, 2.807, 2.797, 2.787, 2.779, 2.771, 2.763, 2.756, 2.750
            };

            // Table index is df - 1 (df=1 is at index 0)
            int tableIdx = Math.min(df, 30) - 1;
            if (tableIdx < 0) return 1.96; // Fallback to z-value

            // Select table based on confidence level (checking higher confidence first)
            if (confintLevel >= 0.99) {
                return tableIdx < tTable99.length ? tTable99[tableIdx] : 2.576;
            } else if (confintLevel >= 0.95) {
                return tableIdx < tTable95.length ? tTable95[tableIdx] : 1.96;
            } else if (confintLevel >= 0.90) {
                return tableIdx < tTable90.length ? tTable90[tableIdx] : 1.645;
            } else {
                return 1.96; // Default to 95% CI z-value
            }
        }

        /** Average queue length (per-class). Kotlin line 13038. */
        private double getAvgQueueLength(int queueIdx, int classId) {
            // Use MSER-5 truncated cumulative queue time if available. The
            // estimate must stay time-weighted: MSER samples are event-spaced,
            // so an unweighted mean of the per-interval averages overweights
            // congested epochs (many short intervals) and is biased upward on
            // bursty workloads. Mirror the utilization estimator instead:
            // (queue-time integral - integral at truncation) / elapsed time.
            if (mserEnabled && queueTimeObservations != null
                    && queueIdx < queueTimeObservations.length
                    && classId < queueTimeObservations[queueIdx].length) {
                List<Double> qtObs = queueTimeObservations[queueIdx][classId];
                int truncationIdx = mserTruncationBatch * effectiveMserBatchSize;
                if (qtObs != null && truncationIdx < qtObs.size()
                        && observationTimes != null
                        && truncationIdx < observationTimes.size()) {
                    double startQueueTime = qtObs.get(truncationIdx);
                    double endQueueTime = totalQueueTime[queueIdx][classId];
                    double startTime = observationTimes.get(truncationIdx);
                    double elapsed = ssjSim.time() - startTime;
                    if (elapsed > 0) {
                        return (endQueueTime - startQueueTime) / elapsed;
                    }
                }
            }
            // Fallback to time-weighted average
            double simTime = getActualSimTime();
            if (simTime > 0) {
                return totalQueueTime[queueIdx][classId] / simTime;
            }
            return 0.0;
        }

        /** Per-class per-server utilization. Kotlin line 13126. */
        private double getUtilization(int queueIdx, int classId) {
            // For Delay nodes (infinite servers), utilization = traffic intensity = λ/μ
            if (isDelayNode.get(queueIdx)) {
                double mu = mus[queueIdx][classId];
                if (mu > 0 && mu < Double.MAX_VALUE) {
                    return getThroughput(queueIdx, classId) / mu;
                }
                return getAvgQueueLength(queueIdx, classId);
            }
            // Load-dependent / class-dependent stations: the scaling multiplies the
            // nominal service rate, so the station's capacity is max_n scaling(n) in
            // units of that nominal rate, not numServers. The measured busy time
            // counts the single server as busy whenever n >= 1, which yields
            // P(busy), not the fraction of capacity in use: a scaling emulating c
            // servers then reports ~1 at saturation instead of E[busy]/c, and
            // disagrees with CTMC/NC on the same model. Use the utilization law
            // against the measured throughput and normalize by the peak scaling,
            // matching solver_ncld (max(lldscaling(ist,:))) and solver_nc_conv. The
            // delay branch above likewise prefers T/mu over the busy-time tally.
            double peakScaling = getPeakScaling(queueIdx);
            if (peakScaling > 0) {
                double muLd = mus[queueIdx][classId];
                if (muLd > 0 && muLd < Double.MAX_VALUE) {
                    return getThroughput(queueIdx, classId) / (muLd * peakScaling);
                }
            }
            // Per-server Utilization = (busy time + blocking time) / (simTime * compatibleServers)
            int compatibleServers = getCompatibleServerCount(queueIdx, classId);
            double blockingTime = 0.0;
            if (totalBlockingTime != null
                    && queueIdx < totalBlockingTime.length
                    && classId < totalBlockingTime[queueIdx].length) {
                blockingTime = totalBlockingTime[queueIdx][classId];
            }
            // Use MSER-5 truncated busy time if available
            if (mserEnabled && busyTimeObservations != null && observationTimes != null
                    && queueIdx < busyTimeObservations.length
                    && classId < busyTimeObservations[queueIdx].length) {
                List<Double> busyObs = busyTimeObservations[queueIdx][classId];
                int truncationIdx = mserTruncationBatch * effectiveMserBatchSize;
                if (busyObs != null && truncationIdx < busyObs.size()
                        && truncationIdx < observationTimes.size()) {
                    double startBusy = busyObs.get(truncationIdx);
                    double endBusy = totalBusyTime[queueIdx][classId];
                    double startTime = observationTimes.get(truncationIdx);
                    double endTime = ssjSim.time();
                    double elapsed = endTime - startTime;
                    if (elapsed > 0 && compatibleServers > 0) {
                        double startBlocking = 0.0;
                        if (blockingTimeObservations != null
                                && queueIdx < blockingTimeObservations.length
                                && classId < blockingTimeObservations[queueIdx].length
                                && blockingTimeObservations[queueIdx][classId] != null
                                && truncationIdx < blockingTimeObservations[queueIdx][classId].size()) {
                            startBlocking = blockingTimeObservations[queueIdx][classId].get(truncationIdx);
                        }
                        return ((endBusy - startBusy) + (blockingTime - startBlocking))
                                / (elapsed * compatibleServers);
                    }
                }
            }
            // Fallback
            double simTime = getActualSimTime();
            if (simTime > 0 && compatibleServers > 0) {
                return (totalBusyTime[queueIdx][classId] + blockingTime)
                        / (simTime * compatibleServers);
            }
            return 0.0;
        }

        /** Average response time. Kotlin line 13183. */
        private double getAvgResponseTime(int queueIdx, int classId) {
            if (responseTimeTally[queueIdx][classId].numberObs() > 0) {
                return responseTimeTally[queueIdx][classId].average();
            }
            if (slotted) {
                // The M/M/1 fallback below assumes exponential interarrival and
                // service times, which a slotted model does not have. Reporting it
                // would silently substitute a continuous-time approximation for a
                // measurement that simply was not taken.
                return 0.0;
            }
            // Use M/M/1 formula as fallback
            double mu = mus[queueIdx][classId];
            double lambda = getThroughput(queueIdx, classId);
            if (mu > lambda && mu < Double.MAX_VALUE) {
                return 1.0 / (mu - lambda);
            }
            return 0.0;
        }

        /** Per-class throughput at queueIdx. Kotlin line 13194. */
        private double getThroughput(int queueIdx, int classId) {
            // Use MSER-5 truncated throughput if available
            if (mserEnabled && throughputObservations != null && observationTimes != null
                    && queueIdx < throughputObservations.length
                    && classId < throughputObservations[queueIdx].length) {
                List<Integer> tputObs = throughputObservations[queueIdx][classId];
                int truncationIdx = mserTruncationBatch * effectiveMserBatchSize;
                if (tputObs != null && truncationIdx < tputObs.size()
                        && truncationIdx < observationTimes.size()) {
                    int startCompletions = tputObs.get(truncationIdx);
                    int endCompletions = completedCustomers[queueIdx][classId];
                    double startTime = observationTimes.get(truncationIdx);
                    double endTime = ssjSim.time();
                    double elapsed = endTime - startTime;
                    if (elapsed > 0) {
                        return ((double) (endCompletions - startCompletions)) / elapsed;
                    }
                }
            }
            // Fallback
            double simTime = getActualSimTime();
            if (simTime > 0) {
                return ((double) completedCustomers[queueIdx][classId]) / simTime;
            }
            return 0.0;
        }

        /** System throughput per class. Kotlin line 13234. */
        private double getSystemThroughput(int classId) {
            double simTime = getActualSimTime();
            if (simTime > 0) {
                return ((double) systemCompletedCustomers[classId]) / simTime;
            }
            return 0.0;
        }

        /** System response time per class. Kotlin line 13244. */
        private double getSystemResponseTime(int classId) {
            if (systemResponseTimeTally[classId].numberObs() > 0) {
                return systemResponseTimeTally[classId].average();
            }
            return 0.0;
        }

        /** Safety bound on immediate (vanishing-marking) firings per epoch. */
        private static final int IMMEDIATE_FIRING_GUARD = 1000000;

        /**
         * Advance the Petri-net marking using atomic-firing race semantics,
         * matching the GSPN continuous-time Markov chain that SolverJMT/SolverCTMC
         * realise. Tokens are <em>not</em> reserved at enabling time; they remain in
         * their input place during the firing delay and are consumed atomically at
         * the firing epoch. This method is invoked after every marking change.
         *
         * <p>Phase 1 resolves the vanishing marking: while any immediate transition
         * is enabled, exactly one is selected (highest firing priority, ties broken
         * by weighted-random draw) and fired atomically in zero simulated time.
         *
         * <p>Phase 2 handles tangible timed transitions: for each timed mode the
         * number of in-flight firings is topped up to {@code min(enablingDegree,
         * servers)} by sampling a firing delay and scheduling a {@link TransitionFiring}
         * event. Single-server modes therefore keep at most one clock running, while
         * infinite/multi-server modes run one clock per concurrently enabled binding,
         * reproducing the enabling-degree rate scaling of the CTMC.
         */
        private void checkAndFireTransitions() {
            // A marking-dependent firing rate is stale once the marking moved, so
            // cancel and free every in-flight dependent clock; the timed top-up
            // below redraws it at the current marking. Exact for exponential firing.
            resampleDependentFirings();
            // Phase 1: fire enabled immediate transitions (vanishing markings), one
            // at a time, until none remain enabled.
            int guard = 0;
            while (true) {
                TransitionModeInfo bestMode = null;
                int bestTrans = -1;
                int maxPriority = Integer.MIN_VALUE;
                double tieWeightTotal = 0.0;
                int tieCount = 0;
                for (int transListIdx = 0; transListIdx < transitionNodes.size(); transListIdx++) {
                    List<TransitionModeInfo> modes = transitionModes[transListIdx];
                    for (int mi = 0; mi < modes.size(); mi++) {
                        TransitionModeInfo m = modes.get(mi);
                        if (m.timingStrategy != TimingStrategy.IMMEDIATE) continue;
                        if (transitionModeEnablingDegree(m) < 1) continue;
                        if (m.priority > maxPriority) {
                            maxPriority = m.priority;
                            bestMode = m;
                            bestTrans = transListIdx;
                            tieWeightTotal = m.weight;
                            tieCount = 1;
                        } else if (m.priority == maxPriority) {
                            // Reservoir-sample among equal-priority modes by weight so
                            // the winning binding is chosen with probability
                            // weight / sum(weights) over the conflict set.
                            tieWeightTotal += m.weight;
                            tieCount++;
                            if (tieWeightTotal > 0
                                    && siroRng.nextDouble() < (m.weight / tieWeightTotal)) {
                                bestMode = m;
                                bestTrans = transListIdx;
                            }
                        }
                    }
                }
                if (bestMode == null) break;
                fireAtomic(bestTrans, bestMode);
                if (++guard > IMMEDIATE_FIRING_GUARD) break;
            }

            // Phase 2: schedule timed firings via a global race that reserves the
            // shared input tokens. committed[p][c] holds the tokens already promised
            // to in-flight firings; the enabling degree of each candidate binding is
            // evaluated against the uncommitted tokens (placeTokens - committed), so
            // competing modes of one transition and different transitions drawing on
            // the same place never over-commit it. Physical place tokens are not
            // touched here - they are consumed atomically at the firing epoch - so
            // the reported place queue length still counts every token that has not
            // yet fired, matching SolverCTMC/SolverJMT. Whichever enabled binding has
            // the earliest sampled firing time wins each slot, reproducing the
            // race/enabling-degree rate scaling of the CTMC.
            int numP = placeNodes.size();
            for (int t = 0; t < transitionNodes.size(); t++) {
                List<TransitionModeInfo> modes = transitionModes[t];
                boolean[] processed = new boolean[modes.size()];
                for (int gi = 0; gi < modes.size(); gi++) {
                    if (processed[gi]) continue;
                    TransitionModeInfo gm = modes.get(gi);
                    if (gm.timingStrategy == TimingStrategy.IMMEDIATE) {
                        processed[gi] = true;
                        continue;
                    }

                    // Group timed modes that share identical enabling conditions.
                    // Such modes are redundant alternatives competing for exactly the
                    // same input tokens, so they must race for a shared firing slot
                    // (e.g. three service-time modes of one self-looping transition).
                    // Modes with distinct enabling conditions draw on different token
                    // sets and are scheduled independently, so they are handled in
                    // their own group.
                    List<Integer> group = new ArrayList<Integer>();
                    for (int mj = gi; mj < modes.size(); mj++) {
                        if (processed[mj]) continue;
                        TransitionModeInfo mm = modes.get(mj);
                        if (mm.timingStrategy == TimingStrategy.IMMEDIATE) continue;
                        if (sameEnablingConditions(gm, mm)) {
                            group.add(mj);
                            processed[mj] = true;
                        }
                    }

                    // committed counts tokens already promised to this group's
                    // in-flight firings; physical place tokens are consumed only at
                    // the firing epoch, so the reported place queue length still
                    // counts every token that has not yet fired (matching
                    // SolverCTMC/SolverJMT). Sample an independent firing time for
                    // each candidate binding (one per free server up to the available
                    // enabling degree), then commit in increasing time order: a
                    // binding whose shared tokens remain uncommitted starts and keeps
                    // its own clock, while one whose tokens were already claimed by an
                    // earlier-firing redundant mode is dropped (it lost the race).
                    // Cross-group and cross-transition token contention is resolved by
                    // the enabling re-check at the firing epoch.
                    int[][] committed = new int[numP][numClasses];
                    for (int idx = 0; idx < group.size(); idx++) {
                        int mi = group.get(idx);
                        int inflight = transitionInService[t][mi];
                        if (inflight <= 0) continue;
                        int[][] enab = modes.get(mi).enablingConditions;
                        for (int p = 0; p < numP; p++) {
                            for (int c = 0; c < numClasses; c++) {
                                committed[p][c] += inflight * enab[p][c];
                            }
                        }
                    }

                    List<double[]> candidates = new ArrayList<double[]>();
                    for (int idx = 0; idx < group.size(); idx++) {
                        int mi = group.get(idx);
                        TransitionModeInfo m = modes.get(mi);
                        long freeServers = (long) m.numServers - transitionInService[t][mi];
                        if (freeServers <= 0) continue;
                        int degree = availableEnablingDegree(m, committed);
                        long nInstances = Math.min(freeServers, (long) degree);
                        for (long s = 0; s < nInstances; s++) {
                            double tm = sampleTransitionFiringTime(t, mi);
                            candidates.add(new double[]{tm, mi});
                        }
                    }
                    candidates.sort(new java.util.Comparator<double[]>() {
                        @Override
                        public int compare(double[] a, double[] b) {
                            return Double.compare(a[0], b[0]);
                        }
                    });
                    for (int i = 0; i < candidates.size(); i++) {
                        double[] cand = candidates.get(i);
                        int mi = (int) cand[1];
                        TransitionModeInfo m = modes.get(mi);
                        if (transitionInService[t][mi] >= m.numServers) continue;
                        if (!isTimedModeSchedulable(m, committed)) continue;
                        transitionInService[t][mi]++;
                        int[][] enab = m.enablingConditions;
                        for (int p = 0; p < numP; p++) {
                            for (int c = 0; c < numClasses; c++) {
                                committed[p][c] += enab[p][c];
                            }
                        }
                        TransitionFiring tf = new TransitionFiring(t, mi);
                        tf.schedule(cand[0]);
                        // Track dependent clocks so they can be resampled on a marking change.
                        if (m.firingDep != null) inflightDependent.add(tf);
                    }
                }
            }
        }

        /** Whether two transition modes have identical enabling (input arc) conditions. */
        private boolean sameEnablingConditions(TransitionModeInfo a, TransitionModeInfo b) {
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    if (a.enablingConditions[placeListIdx][classIdx]
                            != b.enablingConditions[placeListIdx][classIdx]) {
                        return false;
                    }
                }
            }
            return true;
        }

        /**
         * Enabling degree of a timed mode against the tokens still available after
         * subtracting those already committed to in-flight firings of the same
         * transition. Returns 0 when inhibited or when the mode has no input arc.
         */
        private int availableEnablingDegree(TransitionModeInfo mode, int[][] committed) {
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int maxTokens = mode.inhibitingConditions[placeListIdx][classIdx];
                    if (maxTokens < Integer.MAX_VALUE
                            && placeAvail(placeListIdx, classIdx) >= maxTokens) {
                        return 0;
                    }
                }
            }
            int degree = Integer.MAX_VALUE;
            boolean hasInputArc = false;
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int required = mode.enablingConditions[placeListIdx][classIdx];
                    if (required > 0) {
                        hasInputArc = true;
                        int avail = placeAvail(placeListIdx, classIdx)
                                - committed[placeListIdx][classIdx];
                        int d = avail / required;
                        if (d < degree) degree = d;
                    }
                }
            }
            if (!hasInputArc) return 0;
            return (degree == Integer.MAX_VALUE) ? 0 : degree;
        }

        /**
         * Whether a timed transition mode can start another firing given the tokens
         * already committed to in-flight firings: no inhibiting threshold is reached
         * and every input arc is covered by the uncommitted tokens
         * {@code placeTokens - committed}. Returns false for a mode with no input arc.
         */
        private boolean isTimedModeSchedulable(TransitionModeInfo mode, int[][] committed) {
            boolean hasInputArc = false;
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int maxTokens = mode.inhibitingConditions[placeListIdx][classIdx];
                    if (maxTokens < Integer.MAX_VALUE
                            && placeAvail(placeListIdx, classIdx) >= maxTokens) {
                        return false;
                    }
                    int required = mode.enablingConditions[placeListIdx][classIdx];
                    if (required > 0) {
                        hasInputArc = true;
                        if (placeAvail(placeListIdx, classIdx)
                                - committed[placeListIdx][classIdx] < required) {
                            return false;
                        }
                    }
                }
            }
            return hasInputArc;
        }

        /**
         * Enabling degree of a transition mode in the current marking: the largest
         * number of simultaneous bindings, i.e. the floor over input places of
         * available tokens divided by the required arc weight. Returns 0 when an
         * inhibiting condition is met or the mode has no input arc.
         */
        private int transitionModeEnablingDegree(TransitionModeInfo mode) {
            // Inhibiting arcs: disabled if any inhibiting place has reached its
            // threshold token count.
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int maxTokens = mode.inhibitingConditions[placeListIdx][classIdx];
                    if (maxTokens < Integer.MAX_VALUE
                            && placeAvail(placeListIdx, classIdx) >= maxTokens) {
                        return 0;
                    }
                }
            }
            int degree = Integer.MAX_VALUE;
            boolean hasInputArc = false;
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int required = mode.enablingConditions[placeListIdx][classIdx];
                    if (required > 0) {
                        hasInputArc = true;
                        int d = placeAvail(placeListIdx, classIdx) / required;
                        if (d < degree) degree = d;
                    }
                }
            }
            if (!hasInputArc) return 0;
            return (degree == Integer.MAX_VALUE) ? 0 : degree;
        }

        /**
         * Fire a transition mode atomically: consume the required input tokens from
         * the marking and immediately produce the firing outcome at the output nodes.
         * Used for both immediate transitions (Phase 1) and the completion of timed
         * firings.
         */
        private void fireAtomic(int transListIdx, TransitionModeInfo mode) {
            double currentTime = ssjSim.time();
            int transNodeIdx = transitionNodes.get(transListIdx);
            trackEvent();

            // Consume input tokens (updating time-weighted place statistics first).
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                double elapsed = currentTime - lastPlaceUpdateTime[placeListIdx];
                if (elapsed > 0) {
                    for (int k = 0; k < numClasses; k++) {
                        totalPlaceTokenTime[placeListIdx][k] +=
                                placeTokens[placeListIdx][k] * elapsed;
                    }
                    lastPlaceUpdateTime[placeListIdx] = currentTime;
                }
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int required = mode.enablingConditions[placeListIdx][classIdx];
                    if (required > 0) {
                        consumePlaceTokens(placeListIdx, classIdx, required);
                        placeCompletions[placeListIdx][classIdx] += required;
                    }
                }
            }

            checkEventCountStop();

            // Produce output tokens at the firing outcome nodes.
            for (int nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int produced = mode.firingOutcomes[nodeIdx][classIdx];
                    for (int t = 0; t < produced; t++) {
                        routeTokenFromTransition(transNodeIdx, nodeIdx, classIdx);
                    }
                }
            }
        }

        /**
         * Updates time-weighted orbit statistics before changing orbit size.
         */
        private void updateOrbitTimeStats(int queueIdx, int classId) {
            double currentTime = ssjSim.time();
            double elapsed = currentTime - lastOrbitUpdateTime[queueIdx][classId];
            if (elapsed > 0 && warmupDone) {
                totalOrbitTime[queueIdx][classId] += currentOrbitSize[queueIdx][classId] * elapsed;
            }
            lastOrbitUpdateTime[queueIdx][classId] = currentTime;
        }

        /** Body translated in PART 8 (Kotlin line 11149). */
        private void updateBusyStats(int svcIdx, int classId) {
            updateBusyStatsImpl(svcIdx, classId);
        }

        /** Updates time-weighted queue stats. Kotlin line 12645. */
        private void updateQueueStats(int queueIdx, int classId) {
            // Flush the Markov reward integral for the joint state that held over
            // the interval ending now, before this cell's queue length changes.
            updateRewardStats();
            double currentTime = ssjSim.time();
            double elapsed = currentTime - lastQueueUpdateTime[queueIdx][classId];
            if (elapsed > 0) {
                // Effective queue length = waiting/in-service + sync-blocked
                //                          + BAS-blocked at dest + BBS-blocked at dest
                int effectiveQueueLength = currentQueueLength[queueIdx][classId]
                        + currentBlockedServers[queueIdx][classId]
                        + basBlockedAtDest[queueIdx][classId]
                        + bbsBlockedAtDest[queueIdx][classId];
                totalQueueTime[queueIdx][classId] += effectiveQueueLength * elapsed;
                lastQueueUpdateTime[queueIdx][classId] = currentTime;
            }
        }

        /**
         * Builds the aggregated joint state row in the CTMC stateSpaceAggr layout
         * (station-major, class-minor) from the live integer queue lengths, so a
         * reward function registered via model.setReward evaluates identically under
         * LDES and CTMC.
         */
        private Matrix currentRewardRow() {
            Matrix row = new Matrix(1, rewardRowCols);
            row.zero();
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int station = serviceStations.get(svcIdx);
                for (int k = 0; k < numClasses; k++) {
                    int njk = currentQueueLength[svcIdx][k]
                            + currentBlockedServers[svcIdx][k]
                            + basBlockedAtDest[svcIdx][k]
                            + bbsBlockedAtDest[svcIdx][k];
                    if (njk != 0) {
                        row.set(0, station * numClasses + k, njk);
                    }
                }
            }
            return row;
        }

        /**
         * Time-integrates the Markov reward over the interval since the last flush,
         * using the exact integer joint state. Repeated calls within a single event
         * (identical ssjSim.time()) are no-ops. Results feed getAvgReward/getTranReward.
         */
        private void updateRewardStats() {
            if (!hasReward && !exportStateHistogram) {
                return;
            }
            double currentTime = ssjSim.time();
            double dt = currentTime - rewardLastUpdateTime;
            if (dt <= 0) {
                return;
            }
            Matrix row = currentRewardRow();
            // Memoize by joint-state key: small state spaces recur constantly.
            StringBuilder sb = new StringBuilder(rewardRowCols * 3);
            for (int c = 0; c < rewardRowCols; c++) {
                sb.append((int) row.get(0, c)).append(',');
            }
            String key = sb.toString();

            if (hasReward) {
                double[] rvals = rewardStateCache.get(key);
                if (rvals == null) {
                    rvals = new double[rewardFnArr.length];
                    for (int ri = 0; ri < rewardFnArr.length; ri++) {
                        rvals[ri] = rewardFnArr[ri].compute(row, sn);
                    }
                    rewardStateCache.put(key, rvals);
                }
                for (int ri = 0; ri < rewardArea.length; ri++) {
                    rewardArea[ri] += rvals[ri] * dt;
                }
                // Record the transient step value r(X(t)) at the interval start time.
                double[] entry = new double[rewardFnArr.length + 1];
                entry[0] = rewardLastUpdateTime;
                for (int ri = 0; ri < rewardFnArr.length; ri++) {
                    entry[ri + 1] = rvals[ri];
                }
                rewardTranSeries.add(entry);
            }

            if (exportStateHistogram) {
                // Accumulate the exact joint-state residence time so a host language
                // can evaluate its own reward functions on the empirical distribution.
                double[] hentry = stateHistogram.get(key);
                if (hentry == null) {
                    hentry = new double[rewardRowCols + 1];
                    for (int c = 0; c < rewardRowCols; c++) {
                        hentry[c + 1] = row.get(0, c);
                    }
                    stateHistogram.put(key, hentry);
                }
                hentry[0] += dt;
                // Record the integer joint state at the interval start time for the
                // transient reward trajectory r(X(t)).
                double[] tentry = new double[rewardRowCols + 1];
                tentry[0] = rewardLastUpdateTime;
                for (int c = 0; c < rewardRowCols; c++) {
                    tentry[c + 1] = row.get(0, c);
                }
                stateTranSeries.add(tentry);
            }

            rewardTotalTime += dt;
            rewardLastUpdateTime = currentTime;
        }

        /** Body translated in PART 8 (Kotlin line 11195). */
        private void updatePSBusyStats(int svcIdx) {
            updatePSBusyStatsImpl(svcIdx);
        }

        // closeLoggers() — translated in PART 6 (Kotlin lines 6229-6239).

        /** Body translated in PART 8 (Kotlin line 11106). */
        private void resetStatistics() {
            resetStatisticsImpl();
        }

        /** End-of-simulation event used by transient mode to stop {@link Simulator} at
         *  the configured time horizon.  Body translated in PART 8 (Kotlin line 11235). */
        private final class EndOfSimulation extends SimEvent {
            @Override
            protected double slotPhase() {
                return SLOT_PHASE_BOOKKEEPING;
            }

            @Override
            public void actions() {
                endOfSimulationActions();
            }
        }

    // ===================================================================
    // PART 4 — Kotlin lines 4501-6000
    // ===================================================================

        /**
         * Compute MSER-5 truncation point for a list of observations.
         *
         * <p>MSER-5 uses batches of size 5 and selects the truncation point d
         * that minimizes variance(Z_{d+1..N}) / (N - d)^2 where Z are batch means.
         * Returns 0 if there are not enough observations to perform MSER-5.</p>
         */
        private int computeMSER5TruncationPoint(java.util.List<Double> observations) {
            int n = observations.size();
            if (n < effectiveMserBatchSize * 4) {
                return 0;
            }
            int numBatches = n / effectiveMserBatchSize;
            double[] batchMeans = new double[numBatches];
            for (int j = 0; j < numBatches; j++) {
                double sum = 0.0;
                for (int i = 0; i < effectiveMserBatchSize; i++) {
                    sum += observations.get(j * effectiveMserBatchSize + i);
                }
                batchMeans[j] = sum / effectiveMserBatchSize;
            }
            double minMSER = Double.MAX_VALUE;
            int optimalD = 0;
            int maxD = numBatches / 2;
            for (int d = 0; d < maxD; d++) {
                int remainingBatches = numBatches - d;
                if (remainingBatches < 2) break;
                double sum = 0.0;
                for (int j = d; j < numBatches; j++) {
                    sum += batchMeans[j];
                }
                double mean = sum / remainingBatches;
                double variance = 0.0;
                for (int j = d; j < numBatches; j++) {
                    double diff = batchMeans[j] - mean;
                    variance += diff * diff;
                }
                variance /= (remainingBatches - 1);
                double mser = variance / ((double) remainingBatches * remainingBatches);
                if (mser < minMSER) {
                    minMSER = mser;
                    optimalD = d;
                }
            }
            return optimalD;
        }

        /**
         * MSER-5 truncation post-processing (Kotlin lines 4517-4566).
         */
        private void applyMSER5Truncation() {
            if (!mserEnabled || queueLengthObservations == null) return;
            // Bail if no observations recorded
            boolean anyObs = false;
            for (int q = 0; q < numServiceNodes && !anyObs; q++) {
                for (int k = 0; k < numClasses && !anyObs; k++) {
                    if (queueLengthObservations[q][k] != null
                            && !queueLengthObservations[q][k].isEmpty()) {
                        anyObs = true;
                    }
                }
            }
            if (!anyObs) return;

            // First try aggregate queue length (works well for open networks)
            List<Double> aggregateObs = new ArrayList<Double>();
            int obsCount = (observationTimes != null) ? observationTimes.size() : 0;
            for (int i = 0; i < obsCount; i++) {
                double total = 0.0;
                for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                    for (int k = 0; k < numClasses; k++) {
                        List<Double> obs = queueLengthObservations[qIdx][k];
                        if (obs != null && i < obs.size()) {
                            total += obs.get(i);
                        }
                    }
                }
                aggregateObs.add(total);
            }

            int aggTruncation = 0;
            if (aggregateObs.size() >= effectiveMserBatchSize * 4) {
                aggTruncation = computeMSER5TruncationPoint(aggregateObs);
            }

            // If aggregate finds no truncation (e.g., closed network), fall back per queue
            if (aggTruncation == 0) {
                for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                    for (int k = 0; k < numClasses; k++) {
                        List<Double> obs = queueLengthObservations[qIdx][k];
                        if (obs != null && obs.size() >= effectiveMserBatchSize * 4) {
                            int trunc = computeMSER5TruncationPoint(obs);
                            if (trunc > aggTruncation) {
                                aggTruncation = trunc;
                            }
                        }
                    }
                }
            }

            mserTruncationBatch = aggTruncation;
            int truncationObsIdx = mserTruncationBatch * effectiveMserBatchSize;

            if (truncationObsIdx > 0 && observationTimes != null
                    && truncationObsIdx < observationTimes.size()) {
                warmupEndTime = observationTimes.get(truncationObsIdx);
                if (options.verbose == VerboseLevel.DEBUG) {
                    System.out.println("MSER-5: Truncation at batch "
                            + mserTruncationBatch + " (t=" + warmupEndTime + ")");
                }
            } else {
                warmupEndTime = 0.0;
            }
        }

        @SuppressWarnings("unchecked")
        private void initializeGenerators() {
            routingRng = new MRG32k3a();
            pasRng = new MRG32k3a();
            if (seed > 0) {
                routingRng.setSeed(new long[] { seed, seed + 1, seed + 2,
                        seed + 3, seed + 4, seed + 5 });
                pasRng.setSeed(new long[] { seed + 7000, seed + 7001, seed + 7002,
                        seed + 7003, seed + 7004, seed + 7005 });
                siroRng = new Random(seed + 99999);
            } else {
                siroRng = new Random();
            }

            // Determine arrival process types from sn data structure
            arrivalProcessType = new ProcessType[numSources][numClasses];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                int istStation = sourceStations.get(srcIdx);
                jline.lang.nodes.Station station = sn.stations.get(istStation);
                for (int k = 0; k < numClasses; k++) {
                    jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                    Map<jline.lang.JobClass, ProcessType> rmap =
                            (sn.procid != null) ? sn.procid.get(station) : null;
                    ProcessType pt = (rmap != null) ? rmap.get(jobClass) : null;
                    arrivalProcessType[srcIdx][k] = (pt != null) ? pt : ProcessType.DISABLED;
                }
            }

            // Initialize PH process matrices for arrivals from sn.proc
            arrivalProc = new MatrixCell[numSources][numClasses];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                int istStation = sourceStations.get(srcIdx);
                jline.lang.nodes.Station station = sn.stations.get(istStation);
                for (int k = 0; k < numClasses; k++) {
                    jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                    Map<jline.lang.JobClass, MatrixCell> rmap =
                            (sn.proc != null) ? sn.proc.get(station) : null;
                    arrivalProc[srcIdx][k] = (rmap != null) ? rmap.get(jobClass) : null;
                }
            }

            // Initialize random generators for PH/MAP/BMAP/MMAP arrival sampling
            arrivalRng = new Random[numSources][numClasses];
            arrivalMapSampler = new Map_sample.MapSampler[numSources][numClasses];
            arrivalMeSampler = new Me_sample.MeSampler[numSources][numClasses];
            arrivalRapSampler = new Rap_sample.RapSampler[numSources][numClasses];
            arrivalBmapSampler = new Map_sample.BmapSampler[numSources][numClasses];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    ProcessType procType = arrivalProcessType[srcIdx][k];
                    if (procType == ProcessType.PH || procType == ProcessType.APH
                            || procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN
                            || procType == ProcessType.COX2 || procType == ProcessType.MAP
                            || procType == ProcessType.MMPP2 || procType == ProcessType.BMAP
                            || procType == ProcessType.MMAP || procType == ProcessType.ME
                            || procType == ProcessType.RAP
                            || procType == ProcessType.NHPP) {
                        if (seed > 0) {
                            long offset = ((long) srcIdx * numClasses + k) * 10 + 2000;
                            arrivalRng[srcIdx][k] = new Random(seed + offset);
                        } else {
                            arrivalRng[srcIdx][k] = new Random();
                        }
                    } else {
                        arrivalRng[srcIdx][k] = null;
                    }
                }
            }

            // Initialize batch size cache for BMAP arrivals (default batch size = 1)
            arrivalBatchSize = new int[numSources][numClasses];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    arrivalBatchSize[srcIdx][k] = 1;
                }
            }

            // Explicit batch-size law (sn.arrivalbatch), independent of BMAP:
            // the interarrival distribution spaces the epochs and this decides
            // how many jobs each epoch releases. Seeded off the run seed so the
            // batch stream is reproducible and independent of the interarrival
            // stream.
            arrivalBatchDist = new DiscreteDistribution[numSources][numClasses];
            arrivalBatchRng = new Random[numSources][numClasses];
            if (sn.arrivalbatch != null) {
                for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                    for (int k = 0; k < numClasses; k++) {
                        if (k < sn.arrivalbatch.size() && sn.arrivalbatch.get(k) != null) {
                            arrivalBatchDist[srcIdx][k] = sn.arrivalbatch.get(k);
                            arrivalBatchRng[srcIdx][k] =
                                    new Random(seed + 991L * (srcIdx + 1L) + 97L * (k + 1L));
                        }
                    }
                }
            }

            // Initialize MMAP marked-arrival state: shared phase-carrying
            // samplers and the mark->class binding from sn.markidx. The first
            // marked class of a source acts as the carrier: only its
            // ExternalArrival events are scheduled, and each event routes the
            // arriving job to the class bound to the sampled mark.
            arrivalMmapSampler = new jline.api.mam.Mmap_sample.MmapSampler[numSources][numClasses];
            arrivalPendingMark = new int[numSources][numClasses];
            markedGroupClass = new int[numSources][];
            markedCarrierClass = new int[numSources];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                markedCarrierClass[srcIdx] = -1;
                if (sn.markidx == null) {
                    continue;
                }
                int istStation = sourceStations.get(srcIdx);
                int maxMark = 0;
                for (int k = 0; k < numClasses; k++) {
                    int m = (int) sn.markidx.get(istStation, k);
                    if (m > maxMark) {
                        maxMark = m;
                    }
                }
                if (maxMark > 0) {
                    int[] group = new int[maxMark + 1];
                    java.util.Arrays.fill(group, -1);
                    for (int k = 0; k < numClasses; k++) {
                        int m = (int) sn.markidx.get(istStation, k);
                        if (m > 0) {
                            group[m] = k;
                        }
                    }
                    markedGroupClass[srcIdx] = group;
                    markedCarrierClass[srcIdx] = group[1];
                }
            }

            // Initialize arrival rate schedules (NHPP) from arrivalProc
            arrivalSchedule = new double[numSources][numClasses][][];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    ProcessType pt = arrivalProcessType[srcIdx][k];
                    if (pt == ProcessType.NHPP) {
                        arrivalSchedule[srcIdx][k] =
                                unpackRateSchedule(arrivalProc[srcIdx][k]);
                    }
                }
            }

            // Arrival generators for each source and class
            arrivalGens = new RandomVariateGen[numSources][numClasses];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    ProcessType procType = arrivalProcessType[srcIdx][k];
                    MatrixCell proc = arrivalProc[srcIdx][k];

                    RandomVariateGen gen = null;
                    MRG32k3a stream = new MRG32k3a();
                    if (seed > 0) {
                        long offset = ((long) srcIdx * numClasses + k) * 10;
                        stream.setSeed(new long[] { seed + offset, seed + offset + 1,
                                seed + offset + 2, seed + offset + 3,
                                seed + offset + 4, seed + offset + 5 });
                    }

                    if (procType == ProcessType.DET) {
                        if (lambdas[srcIdx][k] > 0) {
                            double mean = 1.0 / lambdas[srcIdx][k];
                            gen = new ConstantGen(stream, mean);
                        } else {
                            throw new RuntimeException("LDES: Deterministic arrival distribution for source "
                                    + srcIdx + ", class " + k + " has invalid rate "
                                    + lambdas[srcIdx][k]);
                        }
                    } else if (procType == ProcessType.EXP && lambdas[srcIdx][k] > 0) {
                        gen = new umontreal.ssj.randvar.ExponentialGen(stream, lambdas[srcIdx][k]);
                    } else {
                        gen = createNonMarkovianArrivalGen(procType, proc,
                                lambdas[srcIdx][k],
                                sn.scv.get(sourceStations.get(srcIdx), k), stream);
                    }
                    if (gen == null && procType == ProcessType.IMMEDIATE) {
                        gen = new ConstantGen(stream, 1.0 / GlobalConstants.Immediate);
                    }
                    arrivalGens[srcIdx][k] = gen;
                }
            }

            // Determine process types for service
            serviceProcessType = new ProcessType[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int istStation = serviceStations.get(svcIdx);
                jline.lang.nodes.Station station = sn.stations.get(istStation);
                for (int k = 0; k < numClasses; k++) {
                    jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                    Map<jline.lang.JobClass, ProcessType> rmap =
                            (sn.procid != null) ? sn.procid.get(station) : null;
                    ProcessType pt = (rmap != null) ? rmap.get(jobClass) : null;
                    serviceProcessType[svcIdx][k] = (pt != null) ? pt : ProcessType.DISABLED;
                }
            }

            // Initialize PH process matrices from sn.proc
            serviceProc = new MatrixCell[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int istStation = serviceStations.get(svcIdx);
                jline.lang.nodes.Station station = sn.stations.get(istStation);
                for (int k = 0; k < numClasses; k++) {
                    jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                    Map<jline.lang.JobClass, MatrixCell> rmap =
                            (sn.proc != null) ? sn.proc.get(station) : null;
                    serviceProc[svcIdx][k] = (rmap != null) ? rmap.get(jobClass) : null;
                }
            }

            // Initialize service rate schedules (NHPP) from serviceProc
            serviceSchedule = new double[numServiceNodes][numClasses][][];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    ProcessType pt = serviceProcessType[svcIdx][k];
                    if (pt == ProcessType.NHPP) {
                        serviceSchedule[svcIdx][k] =
                                unpackRateSchedule(serviceProc[svcIdx][k]);
                    }
                }
            }

            // Initialize random generators for PH/MAP/MMAP service sampling
            serviceRng = new Random[numServiceNodes][numClasses];
            serviceMapSampler = new Map_sample.MapSampler[numServiceNodes][numClasses];
            serviceMeSampler = new Me_sample.MeSampler[numServiceNodes][numClasses];
            serviceRapSampler = new Rap_sample.RapSampler[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    ProcessType procType = serviceProcessType[svcIdx][k];
                    if (procType == ProcessType.PH || procType == ProcessType.APH
                            || procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN
                            || procType == ProcessType.COX2 || procType == ProcessType.MAP
                            || procType == ProcessType.MMPP2 || procType == ProcessType.MMAP
                            || procType == ProcessType.ME || procType == ProcessType.RAP
                            || procType == ProcessType.NHPP) {
                        if (seed > 0) {
                            long offset = ((long) (numSources + svcIdx) * numClasses + k) * 10 + 1000;
                            serviceRng[svcIdx][k] = new Random(seed + offset);
                        } else {
                            serviceRng[svcIdx][k] = new Random();
                        }
                    } else {
                        serviceRng[svcIdx][k] = null;
                    }
                }
            }

            // Service generators
            serviceGens = new RandomVariateGen[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    ProcessType procType = serviceProcessType[svcIdx][k];
                    MatrixCell proc = serviceProc[svcIdx][k];

                    RandomVariateGen gen = null;
                    MRG32k3a stream = new MRG32k3a();
                    if (seed > 0) {
                        long offset = ((long) (numSources + svcIdx) * numClasses + k) * 10 + 1000;
                        stream.setSeed(new long[] { seed + offset, seed + offset + 1,
                                seed + offset + 2, seed + offset + 3,
                                seed + offset + 4, seed + offset + 5 });
                    }

                    if (procType == ProcessType.DET) {
                        if (mus[svcIdx][k] > 0 && mus[svcIdx][k] < Double.MAX_VALUE) {
                            double mean = 1.0 / mus[svcIdx][k];
                            gen = new ConstantGen(stream, mean);
                        } else {
                            throw new RuntimeException("LDES: Deterministic service distribution for station "
                                    + svcIdx + ", class " + k + " has invalid rate "
                                    + mus[svcIdx][k]);
                        }
                    } else if (procType == ProcessType.EXP
                            && mus[svcIdx][k] > 0 && mus[svcIdx][k] < Double.MAX_VALUE) {
                        gen = new umontreal.ssj.randvar.ExponentialGen(stream, mus[svcIdx][k]);
                    } else {
                        gen = createNonMarkovianArrivalGen(procType, proc,
                                mus[svcIdx][k],
                                sn.scv.get(serviceStations.get(svcIdx), k), stream);
                    }
                    if (gen == null && procType == ProcessType.IMMEDIATE) {
                        gen = new ConstantGen(stream, 1.0 / GlobalConstants.Immediate);
                    }
                    serviceGens[svcIdx][k] = gen;
                }
            }

            // Initialize trace samplers for Replayer distributions (arrivals)
            arrivalTraceSamplers = new TraceSampler[numSources][numClasses];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                int istStation = sourceStations.get(srcIdx);
                jline.lang.nodes.Station station = sn.stations.get(istStation);
                for (int k = 0; k < numClasses; k++) {
                    ProcessType procType = arrivalProcessType[srcIdx][k];
                    if (procType == ProcessType.REPLAYER && station instanceof Source) {
                        jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                        jline.lang.processes.Distribution distr =
                                ((Source) station).getArrivalDistribution(jobClass);
                        if (distr instanceof Replayer) {
                            double[] data = ((Replayer) distr).getData();
                            arrivalTraceSamplers[srcIdx][k] = new TraceSampler(data);
                        }
                    }
                }
            }

            // Initialize trace samplers for Replayer distributions (service)
            serviceTraceSamplers = new TraceSampler[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int istStation = serviceStations.get(svcIdx);
                jline.lang.nodes.Station station = sn.stations.get(istStation);
                for (int k = 0; k < numClasses; k++) {
                    ProcessType procType = serviceProcessType[svcIdx][k];
                    if (procType == ProcessType.REPLAYER) {
                        jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                        jline.lang.processes.Distribution distr =
                                (station.getServer() != null)
                                        ? station.getServer().getServiceDistribution(jobClass)
                                        : null;
                        if (distr instanceof Replayer) {
                            double[] data = ((Replayer) distr).getData();
                            serviceTraceSamplers[svcIdx][k] = new TraceSampler(data);
                        }
                    }
                }
            }

            // Initialize setup and delayoff generators
            hasSetupDelayoff = new boolean[numServiceNodes];

            setupGens = new RandomVariateGen[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int istStation = serviceStations.get(svcIdx);
                jline.lang.nodes.Station station = sn.stations.get(istStation);
                for (int k = 0; k < numClasses; k++) {
                    if (station instanceof Queue && ((Queue) station).isDelayOffEnabled()) {
                        hasSetupDelayoff[svcIdx] = true;
                        jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                        jline.lang.processes.Distribution setupDist =
                                ((Queue) station).getSetupTime(jobClass);
                        if (setupDist != null && setupDist.getMean() > 0) {
                            MRG32k3a stream = new MRG32k3a();
                            if (seed > 0) {
                                long offset = ((long) (numSources + numServiceNodes + svcIdx)
                                        * numClasses + k) * 10 + 3000;
                                stream.setSeed(new long[] { seed + offset, seed + offset + 1,
                                        seed + offset + 2, seed + offset + 3,
                                        seed + offset + 4, seed + offset + 5 });
                            }
                            setupGens[svcIdx][k] = createSetupDelayoffGen(setupDist, stream);
                        }
                    }
                }
            }

            delayoffGens = new RandomVariateGen[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int istStation = serviceStations.get(svcIdx);
                jline.lang.nodes.Station station = sn.stations.get(istStation);
                for (int k = 0; k < numClasses; k++) {
                    if (station instanceof Queue && ((Queue) station).isDelayOffEnabled()) {
                        jline.lang.JobClass jobClass = sn.jobclasses.get(k);
                        jline.lang.processes.Distribution delayoffDist =
                                ((Queue) station).getDelayOffTime(jobClass);
                        if (delayoffDist != null && delayoffDist.getMean() > 0) {
                            MRG32k3a stream = new MRG32k3a();
                            if (seed > 0) {
                                long offset = ((long) (numSources + numServiceNodes + svcIdx)
                                        * numClasses + k) * 10 + 4000;
                                stream.setSeed(new long[] { seed + offset, seed + offset + 1,
                                        seed + offset + 2, seed + offset + 3,
                                        seed + offset + 4, seed + offset + 5 });
                            }
                            delayoffGens[svcIdx][k] = createSetupDelayoffGen(delayoffDist, stream);
                        }
                    }
                }
            }

            // Initialize antithetic variates generators if enabled
            initializeAntitheticGenerators();

            // Initialize control variates tracking if enabled
            initializeControlVariates();
        }

        /**
         * Build a {@link RandomVariateGen} for non-Markovian arrival/service distributions.
         * Handles UNIFORM/GAMMA/WEIBULL/LOGNORMAL/PARETO/ERLANG/POISSON/BERNOULLI/BINOMIAL/
         * GEOMETRIC/IMMEDIATE branches shared between arrival and service initialization.
         *
         * <p>UNIFORM, GAMMA, WEIBULL, LOGNORMAL and PARETO are reconstructed from
         * {@code rate} and {@code scv}, NOT from {@code proc}. Those five are exactly
         * the types {@code Network.isRawParameterProcess} lists, and for them
         * {@code Network.refreshService} OVERWRITES the raw {@code getProcess()} pair
         * with an Erlang moment fit (map_erlang), mirroring MATLAB. So by the time sn
         * reaches this engine, {@code proc.get(0)} is a D0 generator whose (0,0) entry
         * is NEGATIVE, not a shape or a lower bound. Reading it as a parameter is what
         * produced "alpha &lt;= 0", "sigma &lt;= 0", and Uniform service times sampled
         * from U(-n, n). Verified on M/Uniform/1 with Uniform(0.5, 1.5): sn.proc is
         * 12-phase Erlang with D0(0,0) = -12, while sn.rates = 1.0 and sn.scv = 1/12
         * are exact, in the JAR and in MATLAB alike.
         *
         * <p>Mean and SCV determine all five uniquely (Weibull via the standard
         * c^-1.086 shape fit), so this recovers the requested distribution rather than
         * simulating its Erlang approximation. The formulas match the two sibling
         * constructions in this file, {@code firingGenFromMeanScv} and
         * {@code createImpatienceGenerator}.
         */
        private RandomVariateGen createNonMarkovianArrivalGen(ProcessType procType,
                                                              MatrixCell proc,
                                                              double momRate,
                                                              double momScv,
                                                              MRG32k3a stream) {
            RandomVariateGen gen = null;
            if (procType == ProcessType.UNIFORM || procType == ProcessType.GAMMA
                    || procType == ProcessType.WEIBULL || procType == ProcessType.LOGNORMAL
                    || procType == ProcessType.PARETO) {
                // Inconsistent moments are a model defect, not something to paper over:
                // returning null here would leave the caller with no generator and a
                // zero service time, and handing the value to SSJ throws from inside
                // the library with no indication of which station or class is at fault.
                if (!(momRate > 0) || !(momRate < Double.MAX_VALUE) || !(momScv > 0)) {
                    throw new RuntimeException("LDES: " + procType
                            + " distribution has unusable moments (momRate=" + momRate
                            + ", momScv=" + momScv + "); expected momRate in (0,inf) and momScv > 0.");
                }
                double mean = 1.0 / momRate;
                if (procType == ProcessType.UNIFORM) {
                    // U[a,b]: mean = (a+b)/2, var = (b-a)^2/12, so momScv <= 1/3 for any
                    // non-negative support, with equality iff a = 0.
                    double halfWidth = mean * FastMath.sqrt(3.0 * momScv);
                    double minVal = mean - halfWidth;
                    double maxVal = mean + halfWidth;
                    if (!(minVal >= 0) || !(maxVal > minVal)) {
                        throw new RuntimeException("LDES: UNIFORM distribution with mean "
                                + mean + " and momScv " + momScv + " implies support ["
                                + minVal + ", " + maxVal + "]; a non-negative uniform "
                                + "requires momScv <= 1/3.");
                    }
                    gen = new umontreal.ssj.randvar.UniformGen(stream, minVal, maxVal);
                } else if (procType == ProcessType.GAMMA) {
                    double shape = 1.0 / momScv;
                    double lambda = momRate / momScv;
                    gen = new umontreal.ssj.randvar.GammaGen(stream, shape, lambda);
                } else if (procType == ProcessType.WEIBULL) {
                    double c = FastMath.sqrt(momScv);
                    double r = FastMath.pow(c, -1.086);
                    double alpha = mean
                            / org.apache.commons.math3.special.Gamma.gamma(1.0 + 1.0 / r);
                    gen = new umontreal.ssj.randvar.WeibullGen(stream, r, 1.0 / alpha, 0.0);
                } else if (procType == ProcessType.LOGNORMAL) {
                    double c2plus1 = momScv + 1.0;
                    double mu = FastMath.log(mean / FastMath.sqrt(c2plus1));
                    double sigma = FastMath.sqrt(FastMath.log(c2plus1));
                    gen = new umontreal.ssj.randvar.LognormalGen(stream, mu, sigma);
                } else {
                    // Pareto: momScv = 1/(a(a-2)) for a > 2, so a = 1 + sqrt(1 + 1/momScv),
                    // and the scale follows from mean = a*m/(a-1).
                    double shape = FastMath.sqrt(1.0 + 1.0 / momScv) + 1.0;
                    double scale = mean * (shape - 1.0) / shape;
                    gen = new umontreal.ssj.randvar.ParetoGen(stream, shape, scale);
                }
            } else if (procType == ProcessType.ERLANG) {
                if (proc != null && proc.size() >= 2
                        && (proc.get(0).getNumRows() > 1 || proc.get(0).get(0, 0) < 0)) {
                    Matrix d0 = proc.get(0);
                    int phases = d0.getNumRows();
                    double lambda = -d0.get(0, 0);
                    if (phases > 0 && lambda > 0.0) {
                        gen = new umontreal.ssj.randvar.ErlangGen(stream, phases, lambda);
                    }
                } else if (proc != null && proc.size() >= 2) {
                    double mean = proc.get(0).get(0, 0);
                    double scv = proc.get(1).get(0, 0);
                    if (mean > 0 && scv > 0) {
                        int shape = Math.max((int) (1.0 / scv), 1);
                        double scale = mean / shape;
                        gen = new umontreal.ssj.randvar.ErlangGen(stream, shape, 1.0 / scale);
                    }
                }
            } else if (procType == ProcessType.POISSON) {
                if (proc != null && proc.size() >= 1) {
                    double mean = proc.get(0).get(0, 0);
                    gen = new umontreal.ssj.randvar.PoissonGen(stream, mean);
                }
            } else if (procType == ProcessType.BERNOULLI) {
                if (proc != null && proc.size() >= 1) {
                    double mean = proc.get(0).get(0, 0);
                    gen = new umontreal.ssj.randvar.BernoulliGen(stream, mean);
                }
            } else if (procType == ProcessType.BINOMIAL) {
                if (proc != null && proc.size() >= 2) {
                    double mean = proc.get(0).get(0, 0);
                    double scv = proc.get(1).get(0, 0);
                    double p = 1.0 - scv * mean;
                    if (p >= 0.0 && p <= 1.0) {
                        int n = (int) (mean / p);
                        gen = new umontreal.ssj.randvar.BinomialGen(stream, n, p);
                    }
                }
            } else if (procType == ProcessType.GEOMETRIC) {
                if (proc != null && proc.size() >= 1) {
                    // Geometric is stored as {mean, SCV} and is supported on
                    // {1,2,...} with mean 1/p, so p = 1/mean. Mapping it onto
                    // SSJ's GeometricGen instead would shift the support to
                    // {0,1,...}: the mean would still match but the SCV would
                    // come out as 1+p rather than 1-p.
                    double mean = proc.get(0).get(0, 0);
                    gen = new ShiftedGeometricGen(stream, 1.0 / mean);
                }
            } else if (procType == ProcessType.DMAP) {
                // Discrete-time MAP. The feature set declares DMAP, so a station
                // carrying one must reach the sampler rather than fall through to
                // a null generator; DmapSampleGen is the same sampler the Petri
                // net transition path uses.
                if (proc != null && proc.size() >= 2) {
                    gen = new DmapSampleGen(stream, proc,
                            new Random((long) stream.nextInt(0, Integer.MAX_VALUE)));
                }
            } else if (procType == ProcessType.IMMEDIATE) {
                double immTime = 1.0 / GlobalConstants.Immediate;
                gen = new ConstantGen(stream, immTime);
            }
            return gen;
        }

        /**
         * Initialize antithetic variates generators (Kotlin lines 5025-5111).
         */
        @SuppressWarnings("unchecked")
        private void initializeAntitheticGenerators() {
            if (!useAntitheticVariates) {
                antitheticArrivalGens = new RandomVariateGen[numSources][numClasses];
                antitheticServiceGens = new RandomVariateGen[numServiceNodes][numClasses];
                antitheticArrivalRng = new Random[numSources][numClasses];
                antitheticServiceRng = new Random[numServiceNodes][numClasses];
                return;
            }

            long antitheticSeedOffset = 50000L;

            antitheticArrivalGens = new RandomVariateGen[numSources][numClasses];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    ProcessType procType = arrivalProcessType[srcIdx][k];
                    RandomVariateGen gen = null;
                    MRG32k3a stream = new MRG32k3a();
                    if (seed > 0) {
                        long offset = ((long) srcIdx * numClasses + k) * 10 + antitheticSeedOffset;
                        stream.setSeed(new long[] { seed + offset, seed + offset + 1,
                                seed + offset + 2, seed + offset + 3,
                                seed + offset + 4, seed + offset + 5 });
                    }
                    umontreal.ssj.rng.AntitheticStream antiStream =
                            new umontreal.ssj.rng.AntitheticStream(stream);
                    if (procType == ProcessType.EXP && lambdas[srcIdx][k] > 0) {
                        gen = new umontreal.ssj.randvar.ExponentialGen(antiStream,
                                lambdas[srcIdx][k]);
                    }
                    antitheticArrivalGens[srcIdx][k] = gen;
                }
            }

            antitheticServiceGens = new RandomVariateGen[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    ProcessType procType = serviceProcessType[svcIdx][k];
                    RandomVariateGen gen = null;
                    MRG32k3a stream = new MRG32k3a();
                    if (seed > 0) {
                        long offset = ((long) numSources * numClasses
                                + (long) svcIdx * numClasses + k) * 10 + antitheticSeedOffset;
                        stream.setSeed(new long[] { seed + offset, seed + offset + 1,
                                seed + offset + 2, seed + offset + 3,
                                seed + offset + 4, seed + offset + 5 });
                    }
                    umontreal.ssj.rng.AntitheticStream antiStream =
                            new umontreal.ssj.rng.AntitheticStream(stream);
                    if (procType == ProcessType.EXP
                            && mus[svcIdx][k] > 0 && mus[svcIdx][k] < Double.MAX_VALUE) {
                        gen = new umontreal.ssj.randvar.ExponentialGen(antiStream,
                                mus[svcIdx][k]);
                    }
                    antitheticServiceGens[svcIdx][k] = gen;
                }
            }

            antitheticArrivalRng = new Random[numSources][numClasses];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    ProcessType procType = arrivalProcessType[srcIdx][k];
                    if (procType == ProcessType.PH || procType == ProcessType.APH
                            || procType == ProcessType.MAP || procType == ProcessType.MMPP2
                            || procType == ProcessType.MMAP) {
                        Random rng = new Random();
                        if (seed > 0) {
                            rng.setSeed(seed + antitheticSeedOffset
                                    + (long) (srcIdx * numClasses + k) * 100);
                        }
                        antitheticArrivalRng[srcIdx][k] = rng;
                    }
                }
            }

            antitheticServiceRng = new Random[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    ProcessType procType = serviceProcessType[svcIdx][k];
                    if (procType == ProcessType.PH || procType == ProcessType.APH
                            || procType == ProcessType.MAP || procType == ProcessType.MMPP2
                            || procType == ProcessType.MMAP) {
                        Random rng = new Random();
                        if (seed > 0) {
                            rng.setSeed(seed + antitheticSeedOffset
                                    + (long) (numSources * numClasses
                                            + svcIdx * numClasses + k) * 100);
                        }
                        antitheticServiceRng[svcIdx][k] = rng;
                    }
                }
            }
        }

        /**
         * Initialize control variates tracking structures (Kotlin lines 5116-5139).
         */
        private void initializeControlVariates() {
            arrivalSampleSum = new double[numSources][numClasses];
            arrivalSampleCount = new long[numSources][numClasses];
            serviceSampleSum = new double[numServiceNodes][numClasses];
            serviceSampleCount = new long[numServiceNodes][numClasses];

            arrivalExpectedMean = new double[numSources][numClasses];
            for (int srcIdx = 0; srcIdx < numSources; srcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    arrivalExpectedMean[srcIdx][k] = (lambdas[srcIdx][k] > 0)
                            ? 1.0 / lambdas[srcIdx][k] : 0.0;
                }
            }

            serviceExpectedMean = new double[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    serviceExpectedMean[svcIdx][k] =
                            (mus[svcIdx][k] > 0 && mus[svcIdx][k] < Double.MAX_VALUE)
                                    ? 1.0 / mus[svcIdx][k] : 0.0;
                }
            }

            initializeHeterogeneousServers();
        }

        /**
         * Initialize heterogeneous server support (Kotlin lines 5146-5262).
         */
        @SuppressWarnings("unchecked")
        private void initializeHeterogeneousServers() {
            numServerTypes = new int[numServiceNodes];
            serversPerType = new int[numServiceNodes][];
            serverCompat = new boolean[numServiceNodes][][];
            busyCountPerType = new int[numServiceNodes][];
            serverToType = new int[numServiceNodes][];
            heteroSchedPolicy = new HeteroSchedPolicy[numServiceNodes];
            heteroServiceGens = new RandomVariateGen[numServiceNodes][][];
            heteroMus = new double[numServiceNodes][][];
            heteroServiceProcType = new ProcessType[numServiceNodes][][];
            heteroServiceProc = new MatrixCell[numServiceNodes][][];
            heteroServiceRng = new Random[numServiceNodes][][];
            heteroServiceMapSampler = new Map_sample.MapSampler[numServiceNodes][][];
            heteroServiceMeSampler = new Me_sample.MeSampler[numServiceNodes][][];
            heteroServiceRapSampler = new Rap_sample.RapSampler[numServiceNodes][][];
            serverTypeOrder = (List<Integer>[]) new List<?>[numServiceNodes];
            alfsOrder = new int[numServiceNodes][];

            // Defaults for homogeneous case
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                serversPerType[svcIdx] = new int[0];
                serverCompat[svcIdx] = new boolean[0][];
                busyCountPerType[svcIdx] = new int[0];
                serverToType[svcIdx] = new int[0];
                heteroServiceGens[svcIdx] = new RandomVariateGen[0][0];
                heteroMus[svcIdx] = new double[0][];
                heteroServiceProcType[svcIdx] = new ProcessType[0][];
                heteroServiceProc[svcIdx] = new MatrixCell[0][];
                heteroServiceRng[svcIdx] = new Random[0][];
                heteroServiceMapSampler[svcIdx] = new Map_sample.MapSampler[0][];
                heteroServiceMeSampler[svcIdx] = new Me_sample.MeSampler[0][];
                heteroServiceRapSampler[svcIdx] = new Rap_sample.RapSampler[0][];
                serverTypeOrder[svcIdx] = new ArrayList<Integer>();
                alfsOrder[svcIdx] = new int[0];
            }

            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                if (isDelayNode.get(svcIdx)) continue;

                int stationIdx = serviceStations.get(svcIdx);
                if (stationIdx >= sn.stations.size()) continue;

                jline.lang.nodes.Station station = sn.stations.get(stationIdx);
                jline.lang.nodeparam.ServiceNodeParam snp = sn.getServiceParam(station);
                if (snp == null) continue;

                int nTypes = snp.nservertypes;
                if (nTypes <= 0) continue;

                numServerTypes[svcIdx] = nTypes;

                // Get servers per type
                Matrix serversMatrix = snp.serverspertype;
                serversPerType[svcIdx] = new int[nTypes];
                for (int typeId = 0; typeId < nTypes; typeId++) {
                    if (serversMatrix != null) {
                        serversPerType[svcIdx][typeId] = (int) serversMatrix.get(typeId);
                    } else {
                        serversPerType[svcIdx][typeId] = 1;
                    }
                }

                // Get compatibility matrix
                Matrix compatMatrix = snp.servercompat;
                serverCompat[svcIdx] = new boolean[nTypes][numClasses];
                for (int typeId = 0; typeId < nTypes; typeId++) {
                    for (int classId = 0; classId < numClasses; classId++) {
                        double compatible = (compatMatrix != null)
                                ? compatMatrix.get(typeId, classId) : 0.0;
                        serverCompat[svcIdx][typeId][classId] = compatible > 0.5;
                    }
                }

                busyCountPerType[svcIdx] = new int[nTypes];

                // Build server ID to type mapping
                int totalServers = 0;
                for (int t = 0; t < nTypes; t++) totalServers += serversPerType[svcIdx][t];
                serverToType[svcIdx] = new int[totalServers];
                int globalServerId = 0;
                for (int typeId = 0; typeId < nTypes; typeId++) {
                    int count = serversPerType[svcIdx][typeId];
                    for (int j = 0; j < count; j++) {
                        if (globalServerId < totalServers) {
                            serverToType[svcIdx][globalServerId] = typeId;
                            globalServerId++;
                        }
                    }
                }

                // Heterogeneous scheduling policy
                HeteroSchedPolicy policy = snp.heteroschedpolicy;
                heteroSchedPolicy[svcIdx] = (policy != null) ? policy : HeteroSchedPolicy.ORDER;

                // Server type order for ALIS/FAIRNESS round-robin
                List<Integer> typeOrder = new ArrayList<Integer>();
                for (int t = 0; t < nTypes; t++) typeOrder.add(t);
                serverTypeOrder[svcIdx] = typeOrder;

                // ALFS order: sorted by number of compatible classes (ascending)
                final int svcIdxFinal = svcIdx;
                Integer[] indices = new Integer[nTypes];
                for (int t = 0; t < nTypes; t++) indices[t] = t;
                Arrays.sort(indices, new Comparator<Integer>() {
                    @Override
                    public int compare(Integer a, Integer b) {
                        int ca = 0, cb = 0;
                        for (boolean v : serverCompat[svcIdxFinal][a]) if (v) ca++;
                        for (boolean v : serverCompat[svcIdxFinal][b]) if (v) cb++;
                        return Integer.compare(ca, cb);
                    }
                });
                alfsOrder[svcIdx] = new int[nTypes];
                for (int t = 0; t < nTypes; t++) alfsOrder[svcIdx][t] = indices[t];

                // Heterogeneous service rates
                Map<Integer, Map<Integer, Double>> stationRates = snp.heterorates;
                heteroMus[svcIdx] = new double[nTypes][numClasses];
                for (int typeId = 0; typeId < nTypes; typeId++) {
                    for (int classId = 0; classId < numClasses; classId++) {
                        Double rate = null;
                        if (stationRates != null) {
                            Map<Integer, Double> typeMap = stationRates.get(typeId);
                            if (typeMap != null) rate = typeMap.get(classId);
                        }
                        heteroMus[svcIdx][typeId][classId] =
                                (rate != null) ? rate : Double.MAX_VALUE;
                    }
                }

                // Heterogeneous process types
                Map<Integer, Map<Integer, ProcessType>> stationProcId =
                        snp.heteroprocid;
                heteroServiceProcType[svcIdx] = new ProcessType[nTypes][numClasses];
                for (int typeId = 0; typeId < nTypes; typeId++) {
                    for (int classId = 0; classId < numClasses; classId++) {
                        ProcessType pt = null;
                        if (stationProcId != null) {
                            Map<Integer, ProcessType> typeMap = stationProcId.get(typeId);
                            if (typeMap != null) pt = typeMap.get(classId);
                        }
                        heteroServiceProcType[svcIdx][typeId][classId] =
                                (pt != null) ? pt : ProcessType.EXP;
                    }
                }

                // Heterogeneous PH processes
                Map<Integer, Map<Integer, MatrixCell>> stationProc =
                        snp.heteroproc;
                heteroServiceProc[svcIdx] = new MatrixCell[nTypes][numClasses];
                for (int typeId = 0; typeId < nTypes; typeId++) {
                    for (int classId = 0; classId < numClasses; classId++) {
                        MatrixCell mc = null;
                        if (stationProc != null) {
                            Map<Integer, MatrixCell> typeMap = stationProc.get(typeId);
                            if (typeMap != null) mc = typeMap.get(classId);
                        }
                        heteroServiceProc[svcIdx][typeId][classId] = mc;
                    }
                }

                // RNG for PH sampling
                heteroServiceRng[svcIdx] = new Random[nTypes][numClasses];
                heteroServiceMapSampler[svcIdx] = new Map_sample.MapSampler[nTypes][numClasses];
                heteroServiceMeSampler[svcIdx] = new Me_sample.MeSampler[nTypes][numClasses];
                heteroServiceRapSampler[svcIdx] = new Rap_sample.RapSampler[nTypes][numClasses];
                for (int typeId = 0; typeId < nTypes; typeId++) {
                    for (int classId = 0; classId < numClasses; classId++) {
                        if (seed > 0) {
                            long offset = ((long) svcIdx * 1000L
                                    + (long) typeId * 100L + classId) + 8000;
                            heteroServiceRng[svcIdx][typeId][classId] = new Random(seed + offset);
                        } else {
                            heteroServiceRng[svcIdx][typeId][classId] = new Random();
                        }
                    }
                }

                // Create heterogeneous service generators
                heteroServiceGens[svcIdx] = new RandomVariateGen[nTypes][numClasses];
                for (int typeId = 0; typeId < nTypes; typeId++) {
                    for (int classId = 0; classId < numClasses; classId++) {
                        heteroServiceGens[svcIdx][typeId][classId] =
                                createHeteroServiceGenerator(svcIdx, typeId, classId);
                    }
                }
            }
        }

        /**
         * Service-time generator for a specific (svcIdx, typeId, classId) triple
         * (Kotlin lines 5267-5308).
         */
        private RandomVariateGen createHeteroServiceGenerator(int svcIdx, int typeId, int classId) {
            double rate = heteroMus[svcIdx][typeId][classId];
            if (rate <= 0 || rate >= Double.MAX_VALUE) return null;

            ProcessType procType = heteroServiceProcType[svcIdx][typeId][classId];
            MatrixCell proc = heteroServiceProc[svcIdx][typeId][classId];

            MRG32k3a stream = new MRG32k3a();
            if (seed > 0) {
                long offset = ((long) svcIdx * 1000L + (long) typeId * 100L + classId) + 7000;
                stream.setSeed(new long[] { seed + offset, seed + offset + 1,
                        seed + offset + 2, seed + offset + 3,
                        seed + offset + 4, seed + offset + 5 });
            }

            if (procType == ProcessType.DET) {
                return new ConstantGen(stream, 1.0 / rate);
            } else if (procType == ProcessType.EXP) {
                return new umontreal.ssj.randvar.ExponentialGen(stream, rate);
            } else if (procType == ProcessType.ERLANG) {
                if (proc != null && proc.size() >= 2) {
                    Matrix d0 = proc.get(0);
                    int nPhases = d0.getNumCols();
                    double phaseRate = rate * nPhases;
                    return new umontreal.ssj.randvar.ErlangGen(stream, nPhases, phaseRate);
                }
                return null;
            } else if (procType == ProcessType.HYPEREXP
                    || procType == ProcessType.PH || procType == ProcessType.APH
                    || procType == ProcessType.COXIAN || procType == ProcessType.COX2
                    || procType == ProcessType.MAP) {
                // Use map_sample at runtime
                return null;
            } else {
                if (rate > 0 && rate < Double.MAX_VALUE) {
                    return new umontreal.ssj.randvar.ExponentialGen(stream, rate);
                }
                return null;
            }
        }

        /**
         * Initialize per-station state variables (Kotlin lines 5310-5463).
         */
        @SuppressWarnings("unchecked")
        private void initializeState() {
            // Wait queues with comparator based on scheduling strategy
            waitQueues = (PriorityQueue<Customer>[]) new PriorityQueue<?>[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                waitQueues[svcIdx] = new PriorityQueue<Customer>(
                        getComparatorForStrategy(schedStrategies[svcIdx], svcIdx));
            }

            // Server busy state
            serverBusy = new boolean[numServiceNodes][];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                serverBusy[svcIdx] = isDelayNode.get(svcIdx)
                        ? new boolean[0] : new boolean[numServers[svcIdx]];
            }
            // Server blocked state for synchronous calls
            serverBlocked = new boolean[numServiceNodes][];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                serverBlocked[svcIdx] = isDelayNode.get(svcIdx)
                        ? new boolean[0] : new boolean[numServers[svcIdx]];
            }
            customersInService = new int[numServiceNodes];

            responseTimeTally = new Tally[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    responseTimeTally[svcIdx][k] =
                            new Tally("Response time S" + svcIdx + " C" + k);
                }
            }
            responseTimeSamples = (List<Double>[][]) new List<?>[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    responseTimeSamples[svcIdx][k] = new ArrayList<Double>();
                }
            }
            completedCustomers = new int[numServiceNodes][numClasses];
            totalQueueTime = new double[numServiceNodes][numClasses];
            lastQueueUpdateTime = new double[numServiceNodes][numClasses];
            currentQueueLength = new int[numServiceNodes][numClasses];

            // Markov reward setup
            this.hasReward = (sn.reward != null && !sn.reward.isEmpty());
            if (hasReward) {
                this.rewardNames = new ArrayList<String>(sn.reward.keySet());
                this.rewardFnArr = new jline.lang.reward.RewardFunction[rewardNames.size()];
                for (int ri = 0; ri < rewardNames.size(); ri++) {
                    this.rewardFnArr[ri] = sn.reward.get(rewardNames.get(ri));
                }
                this.rewardArea = new double[rewardNames.size()];
                this.rewardStateCache = new HashMap<String, double[]>();
                this.rewardTranSeries = new ArrayList<double[]>();
            }
            this.exportStateHistogram = (this.options instanceof LDESOptions)
                    && ((LDESOptions) this.options).exportStateHistogram;
            if (hasReward || exportStateHistogram) {
                this.rewardRowCols = sn.nstations * numClasses;
                this.rewardLastUpdateTime = ssjSim.time();
                this.rewardTotalTime = 0.0;
                if (exportStateHistogram) {
                    this.stateHistogram = new java.util.LinkedHashMap<String, double[]>();
                    this.stateTranSeries = new ArrayList<double[]>();
                }
            }

            totalBusyTime = new double[numServiceNodes][numClasses];
            lastBusyUpdateTime = new double[numServiceNodes][numClasses];
            currentBusyServers = new int[numServiceNodes][numClasses];

            totalBlockingTime = new double[numServiceNodes][numClasses];
            currentBlockedServers = new int[numServiceNodes][numClasses];

            basBlockedAtDest = new int[numServiceNodes][numClasses];
            bbsBlockedAtDest = new int[numServiceNodes][numClasses];
            fcrBlockedAtDest = new int[numServiceNodes][numClasses];

            // Setup/delayoff state — all servers start ACTIVE
            serverState = new ServerState[numServiceNodes][];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                if (isDelayNode.get(svcIdx)) {
                    serverState[svcIdx] = new ServerState[0];
                } else {
                    serverState[svcIdx] = new ServerState[numServers[svcIdx]];
                    for (int s = 0; s < numServers[svcIdx]; s++) {
                        serverState[svcIdx][s] = ServerState.ACTIVE;
                    }
                }
            }
            serverLastClass = new int[numServiceNodes][];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                if (isDelayNode.get(svcIdx)) {
                    serverLastClass[svcIdx] = new int[0];
                } else {
                    serverLastClass[svcIdx] = new int[numServers[svcIdx]];
                    Arrays.fill(serverLastClass[svcIdx], -1);
                }
            }
            pendingDelayoffEvents = new Event[numServiceNodes][];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                pendingDelayoffEvents[svcIdx] = isDelayNode.get(svcIdx)
                        ? new Event[0] : new Event[numServers[svcIdx]];
            }

            totalSetupTime = new double[numServiceNodes][numClasses];
            lastSetupUpdateTime = new double[numServiceNodes][numClasses];
            currentServersInSetup = new int[numServiceNodes][numClasses];

            totalDelayoffTime = new double[numServiceNodes][numClasses];
            lastDelayoffUpdateTime = new double[numServiceNodes][numClasses];
            currentServersInDelayoff = new int[numServiceNodes][numClasses];

            systemResponseTimeTally = new Tally[numClasses];
            systemTardinessTally = new Tally[numClasses];
            for (int k = 0; k < numClasses; k++) {
                systemResponseTimeTally[k] = new Tally("System response time C" + k);
                systemTardinessTally[k] = new Tally("System tardiness C" + k);
            }
            tardinessTally = new Tally[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    tardinessTally[svcIdx][k] = new Tally("Tardiness node" + svcIdx + " C" + k);
                }
            }
            systemCompletedCustomers = new int[numClasses];

            droppedCustomers = new int[numServiceNodes][numClasses];

            arrivedCustomers = new int[numServiceNodes][numClasses];

            // Finite capacity region tracking
            currentJobsInRegion = new int[numRegions][numClasses];
            currentMemInRegion = new double[numRegions];
            droppedByRegion = new int[numRegions][numClasses];

            totalRegionJobTime = new double[numRegions][numClasses];
            totalRegionWeightTime = new double[numRegions][numClasses];
            totalRegionMemTime = new double[numRegions][numClasses];
            lastRegionUpdateTime = new double[numRegions];
            regionCompletions = new int[numRegions][numClasses];
            regionResponseTimeTally = new Tally[numRegions][numClasses];
            for (int f = 0; f < numRegions; f++) {
                for (int k = 0; k < numClasses; k++) {
                    regionResponseTimeTally[f][k] =
                            new Tally("FCR" + f + " response time C" + k);
                }
            }

            // FCR blocking (waiting queue) support
            fcRegionBlockedQueue = (LinkedList<BlockedCustomer>[]) new LinkedList<?>[numRegions];
            for (int f = 0; f < numRegions; f++) {
                fcRegionBlockedQueue[f] = new LinkedList<BlockedCustomer>();
            }
            blockedInRegion = new int[numRegions][numClasses];

            // FCR class weights from NetworkStruct
            fcRegionClassWeights = new double[numRegions][numClasses];
            for (int f = 0; f < numRegions; f++) {
                for (int k = 0; k < numClasses; k++) {
                    fcRegionClassWeights[f][k] = (sn.regionweight != null)
                            ? sn.regionweight.get(f, k) : 1.0;
                }
            }

            // FCR per-class memory footprint (classSize) for the memory constraint.
            // Mirrors JMT's globalMemoryConstraint, which charges classSize units
            // of the memory budget per in-region customer of that class.
            fcRegionClassSize = new double[numRegions][numClasses];
            for (int f = 0; f < numRegions; f++) {
                for (int k = 0; k < numClasses; k++) {
                    fcRegionClassSize[f][k] = (sn.regionsz != null)
                            ? sn.regionsz.get(f, k) : 1.0;
                }
            }

            // FCR arrival rate tracking
            lastRegionArrivalTime = new double[numRegions][numClasses];
            regionArrivalCount = new int[numRegions][numClasses];
            regionInterArrivalTimeSum = new double[numRegions][numClasses];

            // Impatience (reneging, balking, retrial) statistics + configuration
            initializeImpatienceSupport();

            // PS scheduling state
            psJobsInService = (List<PSCustomer>[]) new List<?>[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                psJobsInService[svcIdx] = new ArrayList<PSCustomer>();
            }
            psLastUpdateTime = new double[numServiceNodes];
            psLastBusyUpdateTime = new double[numServiceNodes];

            // PAS (pass-and-swap / order-independent) scheduling state
            isPASStation = new boolean[numServiceNodes];
            pasList = (List<Customer>[]) new List<?>[numServiceNodes];
            pasDepartureEvent = new Event[numServiceNodes];
            pasLastBusyUpdateTime = new double[numServiceNodes];
            pasSvcRateFun = (jline.util.SerializableFunction<jline.util.matrix.Matrix, Double>[])
                    new jline.util.SerializableFunction<?, ?>[numServiceNodes];
            pasSwapGraph = new jline.util.matrix.Matrix[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                pasList[svcIdx] = new ArrayList<Customer>();
                if (schedStrategies[svcIdx] == SchedStrategy.PAS) {
                    isPASStation[svcIdx] = true;
                    int nodeIdx = serviceNodes.get(svcIdx);
                    jline.lang.nodes.Node node = sn.nodes.get(nodeIdx);
                    jline.lang.NodeParam nodeParam =
                            (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
                    if (nodeParam instanceof QueueNodeParam) {
                        QueueNodeParam qnp = (QueueNodeParam) nodeParam;
                        pasSvcRateFun[svcIdx] = qnp.svcRateFun;
                        pasSwapGraph[svcIdx] = qnp.swapGraph;
                    }
                    if (pasSvcRateFun[svcIdx] == null) {
                        throw new RuntimeException(
                                "LDES: PAS station " + svcIdx
                                        + " has no service rate function mu(c).");
                    }
                }
            }

            // Preemptive LCFS / SRPT scheduling state
            preemptiveJobsInService = (List<PreemptiveCustomer>[]) new List<?>[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                preemptiveJobsInService[svcIdx] = new ArrayList<PreemptiveCustomer>();
            }
            isPreemptiveScheduling = new boolean[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                SchedStrategy s = schedStrategies[svcIdx];
                isPreemptiveScheduling[svcIdx] = isPreemptiveLCFSScheduling(s)
                        || isPreemptiveFCFSScheduling(s) || isSizeBasedPreemptiveScheduling(s);
            }
            preemptedJobHistory = new HashMap<PreemptionKey, PreemptionRecord>();

            // Cache DPS/GPS weights from sn.schedparam
            schedWeights = new double[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int stationIdx = serviceStations.get(svcIdx);
                for (int k = 0; k < numClasses; k++) {
                    double weight = sn.schedparam.get(stationIdx, k);
                    schedWeights[svcIdx][k] =
                            (weight > 0 && !Double.isNaN(weight)) ? weight : 1.0;
                }
            }

            // Cache LPS limits from sn.schedparam (uniform across classes)
            lpsLimits = new int[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                if (schedStrategies[svcIdx] == SchedStrategy.LPS) {
                    int stationIdx = serviceStations.get(svcIdx);
                    double limit = sn.schedparam.get(stationIdx, 0);
                    lpsLimits[svcIdx] =
                            (limit > 0 && !Double.isNaN(limit) && !Double.isInfinite(limit))
                                    ? (int) limit : Integer.MAX_VALUE;
                } else {
                    lpsLimits[svcIdx] = 0;
                }
            }

            initializeLoadDependentService();

            initializeClassDependence();

            initializePollingState();

            // Needs isLoadDependent/hasCd, so it must follow the two calls above.
            assertScheduleServiceSupported();
        }

        /**
         * Rejects the one station configuration for which a time-varying service
         * rate mu(t) is not simulable by this engine.
         *
         * <p>The engine measures NHPP service in operational time tau = int mu du
         * ({@link #serviceWorkBetween}, {@link #serviceWallAfterWork}). In tau the
         * process is unit-rate, so a completion is exact whether service runs to
         * completion or is interrupted (processor sharing, preemptive resume) or
         * rescaled (load- or class-dependent rate, folded into the effective rate
         * at delivery). All homogeneous disciplines are therefore supported.
         *
         * <p>The remaining gap is per-server-type heterogeneous service: that path
         * ({@link #generateHeteroServiceTime}) has no NHPP branch and would sample
         * a zero duration, so it is rejected. Failing loudly is deliberate: the
         * alternative is a silently wrong sample path.
         */
        private void assertScheduleServiceSupported() {
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                if (numServerTypes[svcIdx] == 0) {
                    continue;
                }
                for (int k = 0; k < numClasses; k++) {
                    if (serviceProcessType[svcIdx][k] != ProcessType.NHPP) {
                        continue;
                    }
                    String stationName =
                            sn.stations.get(serviceStations.get(svcIdx)).getName();
                    String className = sn.jobclasses.get(k).getName();
                    throw new RuntimeException(
                            "LDES: NHPP service is not supported at station '"
                            + stationName + "' for class '" + className + "' with"
                            + " heterogeneous server types. A per-server-type"
                            + " time-varying rate has no sampler; use homogeneous"
                            + " servers, for which mu(t) is fully supported under"
                            + " any scheduling discipline.");
                }
            }
        }

        /**
         * Initialize impatience support (Kotlin lines 5469-5557).
         */
        @SuppressWarnings("unchecked")
        private void initializeImpatienceSupport() {
            // Reneging statistics
            renegedCustomers = new int[numServiceNodes][numClasses];
            totalRenegingWaitTime = new double[numServiceNodes][numClasses];

            // Balking statistics
            balkedCustomers = new int[numServiceNodes][numClasses];

            // Retrial / orbit statistics
            orbitJobs = (List<OrbitJob>[]) new List<?>[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                orbitJobs[svcIdx] = new ArrayList<OrbitJob>();
            }
            retriedCustomers = new int[numServiceNodes][numClasses];
            maxRetriesExceeded = new int[numServiceNodes][numClasses];
            currentOrbitSize = new int[numServiceNodes][numClasses];
            totalOrbitTime = new double[numServiceNodes][numClasses];
            lastOrbitUpdateTime = new double[numServiceNodes][numClasses];

            // Configuration caches
            hasPatienceConfig = new boolean[numServiceNodes][numClasses];
            patienceGens = new RandomVariateGen[numServiceNodes][numClasses];
            hasBalkingConfig = new boolean[numServiceNodes][numClasses];
            hasRetrialConfig = new boolean[numServiceNodes][numClasses];
            retrialGens = new RandomVariateGen[numServiceNodes][numClasses];
            retrialMaxAttemptsConfig = new int[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                Arrays.fill(retrialMaxAttemptsConfig[svcIdx], -1);
            }

            // Read impatience configuration from NetworkStruct
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int istStation = serviceStations.get(svcIdx);
                jline.lang.nodes.Station station = sn.stations.get(istStation);
                for (int k = 0; k < numClasses; k++) {
                    jline.lang.JobClass jobClass = sn.jobclasses.get(k);

                    // Patience / reneging configuration
                    ProcessType impatienceTypeVal = null;
                    if (sn.impatienceType != null) {
                        Map<jline.lang.JobClass, ProcessType> rmap =
                                sn.impatienceType.get(station);
                        if (rmap != null) impatienceTypeVal = rmap.get(jobClass);
                    }
                    if (impatienceTypeVal != null
                            && impatienceTypeVal != ProcessType.DISABLED) {
                        hasPatienceConfig[svcIdx][k] = true;

                        MRG32k3a stream = new MRG32k3a();
                        if (seed > 0) {
                            long offset = ((long) svcIdx * numClasses + k) * 10 + 50000;
                            stream.setSeed(new long[] { seed + offset, seed + offset + 1,
                                    seed + offset + 2, seed + offset + 3,
                                    seed + offset + 4, seed + offset + 5 });
                        }

                        Matrix muMat = null;
                        if (sn.impatienceMu != null) {
                            Map<jline.lang.JobClass, Matrix> rmap = sn.impatienceMu.get(station);
                            if (rmap != null) muMat = rmap.get(jobClass);
                        }
                        Matrix phiMat = null;
                        if (sn.impatiencePhi != null) {
                            Map<jline.lang.JobClass, Matrix> rmap = sn.impatiencePhi.get(station);
                            if (rmap != null) phiMat = rmap.get(jobClass);
                        }
                        MatrixCell procCell = null;
                        if (sn.impatienceProc != null) {
                            Map<jline.lang.JobClass, MatrixCell> rmap =
                                    sn.impatienceProc.get(station);
                            if (rmap != null) procCell = rmap.get(jobClass);
                        }

                        patienceGens[svcIdx][k] = createImpatienceGenerator(
                                impatienceTypeVal, muMat, phiMat, procCell, stream,
                                svcIdx, k, "patience");
                    }

                    // Balking configuration
                    if (sn.balkingStrategy != null) {
                        Map<jline.lang.JobClass, jline.lang.constant.BalkingStrategy> rmap =
                                sn.balkingStrategy.get(station);
                        if (rmap != null && rmap.get(jobClass) != null) {
                            hasBalkingConfig[svcIdx][k] = true;
                        }
                    }

                    // Retrial configuration
                    ProcessType retrialTypeVal = null;
                    if (sn.retrialType != null) {
                        Map<jline.lang.JobClass, ProcessType> rmap = sn.retrialType.get(station);
                        if (rmap != null) retrialTypeVal = rmap.get(jobClass);
                    }
                    if (retrialTypeVal != null && retrialTypeVal != ProcessType.DISABLED) {
                        hasRetrialConfig[svcIdx][k] = true;

                        MRG32k3a stream = new MRG32k3a();
                        if (seed > 0) {
                            long offset = ((long) svcIdx * numClasses + k) * 10 + 60000;
                            stream.setSeed(new long[] { seed + offset, seed + offset + 1,
                                    seed + offset + 2, seed + offset + 3,
                                    seed + offset + 4, seed + offset + 5 });
                        }

                        Matrix muMat = null;
                        if (sn.retrialMu != null) {
                            Map<jline.lang.JobClass, Matrix> rmap = sn.retrialMu.get(station);
                            if (rmap != null) muMat = rmap.get(jobClass);
                        }
                        Matrix phiMat = null;
                        if (sn.retrialPhi != null) {
                            Map<jline.lang.JobClass, Matrix> rmap = sn.retrialPhi.get(station);
                            if (rmap != null) phiMat = rmap.get(jobClass);
                        }
                        MatrixCell procCell = null;
                        if (sn.retrialProc != null) {
                            Map<jline.lang.JobClass, MatrixCell> rmap =
                                    sn.retrialProc.get(station);
                            if (rmap != null) procCell = rmap.get(jobClass);
                        }

                        retrialGens[svcIdx][k] = createImpatienceGenerator(
                                retrialTypeVal, muMat, phiMat, procCell, stream,
                                svcIdx, k, "retrial");

                        Integer maxAtt = null;
                        if (sn.retrialMaxAttempts != null) {
                            Map<jline.lang.JobClass, Integer> rmap =
                                    sn.retrialMaxAttempts.get(station);
                            if (rmap != null) maxAtt = rmap.get(jobClass);
                        }
                        retrialMaxAttemptsConfig[svcIdx][k] = (maxAtt != null) ? maxAtt : -1;
                    }
                }
            }
        }

        /**
         * Random variate generator for impatience (patience/retrial) distributions
         * (Kotlin lines 5564-5688).
         */
        private RandomVariateGen createImpatienceGenerator(ProcessType procType,
                                                           Matrix muMatrix,
                                                           Matrix phiMatrix,
                                                           MatrixCell proc,
                                                           MRG32k3a stream,
                                                           int svcIdx, int classIdx,
                                                           String label) {
            if (procType == ProcessType.DISABLED) return null;

            double rate = (muMatrix != null) ? muMatrix.get(0, 0) : 0.0;
            double scv = (phiMatrix != null) ? phiMatrix.get(0, 0) : 1.0;

            if (procType == ProcessType.EXP) {
                if (rate > 0 && rate < Double.MAX_VALUE) {
                    return new umontreal.ssj.randvar.ExponentialGen(stream, rate);
                }
                return null;
            } else if (procType == ProcessType.DET) {
                if (rate > 0 && rate < Double.MAX_VALUE) {
                    return new ConstantGen(stream, 1.0 / rate);
                }
                return null;
            } else if (procType == ProcessType.ERLANG) {
                if (proc != null && proc.size() >= 2) {
                    int k = (int) proc.get(0).get(0, 0);
                    double lambda = proc.get(1).get(0, 0);
                    if (k > 0 && lambda > 0) {
                        return new umontreal.ssj.randvar.ErlangGen(stream, k, lambda);
                    }
                }
                return null;
            } else if (procType == ProcessType.HYPEREXP) {
                if (proc != null && proc.size() >= 2) {
                    Matrix d0 = proc.get(0);
                    Matrix d1 = proc.get(1);
                    if (d0 != null && d1 != null) {
                        double lambda1 = -d0.get(0, 0);
                        double lambda2 = (d0.getNumRows() > 1) ? -d0.get(1, 1) : lambda1;
                        double p = (lambda1 > 0) ? d1.get(0, 0) / lambda1 : 0.5;
                        return new HyperExponentialDistGen(stream,
                                new double[] { p, 1.0 - p },
                                new double[] { lambda1, lambda2 });
                    }
                }
                return null;
            } else if (procType == ProcessType.UNIFORM) {
                if (proc != null && proc.size() >= 2) {
                    double minVal = proc.get(0).get(0, 0);
                    double maxVal = proc.get(1).get(0, 0);
                    if (maxVal > minVal) {
                        return new umontreal.ssj.randvar.UniformGen(stream, minVal, maxVal);
                    }
                }
                return null;
            } else if (procType == ProcessType.GAMMA) {
                if (rate > 0 && rate < Double.MAX_VALUE && scv > 0) {
                    double shape = 1.0 / scv;
                    double lambda = rate / scv;
                    return new umontreal.ssj.randvar.GammaGen(stream, shape, lambda);
                }
                return null;
            } else if (procType == ProcessType.PARETO) {
                if (rate > 0 && rate < Double.MAX_VALUE && scv > 0) {
                    double mean = 1.0 / rate;
                    double shape = FastMath.sqrt(1.0 + 1.0 / scv) + 1.0;
                    double scale = mean * (shape - 1.0) / shape;
                    if (shape > 1.0 && scale > 0) {
                        return new umontreal.ssj.randvar.ParetoGen(stream, shape, scale);
                    }
                }
                return null;
            } else if (procType == ProcessType.WEIBULL) {
                if (rate > 0 && rate < Double.MAX_VALUE && scv > 0) {
                    double mean = 1.0 / rate;
                    double c = FastMath.sqrt(scv);
                    double r = FastMath.pow(c, -1.086);
                    double alpha = mean
                            / org.apache.commons.math3.special.Gamma.gamma(1.0 + 1.0 / r);
                    if (r > 0 && alpha > 0) {
                        return new umontreal.ssj.randvar.WeibullGen(stream, r,
                                1.0 / alpha, 0.0);
                    }
                }
                return null;
            } else if (procType == ProcessType.LOGNORMAL) {
                if (rate > 0 && rate < Double.MAX_VALUE && scv > 0) {
                    double mean = 1.0 / rate;
                    double c = FastMath.sqrt(scv);
                    double c2plus1 = c * c + 1.0;
                    double mu = FastMath.log(mean / FastMath.sqrt(c2plus1));
                    double sigma = FastMath.sqrt(FastMath.log(c2plus1));
                    if (sigma > 0) {
                        return new umontreal.ssj.randvar.LognormalGen(stream, mu, sigma);
                    }
                }
                return null;
            } else if (procType == ProcessType.PH || procType == ProcessType.APH
                    || procType == ProcessType.COXIAN || procType == ProcessType.COX2
                    || procType == ProcessType.MAP || procType == ProcessType.MMPP2
                    || procType == ProcessType.ME || procType == ProcessType.RAP) {
                if (proc != null && proc.size() >= 2) {
                    Random rng;
                    if (seed > 0) {
                        long offset = ((long) svcIdx * numClasses + classIdx) * 10 + 70000;
                        rng = new Random(seed + offset);
                    } else {
                        rng = new Random();
                    }
                    if (procType == ProcessType.ME) {
                        return new MeSampleGen(stream, proc, rng);
                    }
                    if (procType == ProcessType.RAP) {
                        return new RapSampleGen(stream, proc, rng);
                    }
                    return new MapSampleGen(stream, proc, rng);
                }
                return null;
            } else {
                if (rate > 0 && rate < Double.MAX_VALUE) {
                    return new umontreal.ssj.randvar.ExponentialGen(stream, rate);
                }
                return null;
            }
        }

        /**
         * Initialize polling scheduling support (Kotlin lines 5694-5750).
         */
        @SuppressWarnings("unchecked")
        private void initializePollingState() {
            isPollingStation = new boolean[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                isPollingStation[svcIdx] = (schedStrategies[svcIdx] == SchedStrategy.POLLING);
            }

            // Detect Batch Markovian Service Process (BMSP) stations: a BMAP
            // assigned as a station's service process turns it into a bulk
            // (batch) server. Supported only in the canonical M/BMSP/1 regime
            // (single server, single class, FCFS); every other combination is
            // rejected here rather than silently producing zero service time.
            isBatchServiceStation = new boolean[numServiceNodes];
            batchServiceClass = new int[numServiceNodes];
            batchServiceSampler = new Map_sample.BmapSampler[numServiceNodes];
            batchServiceRng = new Random[numServiceNodes];
            batchPendingSize = new int[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                batchServiceClass[svcIdx] = -1;
                int bmapClass = -1;
                int enabledCount = 0;
                boolean anyBmap = false;
                for (int k = 0; k < numClasses; k++) {
                    ProcessType pt = serviceProcessType[svcIdx][k];
                    if (pt != ProcessType.DISABLED) {
                        enabledCount++;
                    }
                    if (pt == ProcessType.BMAP) {
                        anyBmap = true;
                        bmapClass = k;
                    }
                }
                if (!anyBmap) {
                    continue;
                }
                int stationIdx = serviceStations.get(svcIdx);
                String stName = sn.stations.get(stationIdx).getName();
                if (isDelayNode.get(svcIdx).booleanValue()) {
                    throw new RuntimeException("LDES: BMAP service (batch Markovian "
                            + "service process) is not supported at the infinite-server "
                            + "station '" + stName + "'. Batch service requires a "
                            + "single-server FCFS queue.");
                }
                if (numServers[svcIdx] != 1) {
                    throw new RuntimeException("LDES: BMAP service (BMSP) at station '"
                            + stName + "' requires a single server; multiserver bulk "
                            + "service is not supported.");
                }
                if (schedStrategies[svcIdx] != SchedStrategy.FCFS) {
                    throw new RuntimeException("LDES: BMAP service (BMSP) at station '"
                            + stName + "' requires FCFS scheduling.");
                }
                if (enabledCount != 1) {
                    throw new RuntimeException("LDES: BMAP service (BMSP) at station '"
                            + stName + "' requires exactly one job class; multiclass "
                            + "bulk service is not supported.");
                }
                if (isPASStation[svcIdx] || isPreemptiveScheduling[svcIdx]
                        || isPollingStation[svcIdx]) {
                    throw new RuntimeException("LDES: BMAP service (BMSP) at station '"
                            + stName + "' is incompatible with PAS/preemptive/polling "
                            + "scheduling.");
                }
                if (hasSetupDelayoff[svcIdx]) {
                    throw new RuntimeException("LDES: BMAP service (BMSP) at station '"
                            + stName + "' is incompatible with server setup/delayoff.");
                }
                if (sn.immfeed != null && stationIdx < sn.immfeed.getNumRows()
                        && bmapClass < sn.immfeed.getNumCols()
                        && sn.immfeed.get(stationIdx, bmapClass) > 0.0) {
                    throw new RuntimeException("LDES: BMAP service (BMSP) at station '"
                            + stName + "' is incompatible with immediate feedback.");
                }
                MatrixCell proc = serviceProc[svcIdx][bmapClass];
                if (proc == null || proc.size() < 3) {
                    throw new RuntimeException("LDES: BMAP service process at station '"
                            + stName + "' is malformed (expected D0, D1_total, D1, ...).");
                }
                isBatchServiceStation[svcIdx] = true;
                batchServiceClass[svcIdx] = bmapClass;
                batchServiceSampler[svcIdx] = new Map_sample.BmapSampler(proc);
                if (seed > 0) {
                    long offset = ((long) svcIdx) * 131L + 700000L;
                    batchServiceRng[svcIdx] = new Random(seed + offset);
                } else {
                    batchServiceRng[svcIdx] = new Random();
                }
                batchPendingSize[svcIdx] = 1;
            }

            pollingType = new PollingType[numServiceNodes];
            pollingK = new int[numServiceNodes];
            pollingCurrentClass = new int[numServiceNodes];
            pollingJobsServedInRound = new int[numServiceNodes];
            pollingGateSize = new int[numServiceNodes];
            pollingInSwitchover = new boolean[numServiceNodes];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                pollingK[svcIdx] = 1;
            }

            // Per-class queues for polling stations
            pollingQueues = (LinkedList<Customer>[][]) new LinkedList<?>[numServiceNodes][numClasses];
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    pollingQueues[svcIdx][k] = new LinkedList<Customer>();
                }
            }

            // Switchover generators
            pollingSwitchoverGens = new RandomVariateGen[numServiceNodes][numClasses];

            // Extract polling parameters from sn.nodeparam
            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                if (!isPollingStation[svcIdx]) continue;

                int nodeIdx = serviceNodes.get(svcIdx);
                jline.lang.nodes.Node node = sn.nodes.get(nodeIdx);

                jline.lang.NodeParam nodeParam =
                        (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
                if (nodeParam instanceof QueueNodeParam) {
                    QueueNodeParam qnp = (QueueNodeParam) nodeParam;
                    pollingType[svcIdx] = (qnp.pollingType != null)
                            ? qnp.pollingType : PollingType.EXHAUSTIVE;
                    pollingK[svcIdx] = (qnp.pollingPar != null) ? qnp.pollingPar : 1;

                    for (int k = 0; k < numClasses; k++) {
                        jline.lang.processes.Distribution switchoverDist =
                                (qnp.switchoverTime != null)
                                        ? qnp.switchoverTime.get(sn.jobclasses.get(k))
                                        : null;
                        // An Immediate switchover is a zero-time leg, not a very fast
                        // one: its rate is GlobalConstants.Immediate (1e8), so drawing
                        // from it would both charge a spurious ~1e-8 delay and, on an
                        // empty station, make the server lap ~1e8 times per unit time.
                        // Leaving the generator null marks the leg as costing no time,
                        // which pollingInitiateSwitchover walks straight through.
                        if (switchoverDist != null && !switchoverDist.isImmediate()) {
                            double rate = switchoverDist.rate();
                            if (rate > 0 && rate < Double.MAX_VALUE) {
                                MRG32k3a stream = new MRG32k3a();
                                stream.setSeed(new long[] {
                                        seed + 5000 + (long) svcIdx * 100L + k,
                                        seed + 5001, seed + 5002,
                                        seed + 5003, seed + 5004, seed + 5005 });
                                pollingSwitchoverGens[svcIdx][k] =
                                        new umontreal.ssj.randvar.ExponentialGen(stream, rate);
                            }
                        }
                    }
                } else {
                    pollingType[svcIdx] = PollingType.EXHAUSTIVE;
                }
            }
        }

        /**
         * Initialize load-dependent service support (Kotlin lines 5757-5789).
         */
        private void initializeLoadDependentService() {
            boolean hasLldConfig = sn.lldscaling != null
                    && !sn.lldscaling.isEmpty()
                    && sn.lldscaling.getNumCols() > 0;

            lldScaling = new double[numServiceNodes][];
            isLoadDependent = new boolean[numServiceNodes];

            if (!hasLldConfig) return;

            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int stationIdx = serviceStations.get(svcIdx);
                int maxN = sn.lldscaling.getNumCols();

                boolean hasScaling = false;
                for (int n = 0; n < maxN; n++) {
                    double scale = sn.lldscaling.get(stationIdx, n);
                    if (scale != 1.0 && scale > 0) {
                        hasScaling = true;
                        break;
                    }
                }

                if (hasScaling) {
                    isLoadDependent[svcIdx] = true;
                    hasLld = true;
                    double[] arr = new double[maxN];
                    for (int n = 0; n < maxN; n++) {
                        double scale = sn.lldscaling.get(stationIdx, n);
                        arr[n] = (scale > 0) ? scale : 1.0;
                    }
                    lldScaling[svcIdx] = arr;
                }
            }

            // A load-dependent (non-PS) station's aggregate service rate changes
            // with its population, so the in-service completion events must be
            // rescheduled whenever the population changes (as for cdscaling). Reuse
            // the departure-event tracking maps for that purpose.
            if (hasLld && sdDepartureEvents == null) {
                // Generic array creation requires an unchecked cast; localize it.
                @SuppressWarnings("unchecked")
                Map<Integer, Event>[] depTmp = (Map<Integer, Event>[]) new Map<?, ?>[numServiceNodes];
                sdDepartureEvents = depTmp;
                @SuppressWarnings("unchecked")
                Map<Integer, Customer>[] inSvcTmp = (Map<Integer, Customer>[]) new Map<?, ?>[numServiceNodes];
                sdInServiceCustomers = inSvcTmp;
                for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                    sdDepartureEvents[svcIdx] = new HashMap<Integer, Event>();
                    sdInServiceCustomers[svcIdx] = new HashMap<Integer, Customer>();
                }
            }
        }

        /**
         * Initialize the per-station class-dependence functions beta_{i,r}(n).
         *
         * <p>A station without class dependence is absent from sn.cdscaling and is
         * left null here; the neutral scaling is applied at the point of use.</p>
         */
        @SuppressWarnings("unchecked")
        private void initializeClassDependence() {
            boolean hasCdMap = sn.cdscaling != null && !sn.cdscaling.isEmpty();
            boolean hasJdMap = sn.jdscaling != null && !sn.jdscaling.isEmpty();
            if (!hasCdMap && !hasJdMap) {
                hasCd = false;
                return;
            }

            cdFunctions = (SerializableFunction<Matrix, Matrix>[]) new SerializableFunction<?, ?>[numServiceNodes];

            for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                int stationIdx = serviceStations.get(svcIdx);
                jline.lang.nodes.Station station = sn.stations.get(stationIdx);
                // Product-form class dependence beta_{i,r}(n) and non-product-form
                // joint dependence eta_i(n) enter the sample-path service rate the
                // same way; when both are declared at a station the effective
                // scaling is their product (eta .* beta), mirroring State.java.
                final SerializableFunction<Matrix, Matrix> beta = hasCdMap ? sn.cdscaling.get(station) : null;
                final SerializableFunction<Matrix, Matrix> eta = hasJdMap ? sn.jdscaling.get(station) : null;
                if (beta != null && eta != null) {
                    hasCd = true;
                    cdFunctions[svcIdx] = new SerializableFunction<Matrix, Matrix>() {
                        public Matrix apply(Matrix ni) {
                            // Elementwise product with scalar broadcasting: either
                            // factor may be 1x1 (shared across classes) or 1xR.
                            Matrix b = beta.apply(ni);
                            Matrix e = eta.apply(ni);
                            int n = Math.max(b.length(), e.length());
                            Matrix out = new Matrix(1, n);
                            for (int j = 0; j < n; j++) {
                                double bv = b.get(b.length() > 1 ? j : 0);
                                double ev = e.get(e.length() > 1 ? j : 0);
                                out.set(0, j, bv * ev);
                            }
                            return out;
                        }
                    };
                } else if (beta != null) {
                    hasCd = true;
                    cdFunctions[svcIdx] = beta;
                } else if (eta != null) {
                    hasCd = true;
                    cdFunctions[svcIdx] = eta;
                }
            }

            // Departure event tracking for dynamic re-scaling: the aggregate rate
            // changes with the population, so in-service completions must be
            // rescheduled whenever it changes.
            if (hasCd && sdDepartureEvents == null) {
                // Generic array creation requires an unchecked cast; localize it.
                @SuppressWarnings("unchecked")
                Map<Integer, Event>[] depTmp = (Map<Integer, Event>[]) new Map<?, ?>[numServiceNodes];
                sdDepartureEvents = depTmp;
                @SuppressWarnings("unchecked")
                Map<Integer, Customer>[] inSvcTmp = (Map<Integer, Customer>[]) new Map<?, ?>[numServiceNodes];
                sdInServiceCustomers = inSvcTmp;
                for (int svcIdx = 0; svcIdx < numServiceNodes; svcIdx++) {
                    sdDepartureEvents[svcIdx] = new HashMap<Integer, Event>();
                    sdInServiceCustomers[svcIdx] = new HashMap<Integer, Customer>();
                }
            }
        }

        /**
         * Re-scale service times for in-service jobs at a state-dependent station,
         * whose aggregate rate changes with the population.
         */
        private void rescaleStateDepInServiceJobs(int queueIdx) {
            if (sdDepartureEvents == null || sdInServiceCustomers == null) return;

            Map<Integer, Event> eventsMap = sdDepartureEvents[queueIdx];
            Map<Integer, Customer> customersMap = sdInServiceCustomers[queueIdx];

            if (eventsMap == null || eventsMap.isEmpty()) return;

            List<Integer> serverIds = new ArrayList<Integer>(eventsMap.keySet());
            for (Integer serverId : serverIds) {
                Event oldEvent = eventsMap.get(serverId);
                if (oldEvent == null) continue;
                Customer customer = customersMap.get(serverId);
                if (customer == null) continue;
                int classId = customer.classId;

                oldEvent.cancel();

                // Resample the (memoryless) completion at the current population.
                // Load-dependent stations use the standard generator, which
                // re-applies the load-dependent factor alpha(n) at the present
                // population; class-dependent stations use the heterogeneous generator.
                double newServiceTime;
                if (isLoadDependent[queueIdx] && !hasCd) {
                    newServiceTime = generateServiceTime(queueIdx, classId);
                } else {
                    newServiceTime = generateHeteroServiceTime(queueIdx, classId,
                            customer.assignedServerType);
                }

                Event newDepartureEvent = new Departure(queueIdx, serverId, customer);
                newDepartureEvent.schedule(serviceWallDelay(queueIdx, classId, newServiceTime));

                eventsMap.put(serverId, newDepartureEvent);

                if (hasRemovalSignals) {
                    inServiceJobs.put(new IntPair(queueIdx, serverId),
                            new InServiceJob(customer, newDepartureEvent));
                }
            }
        }

        /**
         * Initialize trace logging if verbose level is DEBUG (Kotlin lines 5937-5947).
         */
        private void initializeTracing() {
            if (options.verbose == VerboseLevel.DEBUG) {
                traceEnabled = true;
                long timestamp = System.currentTimeMillis();
                String tempDir;
                try {
                    tempDir = lineTempName("ssj");
                } catch (java.io.IOException ex) {
                    throw new RuntimeException("LDES: Cannot allocate temp directory for trace", ex);
                }
                File traceFile = new File(tempDir, "ssj_trace_" + timestamp + ".csv");
                try {
                    traceWriter = new PrintWriter(traceFile);
                } catch (java.io.FileNotFoundException ex) {
                    throw new RuntimeException(
                            "LDES: Cannot create trace file " + traceFile, ex);
                }
                traceWriter.println(
                        "time,event_type,station_idx,class_id,queue_length,busy_servers,phase");
            }
        }

        /**
         * Log a simulation event to the trace file (Kotlin lines 5958-5963).
         */
        private void logEvent(String eventType, int stationIdx, int classId,
                              int queueLength, int busyServers) {
            if (traceEnabled && traceWriter != null) {
                String phase = warmupDone ? "STEADY_STATE" : "WARMUP";
                traceWriter.println(ssjSim.time() + "," + eventType + "," + stationIdx
                        + "," + classId + "," + queueLength + "," + busyServers
                        + "," + phase);
            }
        }

        /**
         * Close the trace writer (Kotlin lines 5968-5974).
         */
        private void closeTracing() {
            if (traceWriter != null) {
                traceWriter.flush();
                traceWriter.close();
                traceWriter = null;
            }
        }

        // ==================== Logger Node Support ====================

        /**
         * Initialize Logger nodes (Kotlin lines 5982- continues into PART 5).
         */
        private void initializeLoggers() {
            if (loggerNodes.isEmpty()) return;

            java.text.SimpleDateFormat dateFormat =
                    new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
            simulationStartTime = dateFormat.format(new java.util.Date());

            for (Integer loggerNodeIdx : loggerNodes) {
                jline.lang.nodes.Node loggerNode = sn.nodes.get(loggerNodeIdx);
                if (!(loggerNode instanceof jline.lang.nodes.Logger)) continue;

                jline.lang.nodes.Logger logger = (jline.lang.nodes.Logger) loggerNode;
                String fileName = (logger.fileName != null)
                        ? logger.fileName : ("logger_" + logger.getName() + ".csv");
                String filePath = (logger.filePath != null) ? logger.filePath : "";
                String loggerName = (logger.getName() != null) ? logger.getName() : "Logger";

                // The remainder of the LoggerConfig construction (boolean flags,
                // PrintWriter setup, header writing) is translated in PART 5
                // (Kotlin lines 6001+).
                initializeLoggerConfigForNode(loggerNodeIdx, fileName, filePath, loggerName);
            }
        }

        /**
         * Builds a LoggerConfig record and writes the CSV header for a Logger
         * node (Kotlin lines 5998-6035, "remainder of {@code initializeLoggers}").
         */
        private void initializeLoggerConfigForNode(int loggerNodeIdx, String fileName,
                                                   String filePath, String loggerName) {
            jline.lang.nodes.Node loggerNode = sn.nodes.get(loggerNodeIdx);
            if (!(loggerNode instanceof jline.lang.nodes.Logger)) return;
            jline.lang.nodes.Logger logger = (jline.lang.nodes.Logger) loggerNode;

            LoggerConfig config = new LoggerConfig(
                    loggerNodeIdx,
                    fileName,
                    filePath,
                    logger.getLoggerName(),
                    logger.getTimestamp(),
                    logger.getJobID(),
                    logger.getJobClass(),
                    logger.getTimeSameClass(),
                    logger.getTimeAnyClass(),
                    logger.getStartTime(),
                    loggerName);
            loggerConfigs.put(loggerNodeIdx, config);

            // Initialize last job time tracking
            double[] perClass = new double[numClasses];
            for (int k = 0; k < numClasses; k++) perClass[k] = 0.0;
            loggerLastJobTimePerClass.put(loggerNodeIdx, perClass);
            loggerLastJobTimeAny.put(loggerNodeIdx, Double.valueOf(0.0));

            // Create output file
            try {
                File outputPath;
                if (filePath != null && filePath.length() > 0) {
                    outputPath = new File(filePath, fileName);
                } else {
                    outputPath = new File(fileName);
                }
                File parent = outputPath.getParentFile();
                if (parent != null) parent.mkdirs();

                BufferedWriter writer =
                        new BufferedWriter(new FileWriter(outputPath, false));
                loggerWriters.put(loggerNodeIdx, writer);

                writer.write("LOGGERNAME,TIMESTAMP,JOB_ID,CLASS_ID,INTERARRIVAL_SAMECLASS,INTERARRIVAL_ANYCLASS,SIMUL_START_TIME");
                writer.newLine();
            } catch (Exception e) {
                line_warning("solver_ssj",
                        "Failed to create logger file for %s: %s",
                        loggerName, (e.getMessage() == null ? "" : e.getMessage()));
            }
        }

        // ==================================================================
        // PART 5 TRANSLATION (Kotlin lines 6001-7500).
        //
        // Routing infrastructure (initialization, destination selection for
        // PROB / RAND / RROBIN / WRROBIN / JSQ / SQ strategies, and
        // routing through pass-through nodes), Logger CSV writing,
        // PS-family rate calculations (PS, DPS, GPS, PSPRIO, DPSPRIO,
        // GPSPRIO), the {@code arriveAtPSQueue} entry point, and the
        // {@code PSDeparture} inner event class.  Forward-reference stubs
        // for methods called by these (resolved in PART 6+) are included
        // at the bottom of this block.
        // ==================================================================

        /**
         * Initialize routing strategies for each node and class
         * (Kotlin lines 6044-6140).
         */
        @SuppressWarnings("unchecked")
        private void initializeRouting() {
            int I = numNodes;
            int R = numClasses;

            // Initialize routing strategies from sn.routing
            nodeRoutingStrategies = new RoutingStrategy[I][R];
            for (int nodeIdx = 0; nodeIdx < I; nodeIdx++) {
                jline.lang.nodes.Node node = (nodeIdx < sn.nodes.size())
                        ? sn.nodes.get(nodeIdx) : null;
                Map<jline.lang.JobClass, RoutingStrategy> classRouting =
                        (node != null && sn.routing != null && sn.routing.containsKey(node))
                                ? sn.routing.get(node) : null;
                for (int classIdx = 0; classIdx < R; classIdx++) {
                    jline.lang.JobClass jobClass =
                            (classIdx < sn.jobclasses.size()) ? sn.jobclasses.get(classIdx) : null;
                    if (classRouting != null && jobClass != null) {
                        nodeRoutingStrategies[nodeIdx][classIdx] = classRouting.get(jobClass);
                    } else {
                        nodeRoutingStrategies[nodeIdx][classIdx] = null;
                    }
                }
            }

            // Initialize round-robin counters for each node
            roundRobinCounters = new int[I][R];

            // Initialize WRROBIN weights and RROBIN/WRROBIN outlinks from nodeparam
            wrrobinWeights = new double[I][R][];
            rroutlinks = new int[I][R][];
            // Initialize per-class SQ K (default to legacy fallback)
            sqDByNodeClass = new int[I][R];
            for (int i = 0; i < I; i++) {
                for (int r = 0; r < R; r++) {
                    sqDByNodeClass[i][r] = sqD;
                }
            }

            for (int nodeIdx = 0; nodeIdx < I; nodeIdx++) {
                jline.lang.nodes.Node node = (nodeIdx < sn.nodes.size())
                        ? sn.nodes.get(nodeIdx) : null;
                if (node == null) continue;
                NodeParam param = (sn.nodeparam != null) ? sn.nodeparam.get(node) : null;
                if (param == null) continue;

                for (int classIdx = 0; classIdx < R; classIdx++) {
                    jline.lang.JobClass jobClass =
                            (classIdx < sn.jobclasses.size()) ? sn.jobclasses.get(classIdx) : null;
                    if (jobClass == null) continue;
                    RoutingStrategy routingStrategy = nodeRoutingStrategies[nodeIdx][classIdx];

                    // Get outlinks for RROBIN/WRROBIN
                    if (routingStrategy == RoutingStrategy.RROBIN
                            || routingStrategy == RoutingStrategy.WRROBIN) {
                        Matrix outlinkMatrix =
                                (param.outlinks != null) ? param.outlinks.get(jobClass) : null;
                        if (outlinkMatrix != null && outlinkMatrix.length() > 0) {
                            int len = (int) outlinkMatrix.length();
                            int[] outlinks = new int[len];
                            for (int i = 0; i < len; i++) {
                                outlinks[i] = (int) outlinkMatrix.get(i);
                            }
                            rroutlinks[nodeIdx][classIdx] = outlinks;
                        }
                    }

                    // Get d for SQ (per-class)
                    if (routingStrategy == RoutingStrategy.SQ) {
                        Map<jline.lang.JobClass, Integer> kMap = param.d;
                        if (kMap != null && kMap.containsKey(jobClass)) {
                            sqDByNodeClass[nodeIdx][classIdx] = kMap.get(jobClass).intValue();
                        }
                    }

                    // Get weights for WRROBIN
                    if (routingStrategy == RoutingStrategy.WRROBIN) {
                        Matrix weightMatrix =
                                (param.weights != null) ? param.weights.get(jobClass) : null;
                        if (weightMatrix != null) {
                            double[] weights = new double[I];
                            int numCols = weightMatrix.getNumCols();
                            for (int destIdx = 0; destIdx < I; destIdx++) {
                                weights[destIdx] = (destIdx < numCols)
                                        ? weightMatrix.get(0, destIdx) : 0.0;
                            }
                            wrrobinWeights[nodeIdx][classIdx] = weights;
                        }
                    }
                }
            }

            // Initialize class switch matrices (already encoded in routing
            // matrix by LINE; retained for potential future use).
            classSwitchMatrices = new double[I][][];
            for (Integer csNodeIdxInt : classSwitchNodes) {
                int csNodeIdx = csNodeIdxInt.intValue();
                jline.lang.nodes.Node node = (csNodeIdx < sn.nodes.size())
                        ? sn.nodes.get(csNodeIdx) : null;
                if (!(node instanceof jline.lang.nodes.ClassSwitch)) continue;
                jline.lang.nodes.ClassSwitch csNode = (jline.lang.nodes.ClassSwitch) node;
                jline.lang.sections.ServiceSection serverSec = csNode.getServer();
                if (!(serverSec instanceof jline.lang.sections.ClassSwitcher)) continue;
                jline.lang.sections.ClassSwitcher csServer =
                        (jline.lang.sections.ClassSwitcher) serverSec;

                double[][] csMatrix = new double[R][R];
                for (int fromClass = 0; fromClass < R; fromClass++) {
                    for (int toClass = 0; toClass < R; toClass++) {
                        csMatrix[fromClass][toClass] = csServer.applyCsFun(fromClass, toClass);
                    }
                }
                classSwitchMatrices[csNodeIdx] = csMatrix;
            }
        }

        /**
         * Log a job passage through a Logger node (Kotlin lines 6150-6212).
         */
        private void logJobPassage(int loggerNodeIdx, int classId, long jobId) {
            LoggerConfig config = loggerConfigs.get(loggerNodeIdx);
            if (config == null) return;
            BufferedWriter writer = loggerWriters.get(loggerNodeIdx);
            if (writer == null) return;

            double currentTime = ssjSim.time();
            String delimiter = config.delimiter;

            // Calculate interarrival times
            double[] lastTimePerClass = loggerLastJobTimePerClass.get(loggerNodeIdx);
            Double lastTimeAnyBoxed = loggerLastJobTimeAny.get(loggerNodeIdx);
            double lastTimeAny = (lastTimeAnyBoxed != null) ? lastTimeAnyBoxed.doubleValue() : 0.0;

            double prevSameClass = (lastTimePerClass != null
                    && classId >= 0 && classId < lastTimePerClass.length)
                    ? lastTimePerClass[classId] : 0.0;
            double interarrivalSameClass = currentTime - prevSameClass;
            double interarrivalAnyClass = currentTime - lastTimeAny;

            // Update last job times
            if (lastTimePerClass != null && classId >= 0 && classId < lastTimePerClass.length) {
                lastTimePerClass[classId] = currentTime;
            }
            loggerLastJobTimeAny.put(loggerNodeIdx, Double.valueOf(currentTime));

            // Get class name
            String className;
            if (classId >= 0 && classId < sn.jobclasses.size()) {
                className = sn.jobclasses.get(classId).getName();
            } else {
                className = Integer.toString(classId);
            }

            StringBuilder row = new StringBuilder();

            // LOGGERNAME
            if (config.logLoggerName) row.append(config.loggerName);
            row.append(delimiter);

            // TIMESTAMP
            if (config.logTimestamp) row.append(formatNumber(currentTime, config.decimalSeparator));
            row.append(delimiter);

            // JOB_ID
            if (config.logJobID) row.append(jobId);
            row.append(delimiter);

            // CLASS_ID
            if (config.logJobClass) row.append(className);
            row.append(delimiter);

            // INTERARRIVAL_SAMECLASS
            if (config.logTimeSameClass)
                row.append(formatNumber(interarrivalSameClass, config.decimalSeparator));
            row.append(delimiter);

            // INTERARRIVAL_ANYCLASS
            if (config.logTimeAnyClass)
                row.append(formatNumber(interarrivalAnyClass, config.decimalSeparator));
            row.append(delimiter);

            // SIMUL_START_TIME
            if (config.logStartTime) row.append(simulationStartTime);

            try {
                writer.write(row.toString());
                writer.newLine();
            } catch (Exception e) {
                // Silently ignore write errors to avoid disrupting simulation
            }
        }

        /**
         * Format a number with the specified decimal separator
         * (Kotlin lines 6217-6224).
         */
        private String formatNumber(double value, String decimalSeparator) {
            String formatted = Double.toString(value);
            if (!".".equals(decimalSeparator)) {
                return formatted.replace(".", decimalSeparator);
            }
            return formatted;
        }

        /**
         * Close all Logger file writers (Kotlin lines 6229-6239).
         */
        private void closeLoggers() {
            for (Map.Entry<Integer, BufferedWriter> e : loggerWriters.entrySet()) {
                BufferedWriter writer = e.getValue();
                try {
                    writer.flush();
                    writer.close();
                } catch (Exception ex) {
                    // Ignore close errors
                }
            }
            loggerWriters.clear();
        }

        /**
         * Route a job through pass-through nodes (Logger, Router,
         * ClassSwitch, Cache) until reaching a service node or sink
         * (Kotlin lines 6253-6296).
         */
        private RoutingResult routeThroughPassthroughNodes(int fromNode, int classId, long jobId) {
            int currentNode = fromNode;
            int currentClass = classId;
            int maxIterations = 100;

            while (maxIterations > 0) {
                RoutingResult result = selectDestinationWithClassSwitch(currentNode, currentClass);

                if (result.destNode < 0) {
                    return new RoutingResult(-1, currentClass);
                }

                if (loggerNodes.contains(Integer.valueOf(result.destNode))) {
                    logJobPassage(result.destNode, result.destClassId, jobId);
                    currentNode = result.destNode;
                    currentClass = result.destClassId;
                    maxIterations--;
                } else if (routerNodes.contains(Integer.valueOf(result.destNode))) {
                    currentNode = result.destNode;
                    currentClass = result.destClassId;
                    maxIterations--;
                } else if (classSwitchNodes.contains(Integer.valueOf(result.destNode))) {
                    currentNode = result.destNode;
                    currentClass = result.destClassId;
                    maxIterations--;
                } else if (cacheNodes.contains(Integer.valueOf(result.destNode))) {
                    int newClass = processCacheAccess(result.destNode, result.destClassId, jobId);
                    if (newClass == CACHE_HELD) {
                        // Request parked as a delayed hit; it will be delivered when the
                        // in-flight fetch completes. Signal "no destination" so the caller
                        // does not route it now.
                        return new RoutingResult(-1, result.destClassId);
                    }
                    currentNode = result.destNode;
                    currentClass = newClass;
                    maxIterations--;
                } else {
                    return result;
                }
            }
            return new RoutingResult(-1, currentClass);
        }

        /**
         * Apply class switch at a ClassSwitch node (Kotlin lines 6306-6325).
         */
        private int applyClassSwitch(int csNodeIdx, int inClass) {
            double[][] csMatrix = classSwitchMatrices[csNodeIdx];
            if (csMatrix == null) return inClass;

            double[] probs = csMatrix[inClass];
            double cumProb = 0.0;
            for (int i = 0; i < probs.length; i++) cumProb += probs[i];
            if (cumProb <= 0) return inClass;

            double rand = routingRng.nextDouble() * cumProb;
            double cumulative = 0.0;
            for (int toClass = 0; toClass < probs.length; toClass++) {
                cumulative += probs[toClass];
                if (rand <= cumulative) {
                    return toClass;
                }
            }
            return inClass;
        }

        /**
         * Select destination node and class, supporting class switching via
         * routing matrix (Kotlin lines 6332-6378).
         */
        private RoutingResult selectDestinationWithClassSwitch(int fromNode, int classId) {
            // Self-looping classes route to reference station node
            if (sn.isslc != null && sn.isslc.get(classId) == 1.0) {
                int refStationIdx = referenceStation[classId];
                return new RoutingResult((int) sn.stationToNode.get(refStationIdx), classId);
            }

            int R = numClasses;
            int I = numNodes;

            List<RoutingResult> destinations = new ArrayList<RoutingResult>();
            List<Double> destProbs = new ArrayList<Double>();

            for (int toNode = 0; toNode < I; toNode++) {
                for (int toClass = 0; toClass < R; toClass++) {
                    double prob = sn.rtnodes.get(fromNode * R + classId, toNode * R + toClass);
                    if (prob > 0) {
                        destinations.add(new RoutingResult(toNode, toClass));
                        destProbs.add(Double.valueOf(prob));
                    }
                }
            }

            if (destinations.isEmpty()) {
                return new RoutingResult(-1, classId);
            }
            if (destinations.size() == 1) {
                return destinations.get(0);
            }

            RoutingStrategy routingStrategy = null;
            if (nodeRoutingStrategies != null
                    && fromNode >= 0 && fromNode < nodeRoutingStrategies.length
                    && classId >= 0 && classId < nodeRoutingStrategies[fromNode].length) {
                routingStrategy = nodeRoutingStrategies[fromNode][classId];
            }

            if (routingStrategy == RoutingStrategy.RAND) {
                return selectRandomDestinationWithClass(destinations);
            } else if (routingStrategy == RoutingStrategy.RROBIN) {
                return selectRoundRobinDestinationWithClass(fromNode, classId, destinations);
            } else if (routingStrategy == RoutingStrategy.WRROBIN) {
                return selectWeightedRoundRobinDestinationWithClass(fromNode, classId, destinations);
            } else if (routingStrategy == RoutingStrategy.JSQ) {
                return selectJSQDestinationWithClass(destinations);
            } else if (routingStrategy == RoutingStrategy.SQ) {
                return selectSQDestinationWithClass(fromNode, classId, destinations);
            } else {
                return selectProbabilisticDestinationWithClass(destinations, destProbs);
            }
        }

        /**
         * Select destination node based on routing strategy (no class switching)
         * (Kotlin lines 6386-6428).
         */
        private int selectDestinationDirect(int fromNode, int classId) {
            if (sn.isslc != null && sn.isslc.get(classId) == 1.0) {
                int refStationIdx = referenceStation[classId];
                return (int) sn.stationToNode.get(refStationIdx);
            }

            int R = numClasses;
            int I = numNodes;

            List<Integer> destinations = new ArrayList<Integer>();
            List<Double> destProbs = new ArrayList<Double>();

            for (int toNode = 0; toNode < I; toNode++) {
                double prob = sn.rtnodes.get(fromNode * R + classId, toNode * R + classId);
                if (prob > 0) {
                    destinations.add(Integer.valueOf(toNode));
                    destProbs.add(Double.valueOf(prob));
                }
            }

            if (destinations.isEmpty()) return -1;
            if (destinations.size() == 1) return destinations.get(0).intValue();

            RoutingStrategy routingStrategy = null;
            if (nodeRoutingStrategies != null
                    && fromNode >= 0 && fromNode < nodeRoutingStrategies.length
                    && classId >= 0 && classId < nodeRoutingStrategies[fromNode].length) {
                routingStrategy = nodeRoutingStrategies[fromNode][classId];
            }

            if (routingStrategy == RoutingStrategy.RAND) {
                return selectRandomDestination(destinations);
            } else if (routingStrategy == RoutingStrategy.RROBIN) {
                return selectRoundRobinDestination(fromNode, classId, destinations);
            } else if (routingStrategy == RoutingStrategy.WRROBIN) {
                return selectWeightedRoundRobinDestination(fromNode, classId, destinations);
            } else if (routingStrategy == RoutingStrategy.JSQ) {
                return selectJSQDestination(destinations);
            } else if (routingStrategy == RoutingStrategy.SQ) {
                return selectSQDestination(fromNode, classId, destinations);
            } else {
                return selectProbabilisticDestination(destinations, destProbs);
            }
        }

        /**
         * Probabilistic routing with class (Kotlin lines 6433-6448).
         */
        private RoutingResult selectProbabilisticDestinationWithClass(
                List<RoutingResult> destinations, List<Double> probs) {
            double cumProb = 0.0;
            for (int i = 0; i < probs.size(); i++) cumProb += probs.get(i).doubleValue();
            if (cumProb <= 0) return destinations.get(0);

            double rand = routingRng.nextDouble() * cumProb;
            double cumulative = 0.0;
            for (int i = 0; i < destinations.size(); i++) {
                cumulative += probs.get(i).doubleValue();
                if (rand <= cumulative) return destinations.get(i);
            }
            return destinations.get(destinations.size() - 1);
        }

        /**
         * Random routing with class (Kotlin lines 6453-6456).
         */
        private RoutingResult selectRandomDestinationWithClass(List<RoutingResult> destinations) {
            int idx = (int) (routingRng.nextDouble() * destinations.size());
            if (idx < 0) idx = 0;
            if (idx > destinations.size() - 1) idx = destinations.size() - 1;
            return destinations.get(idx);
        }

        /**
         * Round-robin routing with class (Kotlin lines 6461-6464).
         */
        private RoutingResult selectRoundRobinDestinationWithClass(
                int fromNode, int classId, List<RoutingResult> destinations) {
            int counter = roundRobinCounters[fromNode][classId]++;
            int idx = counter % destinations.size();
            return destinations.get(idx);
        }

        /**
         * Weighted round-robin routing with class (Kotlin lines 6471-6523).
         */
        private RoutingResult selectWeightedRoundRobinDestinationWithClass(
                int fromNode, int classId, List<RoutingResult> destinations) {
            double[] weights = null;
            if (wrrobinWeights != null
                    && fromNode >= 0 && fromNode < wrrobinWeights.length
                    && classId >= 0 && classId < wrrobinWeights[fromNode].length) {
                weights = wrrobinWeights[fromNode][classId];
            }

            if (weights == null || destinations.size() <= 1) {
                int counter = roundRobinCounters[fromNode][classId]++;
                int idx = counter % destinations.size();
                return destinations.get(idx);
            }

            double[] rawWeights = new double[destinations.size()];
            for (int i = 0; i < destinations.size(); i++) {
                int destNode = destinations.get(i).destNode;
                rawWeights[i] = (destNode >= 0 && destNode < weights.length)
                        ? weights[destNode] : 0.0;
            }

            boolean allIntegers = true;
            for (int i = 0; i < rawWeights.length; i++) {
                double w = rawWeights[i];
                if (w != Math.floor(w) || w < 0) { allIntegers = false; break; }
            }

            int[] intWeights = new int[rawWeights.length];
            if (allIntegers) {
                for (int i = 0; i < rawWeights.length; i++) {
                    double w = rawWeights[i];
                    int v = (int) w;
                    if (w > 0 && v < 1) v = 1;
                    intWeights[i] = v;
                }
            } else {
                for (int i = 0; i < rawWeights.length; i++) {
                    double w = rawWeights[i];
                    int v = (int) (w * 1000);
                    if (w > 0 && v < 1) v = 1;
                    intWeights[i] = v;
                }
            }

            int totalWeight = 0;
            for (int i = 0; i < intWeights.length; i++) totalWeight += intWeights[i];
            if (totalWeight <= 0) {
                int counter = roundRobinCounters[fromNode][classId]++;
                int idx = counter % destinations.size();
                return destinations.get(idx);
            }

            int counter = roundRobinCounters[fromNode][classId] % totalWeight;
            roundRobinCounters[fromNode][classId]++;

            int cumulative = 0;
            for (int i = 0; i < destinations.size(); i++) {
                cumulative += intWeights[i];
                if (counter < cumulative) return destinations.get(i);
            }
            return destinations.get(destinations.size() - 1);
        }

        /**
         * Join-Shortest-Queue routing with class (Kotlin lines 6528-6556).
         */
        private RoutingResult selectJSQDestinationWithClass(List<RoutingResult> destinations) {
            int minQueueLength = Integer.MAX_VALUE;
            List<RoutingResult> tiedDests = new ArrayList<RoutingResult>();

            for (RoutingResult dest : destinations) {
                int destNode = dest.destNode;
                if (destNode < 0 || destNode >= sn.nodeToStation.length()) continue;
                int stationIdx = (int) sn.nodeToStation.get(destNode);
                if (stationIdx < 0) continue;

                int svcIdx = serviceStations.indexOf(Integer.valueOf(stationIdx));
                if (svcIdx < 0) continue;

                int queueLen = 0;
                int[] qarr = currentQueueLength[svcIdx];
                for (int k = 0; k < qarr.length; k++) queueLen += qarr[k];
                if (queueLen < minQueueLength) {
                    minQueueLength = queueLen;
                    tiedDests.clear();
                    tiedDests.add(dest);
                } else if (queueLen == minQueueLength) {
                    tiedDests.add(dest);
                }
            }

            if (tiedDests.isEmpty()) return destinations.get(0);
            if (tiedDests.size() == 1) return tiedDests.get(0);
            int idx = (int) (routingRng.nextDouble() * tiedDests.size());
            if (idx < 0) idx = 0;
            if (idx > tiedDests.size() - 1) idx = tiedDests.size() - 1;
            return tiedDests.get(idx);
        }

        /**
         * Power-of-K-Choices routing with class (Kotlin lines 6561-6621).
         */
        private RoutingResult selectSQDestinationWithClass(
                int fromNode, int classId, List<RoutingResult> destinations) {
            int perClassK;
            if (sqDByNodeClass != null
                    && fromNode >= 0 && fromNode < sqDByNodeClass.length
                    && classId >= 0 && classId < sqDByNodeClass[fromNode].length) {
                perClassK = sqDByNodeClass[fromNode][classId];
            } else {
                perClassK = sqD;
            }
            int k = Math.min(perClassK, destinations.size());

            List<RoutingResult> candidates = new ArrayList<RoutingResult>();
            for (int i = 0; i < k; i++) {
                int idx = (int) (routingRng.nextDouble() * destinations.size());
                if (idx < 0) idx = 0;
                if (idx > destinations.size() - 1) idx = destinations.size() - 1;
                candidates.add(destinations.get(idx));
            }

            int minQueueLength = Integer.MAX_VALUE;
            RoutingResult bestDest = candidates.get(0);

            for (RoutingResult dest : candidates) {
                int destNode = dest.destNode;
                if (destNode < 0 || destNode >= sn.nodeToStation.length()) continue;
                int stationIdx = (int) sn.nodeToStation.get(destNode);
                if (stationIdx < 0) continue;

                int svcIdx = serviceStations.indexOf(Integer.valueOf(stationIdx));
                if (svcIdx < 0) continue;

                int queueLen = 0;
                int[] qarr = currentQueueLength[svcIdx];
                for (int kk = 0; kk < qarr.length; kk++) queueLen += qarr[kk];
                if (queueLen < minQueueLength) {
                    minQueueLength = queueLen;
                    bestDest = dest;
                }
            }

            return bestDest;
        }

        /**
         * Probabilistic routing (Kotlin lines 6627-6642).
         */
        private int selectProbabilisticDestination(List<Integer> destinations, List<Double> probs) {
            double cumProb = 0.0;
            for (int i = 0; i < probs.size(); i++) cumProb += probs.get(i).doubleValue();
            if (cumProb <= 0) return destinations.get(0).intValue();

            double rand = routingRng.nextDouble() * cumProb;
            double cumulative = 0.0;
            for (int i = 0; i < destinations.size(); i++) {
                cumulative += probs.get(i).doubleValue();
                if (rand <= cumulative) return destinations.get(i).intValue();
            }
            return destinations.get(destinations.size() - 1).intValue();
        }

        /**
         * Random uniform routing (Kotlin lines 6647-6650).
         */
        private int selectRandomDestination(List<Integer> destinations) {
            int idx = (int) (routingRng.nextDouble() * destinations.size());
            if (idx < 0) idx = 0;
            if (idx > destinations.size() - 1) idx = destinations.size() - 1;
            return destinations.get(idx).intValue();
        }

        /**
         * Round-robin routing (Kotlin lines 6656-6660).
         */
        private int selectRoundRobinDestination(int fromNode, int classId,
                                                List<Integer> destinations) {
            int idx = roundRobinCounters[fromNode][classId] % destinations.size();
            roundRobinCounters[fromNode][classId]++;
            return destinations.get(idx).intValue();
        }

        /**
         * Weighted round-robin routing (Kotlin lines 6667-6721).
         */
        private int selectWeightedRoundRobinDestination(int fromNode, int classId,
                                                        List<Integer> destinations) {
            double[] weights = null;
            if (wrrobinWeights != null
                    && fromNode >= 0 && fromNode < wrrobinWeights.length
                    && classId >= 0 && classId < wrrobinWeights[fromNode].length) {
                weights = wrrobinWeights[fromNode][classId];
            }

            if (weights == null || destinations.size() <= 1) {
                int idx = roundRobinCounters[fromNode][classId] % destinations.size();
                roundRobinCounters[fromNode][classId]++;
                return destinations.get(idx).intValue();
            }

            double[] rawWeights = new double[destinations.size()];
            for (int i = 0; i < destinations.size(); i++) {
                int destNode = destinations.get(i).intValue();
                rawWeights[i] = (destNode >= 0 && destNode < weights.length)
                        ? weights[destNode] : 0.0;
            }

            boolean allIntegers = true;
            for (int i = 0; i < rawWeights.length; i++) {
                double w = rawWeights[i];
                if (w != Math.floor(w) || w < 0) { allIntegers = false; break; }
            }

            int[] intWeights = new int[rawWeights.length];
            if (allIntegers) {
                for (int i = 0; i < rawWeights.length; i++) {
                    double w = rawWeights[i];
                    int v = (int) w;
                    if (w > 0 && v < 1) v = 1;
                    intWeights[i] = v;
                }
            } else {
                for (int i = 0; i < rawWeights.length; i++) {
                    double w = rawWeights[i];
                    int v = (int) (w * 1000);
                    if (w > 0 && v < 1) v = 1;
                    intWeights[i] = v;
                }
            }

            int totalWeight = 0;
            for (int i = 0; i < intWeights.length; i++) totalWeight += intWeights[i];
            if (totalWeight <= 0) {
                int idx = roundRobinCounters[fromNode][classId] % destinations.size();
                roundRobinCounters[fromNode][classId]++;
                return destinations.get(idx).intValue();
            }

            int counter = roundRobinCounters[fromNode][classId] % totalWeight;
            roundRobinCounters[fromNode][classId]++;

            int cumulative = 0;
            for (int i = 0; i < destinations.size(); i++) {
                cumulative += intWeights[i];
                if (counter < cumulative) return destinations.get(i).intValue();
            }
            return destinations.get(destinations.size() - 1).intValue();
        }

        /**
         * Join-Shortest-Queue routing (Kotlin lines 6727-6757).
         */
        private int selectJSQDestination(List<Integer> destinations) {
            int minQueue = Integer.MAX_VALUE;
            List<Integer> tiedDests = new ArrayList<Integer>();

            for (Integer destNodeBoxed : destinations) {
                int destNode = destNodeBoxed.intValue();
                if (destNode < 0 || destNode >= sn.nodeToStation.length()) continue;
                int stationIdx = (int) sn.nodeToStation.get(destNode);
                if (stationIdx < 0) continue;

                int svcIdx = serviceStations.indexOf(Integer.valueOf(stationIdx));
                if (svcIdx >= 0) {
                    int totalQueue = 0;
                    int[] qarr = currentQueueLength[svcIdx];
                    for (int k = 0; k < qarr.length; k++) totalQueue += qarr[k];
                    if (totalQueue < minQueue) {
                        minQueue = totalQueue;
                        tiedDests.clear();
                        tiedDests.add(Integer.valueOf(destNode));
                    } else if (totalQueue == minQueue) {
                        tiedDests.add(Integer.valueOf(destNode));
                    }
                }
            }

            if (tiedDests.isEmpty()) return destinations.get(0).intValue();
            if (tiedDests.size() == 1) return tiedDests.get(0).intValue();
            int idx = (int) (routingRng.nextDouble() * tiedDests.size());
            if (idx < 0) idx = 0;
            if (idx > tiedDests.size() - 1) idx = tiedDests.size() - 1;
            return tiedDests.get(idx).intValue();
        }

        /**
         * Power-of-K-Choices routing (Kotlin lines 6763-6814).
         */
        private int selectSQDestination(int fromNode, int classId,
                                              List<Integer> destinations) {
            int perClassK;
            if (sqDByNodeClass != null
                    && fromNode >= 0 && fromNode < sqDByNodeClass.length
                    && classId >= 0 && classId < sqDByNodeClass[fromNode].length) {
                perClassK = sqDByNodeClass[fromNode][classId];
            } else {
                perClassK = sqD;
            }

            if (destinations.size() <= perClassK) {
                return selectJSQDestination(destinations);
            }

            List<Integer> candidates = new ArrayList<Integer>();
            List<Integer> available = new ArrayList<Integer>(destinations);
            for (int i = 0; i < perClassK; i++) {
                if (available.isEmpty()) break;
                int idx = (int) (routingRng.nextDouble() * available.size());
                if (idx < 0) idx = 0;
                if (idx > available.size() - 1) idx = available.size() - 1;
                candidates.add(available.remove(idx));
            }

            return selectJSQDestination(candidates);
        }

        // ==================== PS Rate Calculations ====================

        /**
         * Dispatch to appropriate PS rate calculation
         * (Kotlin lines 6828-6846).
         */
        private double[] calculatePSRates(int queueIdx, SchedStrategy strategy,
                                          List<PSCustomer> jobs, double c) {
            int n = jobs.size();
            if (n == 0) return new double[0];

            if (strategy == SchedStrategy.PS || strategy == SchedStrategy.LPS) {
                return calculatePSBasicRates(queueIdx, jobs, c);
            } else if (strategy == SchedStrategy.DPS) {
                return calculateDPSRates(queueIdx, jobs, c);
            } else if (strategy == SchedStrategy.GPS) {
                return calculateGPSRates(queueIdx, jobs, c);
            } else if (strategy == SchedStrategy.PSPRIO) {
                return calculatePSPRIORates(queueIdx, jobs, c);
            } else if (strategy == SchedStrategy.DPSPRIO) {
                return calculateDPSPRIORates(queueIdx, jobs, c);
            } else if (strategy == SchedStrategy.GPSPRIO) {
                return calculateGPSPRIORates(queueIdx, jobs, c);
            } else {
                return calculatePSBasicRates(queueIdx, jobs, c);
            }
        }

        /**
         * PS: equal sharing (Kotlin lines 6855-6869).
         */
        private double[] calculatePSBasicRates(int queueIdx, List<PSCustomer> jobs, double c) {
            int n = jobs.size();
            double[] rates = new double[n];
            double sharePerJob = (n > 0) ? Math.min(1.0, c / n) : 0.0;
            for (int idx = 0; idx < n; idx++) {
                rates[idx] = sharePerJob;
            }
            return rates;
        }

        /**
         * DPS: discriminatory PS - rate proportional to class weight
         * (Kotlin lines 6878-6905).
         */
        private double[] calculateDPSRates(int queueIdx, List<PSCustomer> jobs, double c) {
            int n = jobs.size();
            double[] rates = new double[n];
            if (n == 0) return rates;

            double totalWeightedPop = 0.0;
            for (PSCustomer job : jobs) {
                totalWeightedPop += schedWeights[queueIdx][job.classId];
            }

            if (totalWeightedPop <= 0) {
                return calculatePSBasicRates(queueIdx, jobs, c);
            }

            for (int idx = 0; idx < n; idx++) {
                PSCustomer job = jobs.get(idx);
                double weight = schedWeights[queueIdx][job.classId];
                double share = weight * c / totalWeightedPop;
                rates[idx] = Math.min(1.0, share);
            }
            return rates;
        }

        /**
         * GPS: generalized PS (by class) (Kotlin lines 6915-6951).
         */
        private double[] calculateGPSRates(int queueIdx, List<PSCustomer> jobs, double c) {
            int n = jobs.size();
            double[] rates = new double[n];
            if (n == 0) return rates;

            int[] jobsPerClass = new int[numClasses];
            for (PSCustomer job : jobs) {
                jobsPerClass[job.classId]++;
            }

            double totalWeight = 0.0;
            for (int k = 0; k < numClasses; k++) {
                if (jobsPerClass[k] > 0) {
                    totalWeight += schedWeights[queueIdx][k];
                }
            }

            if (totalWeight <= 0) {
                return calculatePSBasicRates(queueIdx, jobs, c);
            }

            for (int idx = 0; idx < n; idx++) {
                PSCustomer job = jobs.get(idx);
                int classId = job.classId;
                double weight = schedWeights[queueIdx][classId];
                int jobsInClass = jobsPerClass[classId];
                double share = (weight / totalWeight / jobsInClass) * c;
                rates[idx] = Math.min(1.0, share);
            }
            return rates;
        }

        /**
         * PSPRIO: PS with strict priorities (Kotlin lines 6962-7009).
         */
        private double[] calculatePSPRIORates(int queueIdx, List<PSCustomer> jobs, double c) {
            int n = jobs.size();
            double[] rates = new double[n];
            if (n == 0) return rates;

            if (n <= c) {
                for (int idx = 0; idx < n; idx++) rates[idx] = 1.0;
                return rates;
            }

            java.util.TreeSet<Integer> priorities = new java.util.TreeSet<Integer>();
            for (PSCustomer job : jobs) {
                priorities.add(Integer.valueOf(job.priority));
            }

            double remainingCapacity = c;
            for (Integer prioBoxed : priorities) {
                int prio = prioBoxed.intValue();
                if (remainingCapacity <= 0) break;

                int count = 0;
                for (PSCustomer job : jobs) {
                    if (job.priority == prio) count++;
                }

                double allocated = Math.min(remainingCapacity, (double) count);
                double sharePerJob = allocated / count;

                for (int idx = 0; idx < jobs.size(); idx++) {
                    PSCustomer job = jobs.get(idx);
                    if (job.priority == prio) {
                        rates[idx] = sharePerJob;
                    }
                }

                remainingCapacity -= allocated;
            }
            return rates;
        }

        /**
         * DPSPRIO: DPS with strict priorities (Kotlin lines 7020-7079).
         */
        private double[] calculateDPSPRIORates(int queueIdx, List<PSCustomer> jobs, double c) {
            int n = jobs.size();
            double[] rates = new double[n];
            if (n == 0) return rates;

            if (n <= c) {
                for (int idx = 0; idx < n; idx++) rates[idx] = 1.0;
                return rates;
            }

            java.util.TreeSet<Integer> priorities = new java.util.TreeSet<Integer>();
            for (PSCustomer job : jobs) {
                priorities.add(Integer.valueOf(job.priority));
            }

            double remainingCapacity = c;
            for (Integer prioBoxed : priorities) {
                int prio = prioBoxed.intValue();
                if (remainingCapacity <= 0) break;

                int count = 0;
                double totalWeightedPop = 0.0;
                for (PSCustomer job : jobs) {
                    if (job.priority == prio) {
                        count++;
                        totalWeightedPop += schedWeights[queueIdx][job.classId];
                    }
                }

                double allocated = Math.min(remainingCapacity, (double) count);

                if (totalWeightedPop <= 0) {
                    double sharePerJob = allocated / count;
                    for (int idx = 0; idx < jobs.size(); idx++) {
                        PSCustomer job = jobs.get(idx);
                        if (job.priority == prio) {
                            rates[idx] = sharePerJob;
                        }
                    }
                } else {
                    for (int idx = 0; idx < jobs.size(); idx++) {
                        PSCustomer job = jobs.get(idx);
                        if (job.priority == prio) {
                            double weight = schedWeights[queueIdx][job.classId];
                            double share = weight * allocated / totalWeightedPop;
                            rates[idx] = Math.min(1.0, share);
                        }
                    }
                }

                remainingCapacity -= allocated;
            }
            return rates;
        }

        /**
         * GPSPRIO: GPS with strict priorities (Kotlin lines 7090-7160).
         */
        private double[] calculateGPSPRIORates(int queueIdx, List<PSCustomer> jobs, double c) {
            int n = jobs.size();
            double[] rates = new double[n];
            if (n == 0) return rates;

            if (n <= c) {
                for (int idx = 0; idx < n; idx++) rates[idx] = 1.0;
                return rates;
            }

            java.util.TreeSet<Integer> priorities = new java.util.TreeSet<Integer>();
            for (PSCustomer job : jobs) {
                priorities.add(Integer.valueOf(job.priority));
            }

            double remainingCapacity = c;
            for (Integer prioBoxed : priorities) {
                int prio = prioBoxed.intValue();
                if (remainingCapacity <= 0) break;

                int count = 0;
                int[] jobsPerClass = new int[numClasses];
                for (PSCustomer job : jobs) {
                    if (job.priority == prio) {
                        count++;
                        jobsPerClass[job.classId]++;
                    }
                }

                double allocated = Math.min(remainingCapacity, (double) count);

                double totalWeight = 0.0;
                for (int k = 0; k < numClasses; k++) {
                    if (jobsPerClass[k] > 0) {
                        totalWeight += schedWeights[queueIdx][k];
                    }
                }

                if (totalWeight <= 0) {
                    double sharePerJob = allocated / count;
                    for (int idx = 0; idx < jobs.size(); idx++) {
                        PSCustomer job = jobs.get(idx);
                        if (job.priority == prio) {
                            rates[idx] = sharePerJob;
                        }
                    }
                } else {
                    for (int idx = 0; idx < jobs.size(); idx++) {
                        PSCustomer job = jobs.get(idx);
                        if (job.priority == prio) {
                            int classId = job.classId;
                            double weight = schedWeights[queueIdx][classId];
                            int jobsInClass = jobsPerClass[classId];
                            double share = (weight / totalWeight / jobsInClass) * allocated;
                            rates[idx] = Math.min(1.0, share);
                        }
                    }
                }

                remainingCapacity -= allocated;
            }
            return rates;
        }

        // ==================== PS Core Functions ====================

        /**
         * Update remaining service work for all PS jobs based on elapsed time
         * (Kotlin lines 7171-7201).
         */
        private void updatePSRemainingWork(int queueIdx, double currentTime) {
            double elapsedTime = currentTime - psLastUpdateTime[queueIdx];
            if (elapsedTime <= 0) {
                psLastUpdateTime[queueIdx] = currentTime;
                return;
            }

            List<PSCustomer> jobs = psJobsInService[queueIdx];
            if (jobs.isEmpty()) {
                psLastUpdateTime[queueIdx] = currentTime;
                return;
            }

            SchedStrategy strategy = schedStrategies[queueIdx];
            double c = getEffectivePSServerCount(queueIdx);

            double[] rates = calculatePSRates(queueIdx, strategy, jobs, c);

            for (int idx = 0; idx < jobs.size(); idx++) {
                PSCustomer job = jobs.get(idx);
                double effectiveRate = rates[idx];
                if (effectiveRate > 0) {
                    // Operational-time work delivered = share * int mu(u) du over
                    // the elapsed wall interval. The share is piecewise constant
                    // between reschedules and mu piecewise constant between
                    // breakpoints, so integrating the class rate over the whole
                    // interval is exact. serviceWorkBetween is the elapsed wall
                    // time for a station without a time-varying rate.
                    double workDone = effectiveRate * serviceWorkBetween(
                            queueIdx, job.classId, psLastUpdateTime[queueIdx], currentTime);
                    job.remainingServiceWork =
                            Math.max(0.0, job.remainingServiceWork - workDone);
                }
            }

            psLastUpdateTime[queueIdx] = currentTime;
        }

        /**
         * Cancel and reschedule all PS departure events (Kotlin lines 7207-7241).
         */
        private void rescheduleAllPSDepartures(int queueIdx) {
            List<PSCustomer> jobs = psJobsInService[queueIdx];
            if (jobs.isEmpty()) return;

            for (PSCustomer job : jobs) {
                if (job.scheduledDepartureEvent != null) {
                    job.scheduledDepartureEvent.cancel();
                    job.scheduledDepartureEvent = null;
                }
            }

            SchedStrategy strategy = schedStrategies[queueIdx];
            double c = getEffectivePSServerCount(queueIdx);
            double[] rates = calculatePSRates(queueIdx, strategy, jobs, c);

            for (int idx = 0; idx < jobs.size(); idx++) {
                PSCustomer job = jobs.get(idx);
                double rate = rates[idx];

                if (job.remainingServiceWork <= 1e-12) {
                    PSDeparture departureEvent = new PSDeparture(queueIdx, job);
                    departureEvent.schedule(1e-12);
                    job.scheduledDepartureEvent = departureEvent;
                } else if (rate > 0) {
                    // Completion when share * int mu(u) du reaches the residual
                    // work: int mu = remaining/rate, whose wall solution is
                    // serviceWallDelay. Equals remaining/rate for a constant rate.
                    double timeToComplete = serviceWallDelay(
                            queueIdx, job.classId, job.remainingServiceWork / rate);
                    PSDeparture departureEvent = new PSDeparture(queueIdx, job);
                    departureEvent.schedule(timeToComplete);
                    job.scheduledDepartureEvent = departureEvent;
                }
            }
        }

        // ==================== PS Arrival Handler ====================

        /**
         * Customer arrival at a PS-family queue (Kotlin lines 7255-7331).
         * Two-arg overload: forkedJob defaults to null.
         */
        private boolean arriveAtPSQueue(int queueIdx, Customer customer) {
            return arriveAtPSQueue(queueIdx, customer, null);
        }

        /**
         * Customer arrival at a PS-family queue (Kotlin lines 7255-7331).
         */
        private boolean arriveAtPSQueue(int queueIdx, Customer customer, ForkedJob forkedJob) {
            int classId = customer.classId;
            double currentTime = ssjSim.time();

            // Update remaining work for all current jobs based on elapsed time
            updatePSRemainingWork(queueIdx, currentTime);

            // Update queue length statistics
            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]++;

            // Track total customers in service
            customersInService[queueIdx]++;

            // LPS: if at capacity, job waits in FCFS queue
            int lpsLimit = lpsLimits[queueIdx];
            if (lpsLimit > 0 && psJobsInService[queueIdx].size() >= lpsLimit) {
                customer.serviceTime = generateServiceTime(queueIdx, classId);
                customer.forkedJob = forkedJob;
                waitQueues[queueIdx].add(customer);

                int stationIdx = serviceStations.get(queueIdx);
                if (stationIdx < fcRegionIndices.size()
                        && fcRegionIndices.get(stationIdx) >= 0) {
                    int regionIdx = fcRegionIndices.get(stationIdx);
                    updateRegionTimeWeightedStats(regionIdx);
                    regionJobEnter(regionIdx, classId);
                    updateRegionArrivalTracking(regionIdx, classId);
                }

                logEvent("LPS_WAIT", serviceStations.get(queueIdx), classId,
                        currentQueueLength[queueIdx][classId],
                        currentBusyServers[queueIdx][classId]);
                return true;
            }

            // Generate service requirement
            double serviceRequirement = generateServiceTime(queueIdx, classId);

            // Create PS customer
            PSCustomer psCustomer = new PSCustomer(
                    classId,
                    customer.priority,
                    customer.systemArrivalTime,
                    customer.queueArrivalTime,
                    serviceRequirement,
                    serviceRequirement,
                    null,
                    forkedJob,
                    -1);

            // Update busy stats BEFORE adding the new job
            updatePSBusyStats(queueIdx);

            psJobsInService[queueIdx].add(psCustomer);

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size()
                    && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobEnter(regionIdx, classId);
                updateRegionArrivalTracking(regionIdx, classId);
            }

            rescheduleAllPSDepartures(queueIdx);

            logEvent("PS_ARRIVAL", serviceStations.get(queueIdx), classId,
                    currentQueueLength[queueIdx][classId],
                    currentBusyServers[queueIdx][classId]);

            return true;
        }

        // ==================== PS Departure Event ====================

        /**
         * Departure event for PS-family queues (Kotlin lines 7339-7481).
         */
        private final class PSDeparture extends SimEvent {
            private final int queueIdx;
            private final PSCustomer customer;

            PSDeparture(int queueIdx, PSCustomer customer) {
                this.queueIdx = queueIdx;
                this.customer = customer;
            }

            @Override
            public void actions() {
                trackEvent();
                double currentTime = ssjSim.time();
                int classId = customer.classId;

                // Update remaining work
                updatePSRemainingWork(queueIdx, currentTime);

                // Update busy stats BEFORE removing the job
                updatePSBusyStats(queueIdx);

                // Remove departing customer
                psJobsInService[queueIdx].remove(customer);

                // Record queue response time
                double queueResponseTime = currentTime - customer.queueArrivalTime;
                responseTimeTally[queueIdx][classId].add(queueResponseTime);
                responseTimeSamples[queueIdx][classId].add(Double.valueOf(queueResponseTime));
                completedCustomers[queueIdx][classId]++;

                checkEventCountStop();

                updateQueueStats(queueIdx, classId);
                currentQueueLength[queueIdx][classId]--;

                int stationIdx = serviceStations.get(queueIdx);
                if (stationIdx < fcRegionIndices.size()
                        && fcRegionIndices.get(stationIdx) >= 0) {
                    int regionIdx = fcRegionIndices.get(stationIdx);
                    updateRegionTimeWeightedStats(regionIdx);
                    regionJobLeave(regionIdx, classId);
                    // completion count and blocked-FIFO release are deferred to
                    // regionExitCompleted once the routing destination is known:
                    // intra-region hops must not release into the transient slot
                }

                // LQN phase-2: the spawn continuation takes over the freed slot
                maybeSpawnOnCompletion(queueIdx, classId, customer.forkedJob);

                customersInService[queueIdx]--;

                logEvent("PS_DEPARTURE", stationIdx, classId,
                        currentQueueLength[queueIdx][classId],
                        currentBusyServers[queueIdx][classId]);

                // LPS: admit next waiting customer
                int lpsLimit = lpsLimits[queueIdx];
                if (lpsLimit > 0 && !waitQueues[queueIdx].isEmpty()
                        && psJobsInService[queueIdx].size() < lpsLimit) {
                    Customer waitingCustomer = waitQueues[queueIdx].poll();
                    int waitingClassId = waitingCustomer.classId;
                    double waitingServiceReq = (waitingCustomer.serviceTime > 0)
                            ? waitingCustomer.serviceTime
                            : generateServiceTime(queueIdx, waitingClassId);
                    PSCustomer admittedPSCustomer = new PSCustomer(
                            waitingClassId,
                            waitingCustomer.priority,
                            waitingCustomer.systemArrivalTime,
                            waitingCustomer.queueArrivalTime,
                            waitingServiceReq,
                            waitingServiceReq,
                            null,
                            waitingCustomer.forkedJob,
                            -1);
                    updatePSBusyStats(queueIdx);
                    psJobsInService[queueIdx].add(admittedPSCustomer);
                }

                rescheduleAllPSDepartures(queueIdx);

                // Route customer to next destination
                int currentNode = serviceNodes.get(queueIdx);
                RoutingResult routingResult = selectDestination(currentNode, classId);
                int destNode = routingResult.destNode;
                int destClassId = routingResult.destClassId;
                // FCR: on a true region exit, count the completion and release
                // blocked customers; a no-op for intra-region hops
                regionExitCompleted(queueIdx, classId, destNode, destClassId);

                if (destNode >= 0 && sinkNodes.contains(Integer.valueOf(destNode))) {
                    double systemResponseTime = currentTime - customer.systemArrivalTime;
                    systemResponseTimeTally[classId].add(systemResponseTime);
                    systemCompletedCustomers[classId]++;
                } else if (destNode >= 0 && forkNodes.contains(Integer.valueOf(destNode))) {
                    long parentJobId = nextJobId++;
                    handleForkArrival(destNode, parentJobId, destClassId,
                            customer.systemArrivalTime);
                } else if (destNode >= 0 && joinNodes.contains(Integer.valueOf(destNode))) {
                    ForkedJob fj = customer.forkedJob;
                    if (fj != null) {
                        ForkedJob updatedForkedJob = new ForkedJob(
                                fj.forkJobId,
                                fj.parentJobId,
                                destClassId,
                                fj.priority,
                                fj.systemArrivalTime,
                                currentTime,
                                fj.randomRank);
                        handleJoinArrival(destNode, updatedForkedJob);
                    } else {
                        handleUnknownJoinArrival(destNode, destClassId,
                                customer.systemArrivalTime);
                    }
                } else if (destNode >= 0 && placeNodes.contains(Integer.valueOf(destNode))) {
                    handlePlaceArrival(destNode, destClassId, customer.systemArrivalTime);
                } else if (destNode >= 0 && transitionNodes.contains(Integer.valueOf(destNode))) {
                    checkAndFireTransitions();
                } else if (destNode >= 0) {
                    int nextQueueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
                    if (nextQueueIdx >= 0) {
                        ForkedJob fj = customer.forkedJob;
                        Customer nextCustomer = new Customer(
                                destClassId, classPrio[destClassId],
                                customer.systemArrivalTime, currentTime,
                                siroRng.nextDouble());
                        if (fj != null) {
                            ForkedJob nextForkedJob = new ForkedJob(
                                    fj.forkJobId,
                                    fj.parentJobId,
                                    destClassId,
                                    fj.priority,
                                    fj.systemArrivalTime,
                                    currentTime,
                                    fj.randomRank);
                            arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob);
                        } else {
                            arriveAtQueue(nextQueueIdx, nextCustomer);
                        }
                    }
                }
            }
        }

        // ==================== PAS (pass-and-swap) handlers ====================

        /**
         * Result of a pass-and-swap transition: the departing customer and the
         * surviving customers in their new order (oldest first).
         */
        private static final class PasSwapResult {
            final Customer departing;
            final List<Customer> survivors;

            PasSwapResult(Customer departing, List<Customer> survivors) {
                this.departing = departing;
                this.survivors = survivors;
            }
        }

        /**
         * Evaluates the total service rate mu(c) over the first prefixLen jobs
         * of the ordered list (Dorsman and Gardner 2024). Class indices passed
         * to the rate function are 0-based, matching Customer.classId.
         */
        private double pasMu(int queueIdx, List<Customer> jobs, int prefixLen) {
            if (prefixLen <= 0) {
                return 0.0;
            }
            jline.util.matrix.Matrix prefix = new jline.util.matrix.Matrix(1, prefixLen);
            for (int i = 0; i < prefixLen; i++) {
                prefix.set(0, i, jobs.get(i).classId);
            }
            return pasSvcRateFun[queueIdx].apply(prefix);
        }

        /**
         * Applies the pass-and-swap mechanism (Dorsman and Gardner 2024,
         * Sect. 2.3) triggered by the service completion at position p of the
         * ordered job list. The completing job scans towards newer positions
         * for the first swappable class, takes its place and ejects it; the
         * ejected job repeats until one with no swappable successor departs.
         * Class indices (Customer.classId and the swap graph G) are 0-based.
         */
        private PasSwapResult passAndSwap(List<Customer> jobs, int p,
                                          jline.util.matrix.Matrix G) {
            int n = jobs.size();
            List<Integer> chain = new ArrayList<Integer>();
            chain.add(Integer.valueOf(p));
            int movingClass = jobs.get(p).classId;
            int cur = p;
            while (true) {
                int q = -1;
                for (int j = cur + 1; j < n; j++) {
                    boolean swappable = (G == null)
                            || (G.get(movingClass, jobs.get(j).classId) != 0);
                    if (swappable) {
                        q = j;
                        break;
                    }
                }
                if (q < 0) {
                    break;
                }
                chain.add(Integer.valueOf(q));
                movingClass = jobs.get(q).classId;
                cur = q;
            }
            int hole = chain.get(0).intValue();
            int depPos = chain.get(chain.size() - 1).intValue();
            Customer departing = jobs.get(depPos);
            // Each job in the chain moves forward into the next chain position.
            Customer[] arr = jobs.toArray(new Customer[0]);
            for (int i = 0; i < chain.size() - 1; i++) {
                arr[chain.get(i + 1).intValue()] = jobs.get(chain.get(i).intValue());
            }
            List<Customer> survivors = new ArrayList<Customer>();
            for (int i = 0; i < n; i++) {
                if (i == hole) {
                    continue;
                }
                survivors.add(arr[i]);
            }
            return new PasSwapResult(departing, survivors);
        }

        /**
         * Accumulates per-class busy time for a PAS station. A position p is
         * "in service" when its marginal rate delta-mu = mu(c1..cp) -
         * mu(c1..c_{p-1}) is positive; utilization is then E[in-service]/servers,
         * matching the CTMC metric.
         */
        private void updatePASBusyStats(int queueIdx) {
            double currentTime = ssjSim.time();
            double elapsed = currentTime - pasLastBusyUpdateTime[queueIdx];
            if (elapsed <= 0) {
                pasLastBusyUpdateTime[queueIdx] = currentTime;
                return;
            }
            List<Customer> jobs = pasList[queueIdx];
            int n = jobs.size();
            if (n > 0) {
                double muPrev = 0.0;
                for (int p = 0; p < n; p++) {
                    double muCur = pasMu(queueIdx, jobs, p + 1);
                    double dmu = muCur - muPrev;
                    muPrev = muCur;
                    if (dmu > 0) {
                        totalBusyTime[queueIdx][jobs.get(p).classId] += elapsed;
                    }
                }
            }
            pasLastBusyUpdateTime[queueIdx] = currentTime;
        }

        /**
         * Cancels and reschedules the aggregate completion event of a PAS
         * station. The whole station completes at total rate mu(c) (exponential,
         * hence memoryless): on any state change the pending event is replaced
         * by a fresh Exp(mu(c)) draw.
         */
        private void reschedulePASDeparture(int queueIdx) {
            if (pasDepartureEvent[queueIdx] != null) {
                pasDepartureEvent[queueIdx].cancel();
                pasDepartureEvent[queueIdx] = null;
            }
            List<Customer> jobs = pasList[queueIdx];
            int n = jobs.size();
            if (n == 0) {
                return;
            }
            double totalRate = pasMu(queueIdx, jobs, n);
            if (totalRate <= 0) {
                return;
            }
            double t = -Math.log(pasRng.nextDouble()) / totalRate;
            PASDeparture ev = new PASDeparture(queueIdx);
            ev.schedule(t);
            pasDepartureEvent[queueIdx] = ev;
        }

        /**
         * Customer arrival at a PAS (pass-and-swap) queue. The job is appended
         * to the back of the ordered list and the station completion event is
         * rescheduled. Capacity has already been checked in arriveAtQueue.
         */
        private boolean arriveAtPASQueue(int queueIdx, Customer customer) {
            int classId = customer.classId;

            // Flush busy time accrued with the pre-arrival composition.
            updatePASBusyStats(queueIdx);

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]++;
            customersInService[queueIdx]++;

            pasList[queueIdx].add(customer);

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size()
                    && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobEnter(regionIdx, classId);
                updateRegionArrivalTracking(regionIdx, classId);
            }

            reschedulePASDeparture(queueIdx);

            logEvent("PAS_ARRIVAL", stationIdx, classId,
                    currentQueueLength[queueIdx][classId],
                    currentBusyServers[queueIdx][classId]);
            return true;
        }

        /** Aggregate service-completion event for a PAS (pass-and-swap) station. */
        private final class PASDeparture extends SimEvent {
            private final int queueIdx;

            PASDeparture(int queueIdx) {
                this.queueIdx = queueIdx;
            }

            @Override
            public void actions() {
                trackEvent();
                double currentTime = ssjSim.time();
                List<Customer> jobs = pasList[queueIdx];
                int n = jobs.size();
                if (n == 0) {
                    pasDepartureEvent[queueIdx] = null;
                    return;
                }

                // Flush busy time accrued with the pre-departure composition.
                updatePASBusyStats(queueIdx);

                // Choose the completing position p proportionally to delta-mu.
                double[] cum = new double[n];
                double muPrev = 0.0;
                double acc = 0.0;
                for (int p = 0; p < n; p++) {
                    double muCur = pasMu(queueIdx, jobs, p + 1);
                    double dmu = muCur - muPrev;
                    muPrev = muCur;
                    if (dmu < 0) {
                        dmu = 0.0;
                    }
                    acc += dmu;
                    cum[p] = acc;
                }
                double total = acc;
                double u = pasRng.nextDouble() * total;
                int pos = n - 1;
                for (int p = 0; p < n; p++) {
                    if (u <= cum[p]) {
                        pos = p;
                        break;
                    }
                }

                PasSwapResult res = passAndSwap(jobs, pos, pasSwapGraph[queueIdx]);
                Customer departing = res.departing;
                int classId = departing.classId;

                pasList[queueIdx] = res.survivors;
                pasDepartureEvent[queueIdx] = null;

                double queueResponseTime = currentTime - departing.queueArrivalTime;
                responseTimeTally[queueIdx][classId].add(queueResponseTime);
                responseTimeSamples[queueIdx][classId].add(Double.valueOf(queueResponseTime));
                completedCustomers[queueIdx][classId]++;

                checkEventCountStop();

                updateQueueStats(queueIdx, classId);
                currentQueueLength[queueIdx][classId]--;

                int stationIdx = serviceStations.get(queueIdx);
                if (stationIdx < fcRegionIndices.size()
                        && fcRegionIndices.get(stationIdx) >= 0) {
                    int regionIdx = fcRegionIndices.get(stationIdx);
                    updateRegionTimeWeightedStats(regionIdx);
                    regionJobLeave(regionIdx, classId);
                    // completion count and blocked-FIFO release are deferred to
                    // regionExitCompleted once the routing destination is known:
                    // intra-region hops must not release into the transient slot
                }

                // LQN phase-2: the spawn continuation takes over the freed slot
                maybeSpawnOnCompletion(queueIdx, classId, departing);

                customersInService[queueIdx]--;

                logEvent("PAS_DEPARTURE", stationIdx, classId,
                        currentQueueLength[queueIdx][classId],
                        currentBusyServers[queueIdx][classId]);

                // Reschedule the station completion with the new composition.
                reschedulePASDeparture(queueIdx);

                // Route the departing customer onward.
                routePASDeparture(queueIdx, departing, currentTime);
            }
        }

        /** Routes a departing PAS customer to its next destination. */
        private void routePASDeparture(int queueIdx, Customer customer, double currentTime) {
            int classId = customer.classId;
            int currentNode = serviceNodes.get(queueIdx);
            RoutingResult routingResult = selectDestination(currentNode, classId);
            int destNode = routingResult.destNode;
            int destClassId = routingResult.destClassId;
            // FCR: on a true region exit, count the completion and release
            // blocked customers; a no-op for intra-region hops
            regionExitCompleted(queueIdx, classId, destNode, destClassId);

            if (destNode >= 0 && sinkNodes.contains(Integer.valueOf(destNode))) {
                double systemResponseTime = currentTime - customer.systemArrivalTime;
                systemResponseTimeTally[classId].add(systemResponseTime);
                systemCompletedCustomers[classId]++;
            } else if (destNode >= 0 && forkNodes.contains(Integer.valueOf(destNode))) {
                long parentJobId = nextJobId++;
                handleForkArrival(destNode, parentJobId, destClassId,
                        customer.systemArrivalTime);
            } else if (destNode >= 0 && joinNodes.contains(Integer.valueOf(destNode))) {
                ForkedJob fj = customer.forkedJob;
                if (fj != null) {
                    ForkedJob updatedForkedJob = new ForkedJob(
                            fj.forkJobId,
                            fj.parentJobId,
                            destClassId,
                            fj.priority,
                            fj.systemArrivalTime,
                            currentTime,
                            fj.randomRank);
                    handleJoinArrival(destNode, updatedForkedJob);
                } else {
                    handleUnknownJoinArrival(destNode, destClassId,
                            customer.systemArrivalTime);
                }
            } else if (destNode >= 0 && placeNodes.contains(Integer.valueOf(destNode))) {
                handlePlaceArrival(destNode, destClassId, customer.systemArrivalTime);
            } else if (destNode >= 0 && transitionNodes.contains(Integer.valueOf(destNode))) {
                checkAndFireTransitions();
            } else if (destNode >= 0) {
                int nextQueueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
                if (nextQueueIdx >= 0) {
                    ForkedJob fj = customer.forkedJob;
                    Customer nextCustomer = new Customer(
                            destClassId, classPrio[destClassId],
                            customer.systemArrivalTime, currentTime,
                            siroRng.nextDouble());
                    if (fj != null) {
                        ForkedJob nextForkedJob = new ForkedJob(
                                fj.forkJobId,
                                fj.parentJobId,
                                destClassId,
                                fj.priority,
                                fj.systemArrivalTime,
                                currentTime,
                                fj.randomRank);
                        arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob);
                    } else {
                        arriveAtQueue(nextQueueIdx, nextCustomer);
                    }
                }
            }
        }

        // ==================== Forward-reference stubs (PART 6+) ====================
        //
        // The following methods are called from the routing/PS code above
        // but are translated in subsequent chunks (Kotlin lines 7501+).
        // Stubs throw {@link UnsupportedOperationException} so that any
        // accidental invocation prior to translation is detected loudly.

        // selectDestination / handleForkArrival / handleJoinArrival /
        // handleUnknownJoinArrival / arriveAtQueueForked translated in PART 8
        // (Kotlin lines 11437-11935); the stubs are now redirector methods that
        // call the implementations defined later in the file.

        /** Forward stub — Kotlin line 12672. */
        /**
         * Handle a token arriving at a Place node.
         * Adds token to place and triggers check for enabled transitions.
         * Kotlin line 12672.
         */
        private void handlePlaceArrival(int placeNodeIdx, int classId,
                                        double systemArrivalTime) {
            int placeListIdx = placeNodes.indexOf(placeNodeIdx);
            if (placeListIdx < 0) return;

            // Update time-weighted token count before adding new token
            double currentTime = ssjSim.time();
            double elapsed = currentTime - lastPlaceUpdateTime[placeListIdx];
            if (elapsed > 0) {
                for (int k = 0; k < numClasses; k++) {
                    totalPlaceTokenTime[placeListIdx][k] +=
                            placeTokens[placeListIdx][k] * elapsed;
                }
                lastPlaceUpdateTime[placeListIdx] = currentTime;
            }

            // Add token to place (marking). For an ordinary place the token is immediately
            // available to output transitions; for a queueing place it enters the embedded
            // queue and becomes available only after service completion (depository).
            placeTokens[placeListIdx][classId]++;
            if (isQueueingPlace[placeListIdx]) {
                placeWaiting[placeListIdx].addLast(Integer.valueOf(classId));
                tryStartPlaceService(placeListIdx);
            }

            // Check and fire enabled transitions
            checkAndFireTransitions();
        }

        /** Allocates a generic array of ArrayDeque without unchecked-cast warnings at call sites. */
        @SuppressWarnings("unchecked")
        private java.util.ArrayDeque<Integer>[] newDequeArray(int n) {
            return (java.util.ArrayDeque<Integer>[]) new java.util.ArrayDeque<?>[n];
        }

        /** Tokens of a place currently available to output transitions: the depository for a
         *  queueing place, or the whole marking for an ordinary place. */
        private int placeAvail(int placeListIdx, int classId) {
            return isQueueingPlace[placeListIdx]
                    ? placeDepository[placeListIdx][classId]
                    : placeTokens[placeListIdx][classId];
        }

        /** Consumes {@code n} class-{@code classId} tokens from place {@code placeListIdx} when a
         *  transition fires: removes them from the marking and, for a queueing place, from the
         *  depository. */
        private void consumePlaceTokens(int placeListIdx, int classId, int n) {
            placeTokens[placeListIdx][classId] -= n;
            if (isQueueingPlace[placeListIdx]) {
                placeDepository[placeListIdx][classId] -= n;
            }
        }

        /** Flushes the busy-server time integral of a queueing place up to the current
         *  clock, before its in-service count changes. Enables per-class utilization. */
        private void accruePlaceBusy(int placeListIdx) {
            double now = ssjSim.time();
            double dt = now - lastPlaceBusyUpdateTime[placeListIdx];
            if (dt > 0) {
                int[] inSvc = currentPlaceInService[placeListIdx];
                double[] busyTime = totalPlaceBusyTime[placeListIdx];
                for (int k = 0; k < inSvc.length; k++) {
                    if (inSvc[k] > 0) {
                        busyTime[k] += inSvc[k] * dt;
                    }
                }
            }
            lastPlaceBusyUpdateTime[placeListIdx] = now;
        }

        /** Starts service for waiting tokens of a queueing place while a server is free. */
        private void tryStartPlaceService(int placeListIdx) {
            java.util.ArrayDeque<Integer> waiting = placeWaiting[placeListIdx];
            while (placeBusy[placeListIdx] < placeNumServers[placeListIdx] && !waiting.isEmpty()) {
                int classId = waiting.pollFirst().intValue();
                accruePlaceBusy(placeListIdx);
                placeBusy[placeListIdx]++;
                currentPlaceInService[placeListIdx][classId]++;
                double serviceTime = generatePlaceServiceTime(placeListIdx, classId);
                new PlaceDeparture(placeListIdx, classId).schedule(serviceTime);
            }
        }

        /** Service-completion event for a token in a queueing place's embedded queue. */
        private final class PlaceDeparture extends SimEvent {
            final int placeListIdx;
            final int classId;

            PlaceDeparture(int placeListIdx, int classId) {
                this.placeListIdx = placeListIdx;
                this.classId = classId;
            }

            @Override
            public void actions() {
                placeDepartureActions(placeListIdx, classId);
            }
        }

        /** Handle a service completion at a queueing place: move the token to the depository
         *  (now available to output transitions), start the next waiting token's service, and
         *  re-evaluate transition enabling. The place marking is unchanged by a completion. */
        private void placeDepartureActions(int placeListIdx, int classId) {
            accruePlaceBusy(placeListIdx);
            placeBusy[placeListIdx]--;
            currentPlaceInService[placeListIdx][classId]--;
            placeDepository[placeListIdx][classId]++;
            tryStartPlaceService(placeListIdx);
            checkAndFireTransitions();
        }

        /** Generates an embedded-queue service time for a queueing place/class. Supports EXP,
         *  DET, IMMEDIATE, and renewal phase-type families (PH/APH/HyperExp/Coxian/Cox2/ME);
         *  correlated processes (MAP/MMPP2/RAP/MMAP) are rejected loudly. */
        private double generatePlaceServiceTime(int placeListIdx, int classId) {
            ProcessType procType = placeSvcType[placeListIdx][classId];
            if (procType == ProcessType.DISABLED || procType == null) {
                String placeName = sn.nodes.get(placeNodes.get(placeListIdx)).getName();
                String className = sn.jobclasses.get(classId).getName();
                throw new RuntimeException(
                        "LDES: token of class '" + className + "' entered queueing place '"
                        + placeName + "' but no service is defined for this class. "
                        + "Call setService for every token color of a queueing place.");
            }
            RandomVariateGen gen = placeSvcGen[placeListIdx][classId];
            if (gen != null) {
                return gen.nextDouble();
            }
            if (procType == ProcessType.ME) {
                MatrixCell proc = placeSvcProc[placeListIdx][classId];
                java.util.Random rng = placeSvcRng[placeListIdx][classId];
                if (proc != null && rng != null && proc.size() >= 2) {
                    Me_sample.MeSampler s = placeSvcMeSampler[placeListIdx][classId];
                    if (s == null) {
                        s = new Me_sample.MeSampler(meAlphaOf(proc, "service"), proc.get(0));
                        placeSvcMeSampler[placeListIdx][classId] = s;
                    }
                    return s.next(rng);
                }
            }
            if (procType == ProcessType.PH || procType == ProcessType.APH
                    || procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN
                    || procType == ProcessType.COX2) {
                MatrixCell proc = placeSvcProc[placeListIdx][classId];
                java.util.Random rng = placeSvcRng[placeListIdx][classId];
                if (proc != null && rng != null && proc.size() >= 2) {
                    double[] samples = jline.api.mam.Map_sample.map_sample(proc.get(0), proc.get(1), 1L, rng);
                    return samples[0];
                }
            }
            String placeName = sn.nodes.get(placeNodes.get(placeListIdx)).getName();
            throw new RuntimeException(
                    "LDES: queueing place '" + placeName + "' uses service process " + procType
                    + ", which the embedded-queue algorithm does not yet support "
                    + "(supported: EXP, DET, IMMEDIATE, PH/APH/HyperExp/Coxian/Cox2/ME).");
        }

        /** Generates a service time for queue/class. Kotlin line 12195. */
        private double generateServiceTime(int queueIdx, int classId) {
            ProcessType procType = serviceProcessType[queueIdx][classId];

            // Check for disabled service
            if (procType == ProcessType.DISABLED) {
                String stationName = sn.stations.get(serviceStations.get(queueIdx)).getName();
                String className = sn.jobclasses.get(classId).getName();
                throw new RuntimeException(
                        "LDES: Job of class '" + className + "' requested service at station '"
                        + stationName + "' but service is DISABLED for this class-station pair. "
                        + "Check model routing.");
            }

            // Generate base service time
            double baseServiceTime = 0.0;

            // Check for Replayer/trace distribution first
            if (procType == ProcessType.REPLAYER) {
                TraceSampler traceSampler = serviceTraceSamplers[queueIdx][classId];
                if (traceSampler != null) {
                    baseServiceTime = traceSampler.nextSample();
                }
            } else {
                RandomVariateGen gen = serviceGens[queueIdx][classId];
                if (gen != null) {
                    baseServiceTime = gen.nextDouble();
                } else if (procType == ProcessType.MMAP) {
                    MatrixCell proc = serviceProc[queueIdx][classId];
                    Random rng = serviceRng[queueIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 2) {
                        jline.io.Ret.mamMMAPSample mmapResult =
                                jline.api.mam.Mmap_sample.mmap_sample(proc, 1L, rng);
                        baseServiceTime = mmapResult.getSamples()[0];
                    }
                } else if (procType == ProcessType.MAP || procType == ProcessType.MMPP2) {
                    // Correlated service process: carry the modulating phase across
                    // services so autocorrelation is preserved (an M/MAP/1 queue length
                    // differs from the renewal M/PH/1 value with the same marginal).
                    MatrixCell proc = serviceProc[queueIdx][classId];
                    Random rng = serviceRng[queueIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 2) {
                        Map_sample.MapSampler s = serviceMapSampler[queueIdx][classId];
                        if (s == null) {
                            s = new Map_sample.MapSampler(proc.get(0), proc.get(1));
                            serviceMapSampler[queueIdx][classId] = s;
                        }
                        baseServiceTime = s.next(rng);
                    }
                } else if (procType == ProcessType.RAP) {
                    // Correlated rational service process: the conditional vector is
                    // propagated across services, preserving the RAP autocorrelation.
                    MatrixCell proc = serviceProc[queueIdx][classId];
                    Random rng = serviceRng[queueIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 2) {
                        Rap_sample.RapSampler s = serviceRapSampler[queueIdx][classId];
                        if (s == null) {
                            s = new Rap_sample.RapSampler(proc.get(0), proc.get(1));
                            serviceRapSampler[queueIdx][classId] = s;
                        }
                        baseServiceTime = s.next(rng);
                    }
                } else if (procType == ProcessType.ME) {
                    // Matrix-exponential renewal service: inverse-CDF sampling from a
                    // cached table (map_sample's CTMC walk is invalid for a genuine ME).
                    MatrixCell proc = serviceProc[queueIdx][classId];
                    Random rng = serviceRng[queueIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 2) {
                        Me_sample.MeSampler s = serviceMeSampler[queueIdx][classId];
                        if (s == null) {
                            s = new Me_sample.MeSampler(meAlphaOf(proc, "service"), proc.get(0));
                            serviceMeSampler[queueIdx][classId] = s;
                        }
                        baseServiceTime = s.next(rng);
                    }
                } else if (procType == ProcessType.PH || procType == ProcessType.APH
                        || procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN
                        || procType == ProcessType.COX2) {
                    // Renewal phase-type: D1 restarts from a fixed entry vector, so
                    // per-service resampling from the marginal is exact.
                    MatrixCell proc = serviceProc[queueIdx][classId];
                    Random rng = serviceRng[queueIdx][classId];
                    if (proc != null && rng != null && proc.size() >= 2) {
                        Matrix D0 = proc.get(0);
                        Matrix D1 = proc.get(1);
                        double[] samples =
                                jline.api.mam.Map_sample.map_sample(D0, D1, 1L, rng);
                        baseServiceTime = samples[0];
                    }
                } else if (procType == ProcessType.NHPP) {
                    // Time-varying service rate mu(t). The service requirement is
                    // returned in operational time (cumulative intensity tau): a
                    // job completes when int mu(u) du delivered to it reaches an
                    // Exp(1) budget, since the NHPP completions are a unit-rate
                    // Poisson process in tau. Callers convert this tau budget to a
                    // wall-clock completion via serviceWallAfterWork, which makes
                    // the result exact whether or not service is interrupted (PS,
                    // preemption) or rescaled (load/class dependence). A station
                    // with a constant rate has tau == wall, so this coincides with
                    // the ordinary duration everywhere else.
                    double[][] sched = serviceSchedule[queueIdx][classId];
                    if (sched != null) {
                        baseServiceTime =
                                -Math.log(1.0 - serviceRng[queueIdx][classId].nextDouble());
                    }
                }
            }

            // Apply load-dependent scaling if this station has load dependence
            // (PS scaling is handled in the PS rate calculation)
            if (isLoadDependent[queueIdx] && !isPSScheduling(schedStrategies[queueIdx])) {
                int totalJobs = getTotalCustomersAtStation(queueIdx);
                if (totalJobs > 0) {
                    double[] scalingArray = lldScaling[queueIdx];
                    if (scalingArray != null) {
                        int scalingIdx = Math.min(totalJobs - 1, scalingArray.length - 1);
                        double scalingFactor = scalingArray[scalingIdx];
                        if (scalingFactor > 0) {
                            baseServiceTime /= scalingFactor;
                        }
                    }
                }
            }

            // Apply class dependence: beta_{i,r}(n) is the service RATE for class r
            // at the current per-class population n, so the service time is divided
            // by it (Sauer 1983, eq. (40)). This subsumes the former class-dependence tables
            // (which were rates, and likewise divided) and the LJD tables (which
            // were time scalings, and multiplied); the inversion is applied when the
            // table is wrapped as beta. A station with no class dependence has a
            // null entry and is left unscaled.
            if (hasCd && cdFunctions != null && cdFunctions[queueIdx] != null) {
                Matrix nvecCd = new Matrix(1, numClasses);
                for (int c = 0; c < numClasses; c++) {
                    nvecCd.set(0, c, (double) currentQueueLength[queueIdx][c]);
                }
                Matrix bval = cdFunctions[queueIdx].apply(nvecCd);
                double scalingFactor = (bval.length() > 1) ? bval.get(classId) : bval.get(0);
                if (scalingFactor > 1e-10) {
                    baseServiceTime = baseServiceTime / scalingFactor;
                }
            }

            baseServiceTime = resolveZeroAtom(baseServiceTime, procType,
                    "service time at station " + queueIdx + " class " + classId);
            baseServiceTime = slotSnap(baseServiceTime, "service time at station "
                    + queueIdx + " class " + classId);

            // Track sample for control variates if enabled and after warmup
            if (useControlVariates && warmupDone && baseServiceTime > 0) {
                serviceSampleSum[queueIdx][classId] += baseServiceTime;
                serviceSampleCount[queueIdx][classId]++;
            }

            return baseServiceTime;
        }

        // ==================== PART 6 TRANSLATION (Kotlin lines 7548-8973) ====================

        /**
         * Customer arrives at a service node (Queue or Delay).
         * Checks both total station capacity (cap) and per-class capacity
         * (classcap).  Returns true if customer was accepted, false if dropped
         * due to capacity exceeded.
         *
         * <p>Kotlin overload with default {@code skipCapacityCheck = false}.
         */
        private boolean arriveAtQueue(int queueIdx, Customer customer) {
            return arriveAtQueue(queueIdx, customer, false);
        }

        /**
         * Customer arrives at a service node (Queue or Delay).
         * Checks both total station capacity (cap) and per-class capacity
         * (classcap).  Returns true if customer was accepted, false if dropped
         * due to capacity exceeded.
         *
         * @param skipCapacityCheck If true, skip the capacity check
         *                          (used for BAS-admitted jobs).
         */
        private boolean arriveAtQueue(int queueIdx, Customer customer, boolean skipCapacityCheck) {
            int classId = customer.classId;

            // Track arrivals for arrival rate calculation (including jobs that
            // will be dropped)
            if (warmupDone) {
                arrivedCustomers[queueIdx][classId]++;
            }

            // Check for balking BEFORE joining queue
            if (hasBalkingConfig[queueIdx][classId] && shouldBalk(queueIdx, classId)) {
                if (warmupDone) {
                    balkedCustomers[queueIdx][classId]++;
                }
                logEvent("BALK", serviceStations.get(queueIdx), classId,
                        currentQueueLength[queueIdx][classId],
                        currentBusyServers[queueIdx][classId]);
                return false;
            }

            if (!skipCapacityCheck) {
                // Check total station capacity
                int currentTotal = getTotalCustomersAtStation(queueIdx);
                if (currentTotal >= bufferCapacities[queueIdx]) {
                    if (hasRetrialConfig[queueIdx][classId]) {
                        int maxAttempts = retrialMaxAttemptsConfig[queueIdx][classId];
                        scheduleRetrial(queueIdx, customer, 1, maxAttempts);
                        logEvent("RETRIAL_CAPACITY", serviceStations.get(queueIdx), classId,
                                currentQueueLength[queueIdx][classId],
                                currentBusyServers[queueIdx][classId]);
                        return false;
                    }
                    if (warmupDone) {
                        droppedCustomers[queueIdx][classId]++;
                    }
                    logEvent("DROP_CAPACITY", serviceStations.get(queueIdx), classId,
                            currentQueueLength[queueIdx][classId],
                            currentBusyServers[queueIdx][classId]);
                    return false;
                }

                // Check per-class capacity constraint (only if explicitly set)
                if (classCapacities[queueIdx][classId] < Integer.MAX_VALUE) {
                    int currentClassCount = currentQueueLength[queueIdx][classId];
                    if (currentClassCount >= classCapacities[queueIdx][classId]) {
                        if (hasRetrialConfig[queueIdx][classId]) {
                            int maxAttempts = retrialMaxAttemptsConfig[queueIdx][classId];
                            scheduleRetrial(queueIdx, customer, 1, maxAttempts);
                            logEvent("RETRIAL_CLASS_CAPACITY", serviceStations.get(queueIdx),
                                    classId, currentQueueLength[queueIdx][classId],
                                    currentBusyServers[queueIdx][classId]);
                            return false;
                        }
                        if (warmupDone) {
                            droppedCustomers[queueIdx][classId]++;
                        }
                        logEvent("DROP_CLASS_CAPACITY", serviceStations.get(queueIdx),
                                classId, currentQueueLength[queueIdx][classId],
                                currentBusyServers[queueIdx][classId]);
                        return false;
                    }
                }
            }

            // Check finite capacity region constraints
            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);

                // Check global region capacity
                int currentRegionTotal = getTotalCustomersInRegion(regionIdx);
                int globalMax = (regionIdx >= 0 && regionIdx < fcRegionGlobalMax.size())
                        ? fcRegionGlobalMax.get(regionIdx) : Integer.MAX_VALUE;
                if (currentRegionTotal >= globalMax) {
                    return rejectAtRegion(regionIdx, classId, queueIdx, stationIdx, customer, "REGION_CAPACITY");
                }

                // Check global region memory budget (weighted by per-class size).
                // JMT semantics: admit iff usedMemory + incomingSize <= maxMemory.
                if (exceedsRegionMemBudget(regionIdx, classId)) {
                    return rejectAtRegion(regionIdx, classId, queueIdx, stationIdx, customer, "REGION_MEMORY");
                }

                // Check per-class region capacity
                int regionClassMax = (int) fcRegionClassMax.get(regionIdx, classId);
                if (regionClassMax < Integer.MAX_VALUE
                        && getClassJobsInRegion(regionIdx, classId) >= regionClassMax) {
                    return rejectAtRegion(regionIdx, classId, queueIdx, stationIdx, customer, "REGION_CLASS_CAPACITY");
                }

                // Check linear admission constraints: A * n_new <= b
                if (fcRegionLinConA[regionIdx].length > 0) {
                    for (int c = 0; c < fcRegionLinConA[regionIdx].length; c++) {
                        double lhs = fcRegionLinConA[regionIdx][c][classId];
                        for (int r = 0; r < numClasses; r++) {
                            lhs += fcRegionLinConA[regionIdx][c][r]
                                    * currentJobsInRegion[regionIdx][r];
                        }
                        if (lhs > fcRegionLinConb[regionIdx][c]) {
                            return rejectAtRegion(regionIdx, classId, queueIdx, stationIdx, customer, "REGION_LINCON");
                        }
                    }
                }
            }

            // Dispatch to PAS handler for pass-and-swap (order-independent) scheduling
            if (isPASStation[queueIdx]) {
                return arriveAtPASQueue(queueIdx, customer);
            }

            // Dispatch to PS handler for Processor Sharing scheduling
            if (isPSScheduling(schedStrategies[queueIdx])) {
                return arriveAtPSQueue(queueIdx, customer);
            }

            // Dispatch to preemptive LCFS handler
            if (isPreemptiveScheduling[queueIdx]) {
                return arriveAtPreemptiveLCFSQueue(queueIdx, customer);
            }

            // Dispatch to polling handler
            if (isPollingStation[queueIdx]) {
                return arriveAtPollingQueue(queueIdx, customer);
            }

            // For SJF/LJF, we must generate service time upon arrival to sort
            SchedStrategy strategy = schedStrategies[queueIdx];
            if (strategy == SchedStrategy.SJF || strategy == SchedStrategy.LJF) {
                customer.serviceTime = generateServiceTime(queueIdx, classId);
            }

            // Update queue length statistics
            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]++;

            // Re-scale in-service jobs at class-dependent stations when queue length changes
            if ((hasCd || hasLld) && !isDelayNode.get(queueIdx).booleanValue()) {
                rescaleStateDepInServiceJobs(queueIdx);
            }

            // Track max queue length reached
            int totalAtStation = getTotalCustomersAtStation(queueIdx);
            if (totalAtStation > maxQueueLengthReached) {
                maxQueueLengthReached = totalAtStation;
            }

            // Update region job counts if in a region
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobEnter(regionIdx, classId);
                updateRegionArrivalTracking(regionIdx, classId);
            }

            if (isBatchServiceStation[queueIdx]) {
                // Bulk (batch Markovian) server: enqueue and, if the station was
                // empty, start the BMSP clock. Completions are handled in batches
                // by BatchServiceDeparture.
                handleBatchServiceArrival(queueIdx, customer);
                logEvent("ARRIVAL", serviceStations.get(queueIdx), classId,
                        currentQueueLength[queueIdx][classId],
                        currentBusyServers[queueIdx][classId]);
                return true;
            }

            if (isDelayNode.get(queueIdx).booleanValue()) {
                // Delay node (infinite server): always start service immediately
                customersInService[queueIdx]++;

                updateBusyStats(queueIdx, classId);
                currentBusyServers[queueIdx][classId]++;
                lastBusyUpdateTime[queueIdx][classId] = ssjSim.time();

                double serviceTime = generateServiceTime(queueIdx, classId);
                Event departureEvent = new DelayDeparture(queueIdx, customer);
                departureEvent.schedule(serviceWallDelay(queueIdx, classId, serviceTime));
                if (hasRemovalSignals) {
                    long jobId = nextDelayJobId++;
                    delayJobs.put(Long.valueOf(jobId),
                            new DelayJob(queueIdx, customer, departureEvent));
                }
            } else {
                // Queue node: check if server is available (heterogeneous-aware)
                ServerSelection serverSelection = findFreeServerForClass(queueIdx, classId);
                int freeServer = serverSelection.serverId;
                int serverTypeId = serverSelection.serverTypeId;
                if (freeServer >= 0) {
                    if (hasSetupDelayoff[queueIdx]) {
                        ServerState state = serverState[queueIdx][freeServer];
                        if (state == ServerState.ACTIVE) {
                            startService(queueIdx, freeServer, customer, serverTypeId);
                        } else if (state == ServerState.DELAYOFF) {
                            cancelDelayoff(queueIdx, freeServer);
                            startService(queueIdx, freeServer, customer, serverTypeId);
                        } else if (state == ServerState.OFF) {
                            startServerSetup(queueIdx, freeServer, classId);
                            waitQueues[queueIdx].add(customer);
                            scheduleRenegingEvent(queueIdx, customer);
                        } else if (state == ServerState.SETUP) {
                            waitQueues[queueIdx].add(customer);
                            scheduleRenegingEvent(queueIdx, customer);
                        }
                    } else {
                        markServerBusy(queueIdx, freeServer, serverTypeId);
                        customersInService[queueIdx]++;
                        customer.assignedServerType = serverTypeId;

                        updateBusyStats(queueIdx, classId);
                        currentBusyServers[queueIdx][classId]++;
                        lastBusyUpdateTime[queueIdx][classId] = ssjSim.time();

                        double serviceTime;
                        if (customer.serviceTime > 0) {
                            serviceTime = customer.serviceTime;
                        } else {
                            serviceTime = generateHeteroServiceTime(queueIdx, classId, serverTypeId);
                        }
                        if (serviceTime < 0) {
                            throw new RuntimeException(
                                    "LDES: Service time is negative (" + serviceTime
                                            + ") for station " + queueIdx + ", class "
                                            + classId + " at time " + ssjSim.time());
                        }
                        Event departureEvent = new Departure(queueIdx, freeServer, customer);
                        departureEvent.schedule(serviceWallDelay(queueIdx, classId, serviceTime));
                        if ((hasCd || hasLld) && sdDepartureEvents != null) {
                            sdDepartureEvents[queueIdx].put(Integer.valueOf(freeServer), departureEvent);
                            sdInServiceCustomers[queueIdx].put(Integer.valueOf(freeServer), customer);
                        }
                        if (hasRemovalSignals) {
                            inServiceJobs.put(new IntPair(queueIdx, freeServer),
                                    new InServiceJob(customer, departureEvent));
                        }
                    }
                } else {
                    waitQueues[queueIdx].add(customer);
                    scheduleRenegingEvent(queueIdx, customer);
                }
            }

            logEvent("ARRIVAL", serviceStations.get(queueIdx), classId,
                    currentQueueLength[queueIdx][classId],
                    currentBusyServers[queueIdx][classId]);

            return true;
        }

        /**
         * Handles arrival of a negative signal at a service node
         * (G-network semantics).  A negative signal removes jobs from the
         * queue (if present) according to its removal distribution and policy,
         * and is then annihilated: it never joins the station and never
         * continues along its routing chain.  A signal that finds the station
         * empty is lost in the same way.  This is Gelenbe's semantics and the
         * one implemented by SolverCTMC (State.afterEventStationSignal) and by
         * SolverMAM; routing the signal onward would make a single signal fire
         * once per downstream station.
         */
        private void handleNegativeSignalArrival(int queueIdx, int signalClassId,
                                                 double systemArrivalTime) {
            int stationIdx = serviceStations.get(queueIdx);

            int totalJobs = getTotalCustomersAtStation(queueIdx);

            if (totalJobs == 0) {
                logEvent("SIGNAL_MISS", stationIdx, signalClassId, 0, 0);
            } else {
                int numToRemove;
                if (isCatastropheSignal[signalClassId]) {
                    numToRemove = totalJobs;
                } else {
                    DiscreteDistribution dist = signalRemovalDist[signalClassId];
                    if (dist != null) {
                        int sampled = (int) dist.sample(1, siroRng)[0];
                        numToRemove = Math.min(sampled, totalJobs);
                    } else {
                        numToRemove = 1;
                    }
                }

                RemovalPolicy policy = signalRemovalPolicy[signalClassId];
                if (policy == null) {
                    policy = RemovalPolicy.RANDOM;
                }

                int removedCount = 0;
                for (int rep = 0; rep < numToRemove; rep++) {
                    if (getTotalCustomersAtStation(queueIdx) > 0) {
                        int removedClassId = removeJobBySignalWithPolicy(queueIdx, policy);
                        if (removedClassId >= 0) {
                            removedCount++;
                            logEvent("SIGNAL_KILL", stationIdx, signalClassId,
                                    getTotalCustomersAtStation(queueIdx), removedClassId);
                        }
                    }
                }

                if (isCatastropheSignal[signalClassId] && removedCount > 0) {
                    logEvent("CATASTROPHE", stationIdx, signalClassId, removedCount, 0);
                }
            }

            // The signal is annihilated here: account for it as a system
            // completion rather than routing it to the next station.
            if (warmupDone) {
                double respTime = ssjSim.time() - systemArrivalTime;
                systemResponseTimeTally[signalClassId].add(respTime);
                systemCompletedCustomers[signalClassId]++;
            }
        }

        /**
         * Handles arrival of a REPLY signal at a queue for synchronous call
         * semantics.
         */
        private void handleReplySignalArrival(int queueIdx, int signalClassId,
                                              long replyJobId, double systemArrivalTime) {
            handleReplySignalArrival(queueIdx, signalClassId, replyJobId, systemArrivalTime, null);
        }

        private void handleReplySignalArrival(int queueIdx, int signalClassId,
                                              long replyJobId, double systemArrivalTime,
                                              ForkedJob carriedFork) {
            int stationIdx = serviceStations.get(queueIdx);

            // Count the Reply signal as a completion at this queue
            if (warmupDone) {
                completedCustomers[queueIdx][signalClassId]++;
            }

            // Look up pending reply by job ID
            PendingReply pendingReply = pendingReplyMap.remove(Long.valueOf(replyJobId));

            Long parentJobId = replyParentJobId.remove(Long.valueOf(replyJobId));
            long outerJobId = (parentJobId == null) ? -1L : parentJobId.longValue();

            if (pendingReply == null) {
                logEvent("REPLY_MISS", stationIdx, signalClassId, 0, (int) replyJobId);
                routeSignalToNextDestination(queueIdx, signalClassId, systemArrivalTime, outerJobId, carriedFork);
                return;
            }

            int blockedQueueIdx = pendingReply.queueIdx;
            int blockedServerId = pendingReply.serverId;
            int blockedStationIdx = serviceStations.get(blockedQueueIdx);
            int originalClassId = pendingReply.originalClassId;

            // Calculate and accumulate blocking time
            double blockingDuration = ssjSim.time() - pendingReply.blockStartTime;
            totalBlockingTime[blockedQueueIdx][originalClassId] += blockingDuration;

            // Update queue stats before decrementing blocked count
            updateQueueStats(blockedQueueIdx, originalClassId);

            // Unblock server
            serverBlocked[blockedQueueIdx][blockedServerId] = false;
            currentBlockedServers[blockedQueueIdx][originalClassId]--;

            logEvent("REPLY_UNBLOCK", blockedStationIdx, originalClassId,
                    blockedServerId, (int) replyJobId);

            // Unblock server and update stats
            markServerIdle(blockedQueueIdx, blockedServerId);
            customersInService[blockedQueueIdx]--;

            // Start next waiting customer (if any)
            if (!waitQueues[blockedQueueIdx].isEmpty()) {
                Customer nextCustomer = waitQueues[blockedQueueIdx].poll();
                int nextClassId = nextCustomer.classId;

                cancelRenegingIfScheduled(blockedQueueIdx, nextCustomer);

                ServerSelection serverSelection = findFreeServerForClass(blockedQueueIdx, nextClassId);
                if (serverSelection.serverId >= 0) {
                    markServerBusy(blockedQueueIdx, serverSelection.serverId,
                            serverSelection.serverTypeId);
                    customersInService[blockedQueueIdx]++;
                    nextCustomer.assignedServerType = serverSelection.serverTypeId;

                    updateBusyStats(blockedQueueIdx, nextClassId);
                    currentBusyServers[blockedQueueIdx][nextClassId]++;

                    double serviceTime;
                    if (nextCustomer.serviceTime > 0) {
                        serviceTime = nextCustomer.serviceTime;
                    } else {
                        serviceTime = generateHeteroServiceTime(blockedQueueIdx,
                                nextClassId, serverSelection.serverTypeId);
                    }
                    Event departureEvent = new Departure(blockedQueueIdx,
                            serverSelection.serverId, nextCustomer);
                    departureEvent.schedule(serviceWallDelay(blockedQueueIdx, nextClassId, serviceTime));

                    if (hasRemovalSignals) {
                        inServiceJobs.put(new IntPair(blockedQueueIdx, serverSelection.serverId),
                                new InServiceJob(nextCustomer, departureEvent));
                    }

                    logEvent("UNBLOCK_SERVICE", blockedStationIdx, nextClassId,
                            currentQueueLength[blockedQueueIdx][nextClassId],
                            currentBusyServers[blockedQueueIdx][nextClassId]);
                } else {
                    waitQueues[blockedQueueIdx].add(nextCustomer);
                }
            }

            // The REPLY signal continues routing to its destination, carrying
            // the outer call it was nested in so that reply can be paired too.
            routeSignalToNextDestination(queueIdx, signalClassId, systemArrivalTime, outerJobId, carriedFork);
        }

        /**
         * Removes one job uniformly at random from all jobs at the station
         * (waiting + in service).  Implements G-network semantics.
         */
        private int removeJobBySignal(int queueIdx) {
            SchedStrategy strategy = schedStrategies[queueIdx];

            if (isDelayNode.get(queueIdx).booleanValue()) {
                return removeJobFromDelayBySignal(queueIdx);
            } else if (isPSScheduling(strategy)) {
                return removeJobFromPSBySignal(queueIdx);
            } else if (isPreemptiveScheduling[queueIdx]) {
                return removeJobFromPreemptiveBySignal(queueIdx);
            } else if (isPollingStation[queueIdx]) {
                return removeJobFromPollingBySignal(queueIdx);
            } else {
                return removeJobFromStandardQueueBySignal(queueIdx);
            }
        }

        /**
         * Removes a job from the station using the specified removal policy.
         */
        private int removeJobBySignalWithPolicy(int queueIdx, RemovalPolicy policy) {
            SchedStrategy strategy = schedStrategies[queueIdx];

            if (isDelayNode.get(queueIdx).booleanValue()) {
                return removeJobFromDelayBySignal(queueIdx);
            } else if (isPSScheduling(strategy)) {
                return removeJobFromPSBySignal(queueIdx);
            } else if (isPreemptiveScheduling[queueIdx]) {
                return removeJobFromPreemptiveBySignal(queueIdx);
            } else if (isPollingStation[queueIdx]) {
                return removeJobFromPollingBySignal(queueIdx);
            } else {
                return removeJobFromStandardQueueBySignalWithPolicy(queueIdx, policy);
            }
        }

        /**
         * Removes a job from a standard queue using the specified policy.
         */
        private int removeJobFromStandardQueueBySignalWithPolicy(int queueIdx,
                                                                 RemovalPolicy policy) {
            int waitingCount = waitQueues[queueIdx].size();
            int inServiceCount = customersInService[queueIdx];
            int totalJobs = waitingCount + inServiceCount;

            if (totalJobs == 0) return -1;

            if (policy == RemovalPolicy.RANDOM) {
                int victimIdx = (int) (siroRng.nextDouble() * totalJobs);
                if (victimIdx < 0) victimIdx = 0;
                if (victimIdx > totalJobs - 1) victimIdx = totalJobs - 1;
                if (victimIdx < waitingCount) {
                    return removeFromWaitingQueue(queueIdx, victimIdx);
                } else {
                    return removeFromInService(queueIdx, victimIdx - waitingCount);
                }
            } else if (policy == RemovalPolicy.FCFS) {
                if (waitingCount > 0) {
                    return removeFromWaitingQueueOldest(queueIdx);
                } else {
                    return removeFromInServiceOldest(queueIdx);
                }
            } else if (policy == RemovalPolicy.LCFS) {
                if (waitingCount > 0) {
                    return removeFromWaitingQueueNewest(queueIdx);
                } else {
                    return removeFromInServiceNewest(queueIdx);
                }
            }
            return -1;
        }

        /**
         * Removes the oldest job from the waiting queue (FCFS removal policy).
         */
        private int removeFromWaitingQueueOldest(int queueIdx) {
            List<Customer> waitList = new ArrayList<Customer>(waitQueues[queueIdx]);
            if (waitList.isEmpty()) return -1;

            int oldestIdx = 0;
            double oldestTime = waitList.get(0).queueArrivalTime;
            for (int i = 1; i < waitList.size(); i++) {
                if (waitList.get(i).queueArrivalTime < oldestTime) {
                    oldestTime = waitList.get(i).queueArrivalTime;
                    oldestIdx = i;
                }
            }

            return removeFromWaitingQueue(queueIdx, oldestIdx);
        }

        /**
         * Removes the newest job from the waiting queue (LCFS removal policy).
         */
        private int removeFromWaitingQueueNewest(int queueIdx) {
            List<Customer> waitList = new ArrayList<Customer>(waitQueues[queueIdx]);
            if (waitList.isEmpty()) return -1;

            int newestIdx = 0;
            double newestTime = waitList.get(0).queueArrivalTime;
            for (int i = 1; i < waitList.size(); i++) {
                if (waitList.get(i).queueArrivalTime > newestTime) {
                    newestTime = waitList.get(i).queueArrivalTime;
                    newestIdx = i;
                }
            }

            return removeFromWaitingQueue(queueIdx, newestIdx);
        }

        /**
         * Removes the oldest job from in-service (FCFS removal policy).
         */
        private int removeFromInServiceOldest(int queueIdx) {
            double oldestTime = Double.MAX_VALUE;
            int oldestServerId = -1;

            for (int sid = 0; sid < numServers[queueIdx]; sid++) {
                if (serverBusy[queueIdx][sid]) {
                    InServiceJob inServiceJob = inServiceJobs.get(new IntPair(queueIdx, sid));
                    if (inServiceJob != null
                            && inServiceJob.customer.systemArrivalTime < oldestTime) {
                        oldestTime = inServiceJob.customer.systemArrivalTime;
                        oldestServerId = sid;
                    }
                }
            }

            if (oldestServerId >= 0) {
                return removeFromInServiceById(queueIdx, oldestServerId);
            }
            return -1;
        }

        /**
         * Removes the newest job from in-service (LCFS removal policy).
         */
        private int removeFromInServiceNewest(int queueIdx) {
            double newestTime = -Double.MAX_VALUE;
            int newestServerId = -1;

            for (int sid = 0; sid < numServers[queueIdx]; sid++) {
                if (serverBusy[queueIdx][sid]) {
                    InServiceJob inServiceJob = inServiceJobs.get(new IntPair(queueIdx, sid));
                    if (inServiceJob != null
                            && inServiceJob.customer.systemArrivalTime > newestTime) {
                        newestTime = inServiceJob.customer.systemArrivalTime;
                        newestServerId = sid;
                    }
                }
            }

            if (newestServerId >= 0) {
                return removeFromInServiceById(queueIdx, newestServerId);
            }
            return -1;
        }

        /**
         * Removes a job in-service at a specific server ID.
         */
        private int removeFromInServiceById(int queueIdx, int serverId) {
            if (!serverBusy[queueIdx][serverId]) return -1;

            InServiceJob inServiceJob = inServiceJobs.remove(new IntPair(queueIdx, serverId));
            if (inServiceJob == null) return -1;
            inServiceJob.departureEvent.cancel();

            int classId = inServiceJob.customer.classId;

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
                // removed jobs exit the region for real: release blocked customers
                tryReleaseBlockedCustomers(regionIdx);
            }

            updateBusyStats(queueIdx, classId);
            currentBusyServers[queueIdx][classId]--;

            markServerIdle(queueIdx, serverId);
            customersInService[queueIdx]--;

            // Start service for next customer in queue (if any)
            if (!waitQueues[queueIdx].isEmpty()) {
                Customer nextCustomer = waitQueues[queueIdx].poll();
                int nextClassId = nextCustomer.classId;

                cancelRenegingIfScheduled(queueIdx, nextCustomer);

                ServerSelection serverSelection = findFreeServerForClass(queueIdx, nextClassId);
                if (serverSelection.serverId >= 0) {
                    markServerBusy(queueIdx, serverSelection.serverId, serverSelection.serverTypeId);
                    customersInService[queueIdx]++;
                    nextCustomer.assignedServerType = serverSelection.serverTypeId;

                    updateBusyStats(queueIdx, nextClassId);
                    currentBusyServers[queueIdx][nextClassId]++;

                    double serviceTime;
                    if (nextCustomer.serviceTime > 0) {
                        serviceTime = nextCustomer.serviceTime;
                    } else {
                        serviceTime = generateHeteroServiceTime(queueIdx, nextClassId,
                                serverSelection.serverTypeId);
                    }
                    Event departureEvent = new Departure(queueIdx, serverSelection.serverId,
                            nextCustomer);
                    departureEvent.schedule(serviceWallDelay(queueIdx, nextClassId, serviceTime));
                    inServiceJobs.put(new IntPair(queueIdx, serverSelection.serverId),
                            new InServiceJob(nextCustomer, departureEvent));
                } else {
                    waitQueues[queueIdx].add(nextCustomer);
                }
            }

            return classId;
        }

        /**
         * Removes a random job from a standard FCFS/LCFS/SIRO queue.
         */
        private int removeJobFromStandardQueueBySignal(int queueIdx) {
            int waitingCount = waitQueues[queueIdx].size();
            int inServiceCount = customersInService[queueIdx];
            int totalJobs = waitingCount + inServiceCount;

            if (totalJobs == 0) return -1;

            int victimIdx = (int) (siroRng.nextDouble() * totalJobs);
            if (victimIdx < 0) victimIdx = 0;
            if (victimIdx > totalJobs - 1) victimIdx = totalJobs - 1;

            if (victimIdx < waitingCount) {
                return removeFromWaitingQueue(queueIdx, victimIdx);
            } else {
                int serverIdx = victimIdx - waitingCount;
                return removeFromInService(queueIdx, serverIdx);
            }
        }

        /**
         * Removes a job at the given index from the waiting queue.
         */
        private int removeFromWaitingQueue(int queueIdx, int victimIdx) {
            List<Customer> waitList = new ArrayList<Customer>(waitQueues[queueIdx]);
            if (victimIdx >= waitList.size()) return -1;

            Customer victim = waitList.get(victimIdx);
            waitQueues[queueIdx].remove(victim);

            int classId = victim.classId;

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
                // removed jobs exit the region for real: release blocked customers
                tryReleaseBlockedCustomers(regionIdx);
            }

            return classId;
        }

        /**
         * Removes a job in service at the given server index and cancels its
         * departure event.
         */
        private int removeFromInService(int queueIdx, int serverIdx) {
            int busyCount = 0;
            int targetServerId = -1;
            for (int sid = 0; sid < numServers[queueIdx]; sid++) {
                if (serverBusy[queueIdx][sid]) {
                    if (busyCount == serverIdx) {
                        targetServerId = sid;
                        break;
                    }
                    busyCount++;
                }
            }

            if (targetServerId < 0) return -1;

            InServiceJob inServiceJob = inServiceJobs.remove(new IntPair(queueIdx, targetServerId));
            if (inServiceJob == null) return -1;
            inServiceJob.departureEvent.cancel();

            int classId = inServiceJob.customer.classId;

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
                // removed jobs exit the region for real: release blocked customers
                tryReleaseBlockedCustomers(regionIdx);
            }

            updateBusyStats(queueIdx, classId);
            currentBusyServers[queueIdx][classId]--;

            markServerIdle(queueIdx, targetServerId);
            customersInService[queueIdx]--;

            // Start service for next customer in queue (if any)
            if (!waitQueues[queueIdx].isEmpty()) {
                Customer nextCustomer = waitQueues[queueIdx].poll();
                int nextClassId = nextCustomer.classId;

                cancelRenegingIfScheduled(queueIdx, nextCustomer);

                ServerSelection serverSelection = findFreeServerForClass(queueIdx, nextClassId);
                if (serverSelection.serverId >= 0) {
                    markServerBusy(queueIdx, serverSelection.serverId, serverSelection.serverTypeId);
                    customersInService[queueIdx]++;
                    nextCustomer.assignedServerType = serverSelection.serverTypeId;

                    updateBusyStats(queueIdx, nextClassId);
                    currentBusyServers[queueIdx][nextClassId]++;

                    double serviceTime;
                    if (nextCustomer.serviceTime > 0) {
                        serviceTime = nextCustomer.serviceTime;
                    } else {
                        serviceTime = generateHeteroServiceTime(queueIdx, nextClassId,
                                serverSelection.serverTypeId);
                    }
                    Event departureEvent = new Departure(queueIdx, serverSelection.serverId,
                            nextCustomer);
                    departureEvent.schedule(serviceWallDelay(queueIdx, nextClassId, serviceTime));
                    inServiceJobs.put(new IntPair(queueIdx, serverSelection.serverId),
                            new InServiceJob(nextCustomer, departureEvent));
                } else {
                    waitQueues[queueIdx].add(nextCustomer);
                }
            }

            return classId;
        }

        /**
         * Removes a random job from a PS (Processor Sharing) queue.
         */
        private int removeJobFromPSBySignal(int queueIdx) {
            List<PSCustomer> jobs = psJobsInService[queueIdx];
            if (jobs.isEmpty()) return -1;

            updatePSRemainingWork(queueIdx, ssjSim.time());

            int victimIdx = (int) (siroRng.nextDouble() * jobs.size());
            if (victimIdx < 0) victimIdx = 0;
            if (victimIdx > jobs.size() - 1) victimIdx = jobs.size() - 1;
            PSCustomer victim = jobs.get(victimIdx);

            if (victim.scheduledDepartureEvent != null) {
                victim.scheduledDepartureEvent.cancel();
            }

            jobs.remove(victimIdx);

            int classId = victim.classId;

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
                // removed jobs exit the region for real: release blocked customers
                tryReleaseBlockedCustomers(regionIdx);
            }

            updatePSBusyStats(queueIdx);

            rescheduleAllPSDepartures(queueIdx);

            return classId;
        }

        /**
         * Removes a random job from a Delay node (infinite servers).
         */
        private int removeJobFromDelayBySignal(int queueIdx) {
            List<Map.Entry<Long, DelayJob>> queueDelayJobs = new ArrayList<Map.Entry<Long, DelayJob>>();
            for (Map.Entry<Long, DelayJob> entry : delayJobs.entrySet()) {
                if (entry.getValue().queueIdx == queueIdx) {
                    queueDelayJobs.add(entry);
                }
            }
            if (queueDelayJobs.isEmpty()) return -1;

            int victimIdx = (int) (siroRng.nextDouble() * queueDelayJobs.size());
            if (victimIdx < 0) victimIdx = 0;
            if (victimIdx > queueDelayJobs.size() - 1) victimIdx = queueDelayJobs.size() - 1;
            Map.Entry<Long, DelayJob> victimEntry = queueDelayJobs.get(victimIdx);
            DelayJob victim = victimEntry.getValue();

            victim.departureEvent.cancel();

            delayJobs.remove(victimEntry.getKey());

            int classId = victim.customer.classId;

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
                // removed jobs exit the region for real: release blocked customers
                tryReleaseBlockedCustomers(regionIdx);
            }

            updateBusyStats(queueIdx, classId);
            currentBusyServers[queueIdx][classId]--;

            customersInService[queueIdx]--;

            return classId;
        }

        /**
         * Removes a random job from a preemptive LCFS queue.
         */
        private int removeJobFromPreemptiveBySignal(int queueIdx) {
            List<PreemptiveCustomer> jobs = preemptiveJobsInService[queueIdx];
            if (jobs.isEmpty()) return -1;

            int victimIdx = (int) (siroRng.nextDouble() * jobs.size());
            if (victimIdx < 0) victimIdx = 0;
            if (victimIdx > jobs.size() - 1) victimIdx = jobs.size() - 1;
            PreemptiveCustomer victim = jobs.get(victimIdx);

            if (victim.scheduledDepartureEvent != null) {
                victim.scheduledDepartureEvent.cancel();
            }

            jobs.remove(victimIdx);

            int classId = victim.classId;

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
                // removed jobs exit the region for real: release blocked customers
                tryReleaseBlockedCustomers(regionIdx);
            }

            updateBusyStats(queueIdx, classId);
            currentBusyServers[queueIdx][classId]--;

            markServerIdle(queueIdx, victim.serverId);
            customersInService[queueIdx]--;

            return classId;
        }

        /**
         * Removes a random job from a polling queue.
         */
        private int removeJobFromPollingBySignal(int queueIdx) {
            // Collect all jobs across all class queues
            List<int[]> allClassIdxs = new ArrayList<int[]>();
            List<Customer> allJobs = new ArrayList<Customer>();
            for (int k = 0; k < numClasses; k++) {
                for (Customer customer : pollingQueues[queueIdx][k]) {
                    allClassIdxs.add(new int[] { k });
                    allJobs.add(customer);
                }
            }

            if (allJobs.isEmpty()) return -1;

            int victimIdx = (int) (siroRng.nextDouble() * allJobs.size());
            if (victimIdx < 0) victimIdx = 0;
            if (victimIdx > allJobs.size() - 1) victimIdx = allJobs.size() - 1;
            int classId = allClassIdxs.get(victimIdx)[0];
            Customer victim = allJobs.get(victimIdx);

            pollingQueues[queueIdx][classId].remove(victim);

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
                // removed jobs exit the region for real: release blocked customers
                tryReleaseBlockedCustomers(regionIdx);
            }

            return classId;
        }

        /**
         * Routes a signal to its next destination after processing at the
         * current node.
         */
        /**
         * Job id to give a customer entering class destClassId while carrying
         * inherited. A class that expects a reply starts a new synchronous
         * call and so needs an id of its own; an id already carried belongs to
         * an outer call still in flight and is remembered as its parent.
         */
        private long mintCallJobId(int destClassId, long inherited, int srcClassId) {
            if (destClassId < 0 || synchCallReplyClass[destClassId] < 0) {
                return inherited;
            }
            if (destClassId == srcClassId && inherited >= 0) {
                // Same class at the next station: the job is still inside the
                // call it already carries, not issuing a new one. Only a class
                // switch into a blocking class starts a call.
                return inherited;
            }
            long minted = nextJobId++;
            if (inherited >= 0) {
                replyParentJobId.put(Long.valueOf(minted), Long.valueOf(inherited));
            }
            return minted;
        }

        /**
         * Injects the spawn-on-completion continuation of a completing job
         * (sn.classspawn, LQN phase-2): a fresh job of the mapped class
         * arrives at the same station. Called after regionJobLeave and
         * before regionExitCompleted, so the continuation takes over the
         * region slot the completing job just freed, ahead of any customer
         * waiting in the region blocked FIFO.
         */
        private void maybeSpawnOnCompletion(int queueIdx, int classId, Customer customer) {
            int spawnCls = spawnTargetOf(classId);
            if (spawnCls < 0) {
                return;
            }
            ForkedJob inherited = null;
            if (customer != null && spawnReachesJoin(spawnCls)) {
                // A continuation aimed at a Join stands in for the trigger job
                // at the fork it belongs to (LQN phase-2 at an AND-join branch
                // tail): hand the fork identity to the clone and strip it from
                // the trigger, whose reply signal must not carry it upstream.
                IntDoubleKey key = new IntDoubleKey(queueIdx, customer.queueArrivalTime,
                        customer.randomRank);
                inherited = forkedCustomerMap.remove(key);
                if (inherited == null) {
                    inherited = customer.forkedJob;
                }
                customer.forkedJob = null;
            }
            spawnContinuation(queueIdx, spawnCls, inherited);
        }

        /** Key-based overload for customer types whose fork record lives only
         *  in forkedCustomerMap (preemptive queues). */
        private void maybeSpawnOnCompletion(int queueIdx, int classId,
                                            double queueArrivalTime, double randomRank) {
            int spawnCls = spawnTargetOf(classId);
            if (spawnCls < 0) {
                return;
            }
            ForkedJob inherited = null;
            if (spawnReachesJoin(spawnCls)) {
                inherited = forkedCustomerMap.remove(
                        new IntDoubleKey(queueIdx, queueArrivalTime, randomRank));
            }
            spawnContinuation(queueIdx, spawnCls, inherited);
        }

        /** PS-family overload: the fork record sits in a final field, so it is
         *  handed over without being stripped from the departing customer. */
        private void maybeSpawnOnCompletion(int queueIdx, int classId, ForkedJob candidateFork) {
            int spawnCls = spawnTargetOf(classId);
            if (spawnCls < 0) {
                return;
            }
            ForkedJob inherited = spawnReachesJoin(spawnCls) ? candidateFork : null;
            spawnContinuation(queueIdx, spawnCls, inherited);
        }

        private int spawnTargetOf(int classId) {
            if (spawnClassOf == null || classId < 0 || classId >= spawnClassOf.length) {
                return -1;
            }
            int spawnCls = spawnClassOf[classId];
            return (spawnCls < 0 || spawnCls >= numClasses) ? -1 : spawnCls;
        }

        private boolean spawnReachesJoin(int spawnCls) {
            if (joinNodes.isEmpty() || sn.rtnodes == null) {
                return false;
            }
            if (spawnJoinReach == null) {
                spawnJoinReach = new byte[numClasses];
            }
            byte v = spawnJoinReach[spawnCls];
            if (v == 0) {
                boolean reach = false;
                int R = numClasses;
                outer:
                for (Integer jn : joinNodes) {
                    for (int from = 0; from < numNodes; from++) {
                        for (int tc = 0; tc < R; tc++) {
                            if (sn.rtnodes.get(from * R + spawnCls, jn.intValue() * R + tc) > 0) {
                                reach = true;
                                break outer;
                            }
                        }
                    }
                }
                v = reach ? (byte) 2 : (byte) 1;
                spawnJoinReach[spawnCls] = v;
            }
            return v == 2;
        }

        private void spawnContinuation(int queueIdx, int spawnCls, ForkedJob inherited) {
            double now = ssjSim.time();
            long jobId = mintCallJobId(spawnCls, -1L, -1);
            Customer spawned = new Customer(spawnCls, classPrio[spawnCls], now, now,
                    siroRng.nextDouble(), -1.0, jobId, Double.POSITIVE_INFINITY, -1, null);
            logEvent("SPAWN", serviceStations.get(queueIdx), spawnCls,
                    currentQueueLength[queueIdx][spawnCls], 0);
            if (inherited != null) {
                ForkedJob nextFork = new ForkedJob(inherited.forkJobId, inherited.parentJobId,
                        spawnCls, inherited.priority, inherited.systemArrivalTime,
                        now, inherited.randomRank);
                arriveAtQueueForked(queueIdx, spawned, nextFork);
            } else {
                arriveAtQueue(queueIdx, spawned);
            }
        }

        private void routeSignalToNextDestination(int queueIdx, int signalClassId,
                                                  double systemArrivalTime) {
            routeSignalToNextDestination(queueIdx, signalClassId, systemArrivalTime, -1L);
        }

        private void routeSignalToNextDestination(int queueIdx, int signalClassId,
                                                  double systemArrivalTime, long carriedJobId) {
            routeSignalToNextDestination(queueIdx, signalClassId, systemArrivalTime, carriedJobId, null);
        }

        private void routeSignalToNextDestination(int queueIdx, int signalClassId,
                                                  double systemArrivalTime, long carriedJobId,
                                                  ForkedJob carriedFork) {
            int currentNode = serviceNodes.get(queueIdx).intValue();
            RoutingResult routingResult = selectDestination(currentNode, signalClassId);
            int destNode = routingResult.destNode;
            int destClassId = routingResult.destClassId;

            if (destNode < 0) {
                // No destination - signal ends here
                return;
            }
            if (sinkNodes.contains(Integer.valueOf(destNode))) {
                if (warmupDone) {
                    double respTime = ssjSim.time() - systemArrivalTime;
                    systemResponseTimeTally[signalClassId].add(respTime);
                    systemCompletedCustomers[signalClassId]++;
                }
                return;
            }
            if (forkNodes.contains(Integer.valueOf(destNode))) {
                long parentJobId = nextJobId++;
                handleForkArrival(destNode, parentJobId, destClassId, systemArrivalTime);
                return;
            }
            // Route to next queue
            int nextQueueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
            if (nextQueueIdx >= 0) {
                if (isNegativeSignal[destClassId] || isCatastropheSignal[destClassId]) {
                    handleNegativeSignalArrival(nextQueueIdx, destClassId, systemArrivalTime);
                } else if (hasReplySignals && isReplySignal[destClassId]) {
                    // A reply that itself class-switches into another reply is
                    // the return leg of a nested call: it must unblock the
                    // outer caller rather than queue there as a customer.
                    handleReplySignalArrival(nextQueueIdx, destClassId, carriedJobId, systemArrivalTime, carriedFork);
                } else {
                    long jobId = mintCallJobId(destClassId, carriedJobId, signalClassId);
                    Customer customer = new Customer(
                            destClassId, classPrio[destClassId],
                            systemArrivalTime, ssjSim.time(),
                            siroRng.nextDouble(),
                            -1.0, jobId,
                            Double.POSITIVE_INFINITY, -1, null);
                    if (carriedFork != null) {
                        // The signal was the return leg of a fork branch: the
                        // continuation must reach the Join as a forked task or
                        // the siblings are never matched.
                        ForkedJob nextFork = new ForkedJob(
                                carriedFork.forkJobId,
                                carriedFork.parentJobId,
                                destClassId,
                                carriedFork.priority,
                                carriedFork.systemArrivalTime,
                                ssjSim.time(),
                                carriedFork.randomRank);
                        arriveAtQueueForked(nextQueueIdx, customer, nextFork);
                    } else {
                        arriveAtQueue(nextQueueIdx, customer);
                    }
                }
            }
        }

        /**
         * Customer arrives at a preemptive service node.  Handles preemption
         * logic for LCFSPR/LCFSPI, FCFSPR/FCFSPI, and size-based preemptive
         * scheduling.
         */
        private boolean arriveAtPreemptiveLCFSQueue(int queueIdx, Customer customer) {
            int classId = customer.classId;
            SchedStrategy strategy = schedStrategies[queueIdx];

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]++;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobEnter(regionIdx, classId);
                updateRegionArrivalTracking(regionIdx, classId);
            }

            int totalAtStation = getTotalCustomersAtStation(queueIdx);
            if (totalAtStation > maxQueueLengthReached) {
                maxQueueLengthReached = totalAtStation;
            }

            if (isSizeBasedPreemptiveScheduling(strategy)) {
                customer.serviceTime = generateServiceTime(queueIdx, classId);

                boolean shouldPreempt = shouldPreemptForSizeBasedPolicy(queueIdx, customer, strategy);
                if (shouldPreempt) {
                    PreemptiveCustomer victimJob =
                            findJobToPreemptForSizeBasedPolicy(queueIdx, customer, strategy);
                    if (victimJob != null) {
                        preemptJob(queueIdx, victimJob);
                        startPreemptiveService(queueIdx, customer, victimJob.serverId,
                                victimJob.assignedServerType);
                    } else {
                        waitQueues[queueIdx].add(customer);
                    }
                } else {
                    ServerSelection serverSelection = findFreeServerForClass(queueIdx, classId);
                    if (serverSelection.serverId >= 0) {
                        startPreemptiveService(queueIdx, customer, serverSelection.serverId,
                                serverSelection.serverTypeId);
                    } else {
                        waitQueues[queueIdx].add(customer);
                    }
                }
            } else if (isPreemptiveFCFSScheduling(strategy)) {
                boolean shouldPreempt = shouldPreemptForFCFS(queueIdx, customer);

                if (shouldPreempt) {
                    PreemptiveCustomer victimJob = findJobToPreemptFCFS(queueIdx, customer, strategy);
                    if (victimJob != null) {
                        preemptJob(queueIdx, victimJob);
                        startPreemptiveService(queueIdx, customer, victimJob.serverId,
                                victimJob.assignedServerType);
                    } else {
                        waitQueues[queueIdx].add(customer);
                    }
                } else {
                    ServerSelection serverSelection = findFreeServerForClass(queueIdx, classId);
                    if (serverSelection.serverId >= 0) {
                        startPreemptiveService(queueIdx, customer, serverSelection.serverId,
                                serverSelection.serverTypeId);
                    } else {
                        waitQueues[queueIdx].add(customer);
                    }
                }
            } else {
                // LCFS preemptive handling
                boolean shouldPreempt = shouldPreemptForLCFS(queueIdx, customer);

                if (shouldPreempt) {
                    PreemptiveCustomer victimJob = findJobToPreempt(queueIdx, customer, strategy);
                    if (victimJob != null) {
                        preemptJob(queueIdx, victimJob);
                        startPreemptiveService(queueIdx, customer, victimJob.serverId,
                                victimJob.assignedServerType);
                    } else {
                        waitQueues[queueIdx].add(customer);
                    }
                } else {
                    ServerSelection serverSelection = findFreeServerForClass(queueIdx, classId);
                    if (serverSelection.serverId >= 0) {
                        startPreemptiveService(queueIdx, customer, serverSelection.serverId,
                                serverSelection.serverTypeId);
                    } else {
                        waitQueues[queueIdx].add(customer);
                    }
                }
            }

            logEvent("ARRIVAL", stationIdx, classId,
                    currentQueueLength[queueIdx][classId],
                    currentBusyServers[queueIdx][classId]);

            return true;
        }

        /**
         * Determines if an arriving customer should preempt a currently
         * serving job.
         */
        private boolean shouldPreemptForLCFS(int queueIdx, Customer arrivingCustomer) {
            ServerSelection serverSelection = findFreeServerForClass(queueIdx, arrivingCustomer.classId);
            if (serverSelection.serverId >= 0) return false;

            List<PreemptiveCustomer> jobsInService = preemptiveJobsInService[queueIdx];
            if (jobsInService.isEmpty()) return false;

            SchedStrategy strategy = schedStrategies[queueIdx];
            if (strategy == SchedStrategy.LCFSPRPRIO || strategy == SchedStrategy.LCFSPIPRIO) {
                for (PreemptiveCustomer it : jobsInService) {
                    if (it.priority >= arrivingCustomer.priority) return true;
                }
                return false;
            }

            return true;
        }

        /**
         * Finds the job to preempt based on scheduling strategy.
         */
        private PreemptiveCustomer findJobToPreempt(int queueIdx, Customer arrivingCustomer,
                                                    SchedStrategy strategy) {
            List<PreemptiveCustomer> jobsInService = preemptiveJobsInService[queueIdx];

            if (strategy == SchedStrategy.LCFSPR || strategy == SchedStrategy.LCFSPI) {
                PreemptiveCustomer best = null;
                double bestTime = Double.POSITIVE_INFINITY;
                for (PreemptiveCustomer it : jobsInService) {
                    if (it.queueArrivalTime < bestTime) {
                        bestTime = it.queueArrivalTime;
                        best = it;
                    }
                }
                return best;
            }
            if (strategy == SchedStrategy.LCFSPRPRIO || strategy == SchedStrategy.LCFSPIPRIO) {
                List<PreemptiveCustomer> lowerPriorityJobs = new ArrayList<PreemptiveCustomer>();
                for (PreemptiveCustomer it : jobsInService) {
                    if (it.priority > arrivingCustomer.priority) {
                        lowerPriorityJobs.add(it);
                    }
                }
                if (!lowerPriorityJobs.isEmpty()) {
                    PreemptiveCustomer best = null;
                    for (PreemptiveCustomer it : lowerPriorityJobs) {
                        if (best == null) {
                            best = it;
                            continue;
                        }
                        // compareBy({-priority}, {queueArrivalTime}): primary
                        // is descending priority, then ascending arrival time.
                        if (-it.priority < -best.priority) {
                            best = it;
                        } else if (-it.priority == -best.priority
                                && it.queueArrivalTime < best.queueArrivalTime) {
                            best = it;
                        }
                    }
                    return best;
                }
                List<PreemptiveCustomer> samePriorityJobs = new ArrayList<PreemptiveCustomer>();
                for (PreemptiveCustomer it : jobsInService) {
                    if (it.priority == arrivingCustomer.priority) {
                        samePriorityJobs.add(it);
                    }
                }
                if (!samePriorityJobs.isEmpty()) {
                    PreemptiveCustomer best = null;
                    double bestStart = Double.POSITIVE_INFINITY;
                    for (PreemptiveCustomer it : samePriorityJobs) {
                        if (it.serviceStartTime < bestStart) {
                            bestStart = it.serviceStartTime;
                            best = it;
                        }
                    }
                    return best;
                }
                return null;
            }
            return null;
        }

        /**
         * Determines if an arriving customer should preempt a currently
         * serving job under FCFS preemptive scheduling.
         */
        private boolean shouldPreemptForFCFS(int queueIdx, Customer arrivingCustomer) {
            ServerSelection serverSelection = findFreeServerForClass(queueIdx, arrivingCustomer.classId);
            if (serverSelection.serverId >= 0) return false;

            List<PreemptiveCustomer> jobsInService = preemptiveJobsInService[queueIdx];
            if (jobsInService.isEmpty()) return false;

            for (PreemptiveCustomer it : jobsInService) {
                if (it.priority > arrivingCustomer.priority) return true;
            }
            return false;
        }

        /**
         * Finds the job to preempt under FCFS preemptive scheduling.
         */
        private PreemptiveCustomer findJobToPreemptFCFS(int queueIdx, Customer arrivingCustomer,
                                                        SchedStrategy strategy) {
            List<PreemptiveCustomer> jobsInService = preemptiveJobsInService[queueIdx];

            List<PreemptiveCustomer> lowerPriorityJobs = new ArrayList<PreemptiveCustomer>();
            for (PreemptiveCustomer it : jobsInService) {
                if (it.priority > arrivingCustomer.priority) {
                    lowerPriorityJobs.add(it);
                }
            }
            if (lowerPriorityJobs.isEmpty()) return null;

            // maxWithOrNull(compareBy(priority, serviceStartTime)): pick the
            // largest priority (lowest priority class), tie-broken by latest
            // service start time.
            PreemptiveCustomer best = null;
            for (PreemptiveCustomer it : lowerPriorityJobs) {
                if (best == null) {
                    best = it;
                    continue;
                }
                if (it.priority > best.priority) {
                    best = it;
                } else if (it.priority == best.priority
                        && it.serviceStartTime > best.serviceStartTime) {
                    best = it;
                }
            }
            return best;
        }

        /**
         * Checks if arriving customer should trigger SRPT preemption.
         */
        private boolean shouldPreemptForSRPT(int queueIdx, Customer arrivingCustomer) {
            SchedStrategy strategy = schedStrategies[queueIdx];
            if (strategy != SchedStrategy.SRPT && strategy != SchedStrategy.SRPTPRIO) {
                return false;
            }

            ServerSelection serverSelection = findFreeServerForClass(queueIdx, arrivingCustomer.classId);
            if (serverSelection.serverId >= 0) return false;

            List<PreemptiveCustomer> jobsInService = preemptiveJobsInService[queueIdx];
            if (jobsInService.isEmpty()) return false;

            if (strategy == SchedStrategy.SRPTPRIO) {
                List<PreemptiveCustomer> eligibleJobs = new ArrayList<PreemptiveCustomer>();
                for (PreemptiveCustomer it : jobsInService) {
                    if (it.priority >= arrivingCustomer.priority) {
                        eligibleJobs.add(it);
                    }
                }
                if (eligibleJobs.isEmpty()) return false;

                PreemptiveCustomer jobWithMostWork = null;
                double mostWork = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : eligibleJobs) {
                    if (it.remainingServiceWork > mostWork) {
                        mostWork = it.remainingServiceWork;
                        jobWithMostWork = it;
                    }
                }
                return jobWithMostWork != null
                        && arrivingCustomer.serviceTime < jobWithMostWork.remainingServiceWork;
            }

            PreemptiveCustomer jobWithMostWork = null;
            double mostWork = -Double.MAX_VALUE;
            for (PreemptiveCustomer it : jobsInService) {
                if (it.remainingServiceWork > mostWork) {
                    mostWork = it.remainingServiceWork;
                    jobWithMostWork = it;
                }
            }
            return jobWithMostWork != null
                    && arrivingCustomer.serviceTime < jobWithMostWork.remainingServiceWork;
        }

        /**
         * Finds the job to preempt for SRPT - job with maximum remaining work.
         */
        private PreemptiveCustomer findJobToPreemptForSRPT(int queueIdx, Customer arrivingCustomer) {
            SchedStrategy strategy = schedStrategies[queueIdx];
            List<PreemptiveCustomer> jobsInService = preemptiveJobsInService[queueIdx];

            if (strategy == SchedStrategy.SRPTPRIO) {
                PreemptiveCustomer best = null;
                double bestWork = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    if (it.priority >= arrivingCustomer.priority
                            && it.remainingServiceWork > bestWork) {
                        bestWork = it.remainingServiceWork;
                        best = it;
                    }
                }
                return best;
            }

            PreemptiveCustomer best = null;
            double bestWork = -Double.MAX_VALUE;
            for (PreemptiveCustomer it : jobsInService) {
                if (it.remainingServiceWork > bestWork) {
                    bestWork = it.remainingServiceWork;
                    best = it;
                }
            }
            return best;
        }

        /**
         * Determines if an arriving customer should preempt a currently
         * serving job for size-based preemptive policies (SRPT, PSJF, FB,
         * LRPT).
         */
        private boolean shouldPreemptForSizeBasedPolicy(int queueIdx, Customer arrivingCustomer,
                                                        SchedStrategy strategy) {
            ServerSelection serverSelection = findFreeServerForClass(queueIdx, arrivingCustomer.classId);
            if (serverSelection.serverId >= 0) return false;

            List<PreemptiveCustomer> jobsInService = preemptiveJobsInService[queueIdx];
            if (jobsInService.isEmpty()) return false;

            double currentTime = ssjSim.time();

            if (isSRPTScheduling(strategy)) {
                PreemptiveCustomer jobWithMostWork = null;
                double mostWork = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    double key = it.remainingServiceWork - serviceWorkBetween(queueIdx, it.classId, it.serviceStartTime, currentTime);
                    if (key > mostWork) {
                        mostWork = key;
                        jobWithMostWork = it;
                    }
                }
                if (jobWithMostWork != null) {
                    double currentRemaining = jobWithMostWork.remainingServiceWork
                            - serviceWorkBetween(queueIdx, jobWithMostWork.classId, jobWithMostWork.serviceStartTime, currentTime);
                    return arrivingCustomer.serviceTime < currentRemaining;
                }
                return false;
            }
            if (isPSJFScheduling(strategy)) {
                PreemptiveCustomer jobWithLargestOriginal = null;
                double largest = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    if (it.totalServiceRequirement > largest) {
                        largest = it.totalServiceRequirement;
                        jobWithLargestOriginal = it;
                    }
                }
                return jobWithLargestOriginal != null
                        && arrivingCustomer.serviceTime < jobWithLargestOriginal.totalServiceRequirement;
            }
            if (isFBScheduling(strategy)) {
                PreemptiveCustomer jobWithMostProgress = null;
                double mostProgress = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    double key = it.elapsedServiceTime + serviceWorkBetween(queueIdx, it.classId, it.serviceStartTime, currentTime);
                    if (key > mostProgress) {
                        mostProgress = key;
                        jobWithMostProgress = it;
                    }
                }
                if (jobWithMostProgress != null) {
                    double currentAttained = jobWithMostProgress.elapsedServiceTime
                            + serviceWorkBetween(queueIdx, jobWithMostProgress.classId, jobWithMostProgress.serviceStartTime, currentTime);
                    return currentAttained > 0.0;
                }
                return false;
            }
            if (isLRPTScheduling(strategy)) {
                PreemptiveCustomer jobWithLeastWork = null;
                double leastWork = Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    double key = it.remainingServiceWork - serviceWorkBetween(queueIdx, it.classId, it.serviceStartTime, currentTime);
                    if (key < leastWork) {
                        leastWork = key;
                        jobWithLeastWork = it;
                    }
                }
                if (jobWithLeastWork != null) {
                    double currentRemaining = jobWithLeastWork.remainingServiceWork
                            - serviceWorkBetween(queueIdx, jobWithLeastWork.classId, jobWithLeastWork.serviceStartTime, currentTime);
                    return arrivingCustomer.serviceTime > currentRemaining;
                }
                return false;
            }
            if (isFSPScheduling(strategy)) {
                double arrivalVft = computeFSPVirtualFinishTime(queueIdx, arrivingCustomer);
                PreemptiveCustomer worstVftJob = null;
                double worstVft = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    double vft = computeFSPVirtualFinishTimeInService(queueIdx, it);
                    if (vft > worstVft) {
                        worstVft = vft;
                        worstVftJob = it;
                    }
                }
                return worstVftJob != null && arrivalVft < worstVft;
            }
            if (isEDFScheduling(strategy)) {
                // EDF preempts the in-service job with the latest (largest) absolute deadline
                // when an arriving job has an earlier deadline.
                PreemptiveCustomer jobWithLatestDeadline = null;
                double latestDeadline = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    if (it.absoluteDeadline > latestDeadline) {
                        latestDeadline = it.absoluteDeadline;
                        jobWithLatestDeadline = it;
                    }
                }
                return jobWithLatestDeadline != null
                        && arrivingCustomer.absoluteDeadline < jobWithLatestDeadline.absoluteDeadline;
            }
            return false;
        }

        /**
         * Finds the job to preempt for size-based preemptive policies.
         */
        private PreemptiveCustomer findJobToPreemptForSizeBasedPolicy(int queueIdx,
                                                                      Customer arrivingCustomer,
                                                                      SchedStrategy strategy) {
            List<PreemptiveCustomer> jobsInService = preemptiveJobsInService[queueIdx];
            double currentTime = ssjSim.time();

            if (isSRPTScheduling(strategy)) {
                PreemptiveCustomer best = null;
                double bestKey = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    double key = it.remainingServiceWork - serviceWorkBetween(queueIdx, it.classId, it.serviceStartTime, currentTime);
                    if (key > bestKey) {
                        bestKey = key;
                        best = it;
                    }
                }
                return best;
            }
            if (isPSJFScheduling(strategy)) {
                PreemptiveCustomer best = null;
                double bestKey = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    if (it.totalServiceRequirement > bestKey) {
                        bestKey = it.totalServiceRequirement;
                        best = it;
                    }
                }
                return best;
            }
            if (isFBScheduling(strategy)) {
                PreemptiveCustomer best = null;
                double bestKey = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    double key = it.elapsedServiceTime + serviceWorkBetween(queueIdx, it.classId, it.serviceStartTime, currentTime);
                    if (key > bestKey) {
                        bestKey = key;
                        best = it;
                    }
                }
                return best;
            }
            if (isLRPTScheduling(strategy)) {
                PreemptiveCustomer best = null;
                double bestKey = Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    double key = it.remainingServiceWork - serviceWorkBetween(queueIdx, it.classId, it.serviceStartTime, currentTime);
                    if (key < bestKey) {
                        bestKey = key;
                        best = it;
                    }
                }
                return best;
            }
            if (isFSPScheduling(strategy)) {
                PreemptiveCustomer best = null;
                double bestKey = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    double vft = computeFSPVirtualFinishTimeInService(queueIdx, it);
                    if (vft > bestKey) {
                        bestKey = vft;
                        best = it;
                    }
                }
                return best;
            }
            if (isEDFScheduling(strategy)) {
                PreemptiveCustomer best = null;
                double latestDeadline = -Double.MAX_VALUE;
                for (PreemptiveCustomer it : jobsInService) {
                    if (it.absoluteDeadline > latestDeadline) {
                        latestDeadline = it.absoluteDeadline;
                        best = it;
                    }
                }
                return best;
            }
            return null;
        }

        // ==================== Forward-reference stubs (PART 7+) ====================
        //
        // The following methods/classes are called from the routing/preemptive
        // code above but are translated in subsequent chunks (Kotlin lines
        // 9000+).

        /** Forward stub — Kotlin line 11984. */
        private static final class ServerSelection {
            final int serverId;
            final int serverTypeId;
            ServerSelection(int serverId, int serverTypeId) {
                this.serverId = serverId;
                this.serverTypeId = serverTypeId;
            }
        }

        /** Body translated in PART 8 (see {@link #findFreeServerForClassImpl(int, int)}). */
        private ServerSelection findFreeServerForClass(int queueIdx, int classId) {
            return findFreeServerForClassImpl(queueIdx, classId);
        }

        /** Marks a server as busy and updates heterogeneous tracking. Kotlin line 12170. */
        private void markServerBusy(int queueIdx, int serverId, int typeId) {
            serverBusy[queueIdx][serverId] = true;
            if (numServerTypes[queueIdx] > 0) {
                int effectiveTypeId = (typeId >= 0) ? typeId : serverToType[queueIdx][serverId];
                if (effectiveTypeId >= 0
                        && effectiveTypeId < busyCountPerType[queueIdx].length) {
                    busyCountPerType[queueIdx][effectiveTypeId]++;
                }
            }
        }

        /** Marks a server as idle and updates heterogeneous tracking. Kotlin line 12183. */
        private void markServerIdle(int queueIdx, int serverId) {
            serverBusy[queueIdx][serverId] = false;
            if (numServerTypes[queueIdx] > 0
                    && serverId < serverToType[queueIdx].length) {
                int typeId = serverToType[queueIdx][serverId];
                if (typeId >= 0 && typeId < busyCountPerType[queueIdx].length) {
                    busyCountPerType[queueIdx][typeId]--;
                }
            }
        }

        /**
         * Generates service time for heterogeneous servers.
         * Uses server-type-specific distribution if available.
         * Falls back to homogeneous generation if not heterogeneous or typeId is -1.
         * Kotlin line 12338.
         */
        private double generateHeteroServiceTime(int queueIdx, int classId, int typeId) {
            // If not heterogeneous or typeId is -1, use homogeneous generation
            if (numServerTypes[queueIdx] == 0 || typeId < 0) {
                return generateServiceTime(queueIdx, classId);
            }

            ProcessType procType = heteroServiceProcType[queueIdx][typeId][classId];

            // Check for disabled service
            if (procType == ProcessType.DISABLED) {
                // Fall back to homogeneous (which will throw if also disabled)
                return generateServiceTime(queueIdx, classId);
            }

            double baseServiceTime = 0.0;

            // Try to use the heterogeneous generator
            RandomVariateGen gen = heteroServiceGens[queueIdx][typeId][classId];
            if (gen != null) {
                baseServiceTime = gen.nextDouble();
            } else if (procType == ProcessType.MMAP) {
                // Use mmap_sample for MMAP distributions
                MatrixCell proc = heteroServiceProc[queueIdx][typeId][classId];
                Random rng = heteroServiceRng[queueIdx][typeId][classId];
                if (proc != null && rng != null && proc.size() >= 2) {
                    jline.io.Ret.mamMMAPSample mmapResult =
                            jline.api.mam.Mmap_sample.mmap_sample(proc, 1L, rng);
                    baseServiceTime = mmapResult.getSamples()[0];
                }
            } else if (procType == ProcessType.MAP || procType == ProcessType.MMPP2) {
                // Correlated heterogeneous service: carry the modulating phase across services.
                MatrixCell proc = heteroServiceProc[queueIdx][typeId][classId];
                Random rng = heteroServiceRng[queueIdx][typeId][classId];
                if (proc != null && rng != null && proc.size() >= 2) {
                    Map_sample.MapSampler s = heteroServiceMapSampler[queueIdx][typeId][classId];
                    if (s == null) {
                        s = new Map_sample.MapSampler(proc.get(0), proc.get(1));
                        heteroServiceMapSampler[queueIdx][typeId][classId] = s;
                    }
                    baseServiceTime = s.next(rng);
                }
            } else if (procType == ProcessType.RAP) {
                // Correlated heterogeneous rational service: propagate the conditional vector.
                MatrixCell proc = heteroServiceProc[queueIdx][typeId][classId];
                Random rng = heteroServiceRng[queueIdx][typeId][classId];
                if (proc != null && rng != null && proc.size() >= 2) {
                    Rap_sample.RapSampler s = heteroServiceRapSampler[queueIdx][typeId][classId];
                    if (s == null) {
                        s = new Rap_sample.RapSampler(proc.get(0), proc.get(1));
                        heteroServiceRapSampler[queueIdx][typeId][classId] = s;
                    }
                    baseServiceTime = s.next(rng);
                }
            } else if (procType == ProcessType.ME) {
                // Matrix-exponential heterogeneous service: cached inverse-CDF sampling.
                MatrixCell proc = heteroServiceProc[queueIdx][typeId][classId];
                Random rng = heteroServiceRng[queueIdx][typeId][classId];
                if (proc != null && rng != null && proc.size() >= 2) {
                    Me_sample.MeSampler s = heteroServiceMeSampler[queueIdx][typeId][classId];
                    if (s == null) {
                        s = new Me_sample.MeSampler(meAlphaOf(proc, "service"), proc.get(0));
                        heteroServiceMeSampler[queueIdx][typeId][classId] = s;
                    }
                    baseServiceTime = s.next(rng);
                }
            } else if (procType == ProcessType.PH || procType == ProcessType.APH
                    || procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN
                    || procType == ProcessType.COX2) {
                // Use map_sample for renewal phase-type distributions
                MatrixCell proc = heteroServiceProc[queueIdx][typeId][classId];
                Random rng = heteroServiceRng[queueIdx][typeId][classId];
                if (proc != null && rng != null && proc.size() >= 2) {
                    Matrix D0 = proc.get(0);
                    Matrix D1 = proc.get(1);
                    double[] samples = jline.api.mam.Map_sample.map_sample(D0, D1, 1L, rng);
                    baseServiceTime = samples[0];
                }
            } else {
                // Fallback to rate-based exponential
                double rate = heteroMus[queueIdx][typeId][classId];
                if (rate > 0 && rate < Double.MAX_VALUE) {
                    baseServiceTime = -Math.log(routingRng.nextDouble()) / rate;
                }
            }

            // Apply load-dependent scaling if this station has load dependence
            if (isLoadDependent[queueIdx] && !isPSScheduling(schedStrategies[queueIdx])) {
                int totalJobs = getTotalCustomersAtStation(queueIdx);
                if (totalJobs > 0) {
                    double[] scalingArray = lldScaling[queueIdx];
                    if (scalingArray != null) {
                        int scalingIdx = Math.min(totalJobs - 1, scalingArray.length - 1);
                        double scalingFactor = scalingArray[scalingIdx];
                        if (scalingFactor > 0) {
                            baseServiceTime /= scalingFactor;
                        }
                    }
                }
            }

            baseServiceTime = resolveZeroAtom(baseServiceTime, procType,
                    "service time at station " + queueIdx + " class " + classId);
            baseServiceTime = slotSnap(baseServiceTime, "service time at station "
                    + queueIdx + " class " + classId);

            // Track sample for control variates if enabled and after warmup
            if (useControlVariates && warmupDone && baseServiceTime > 0) {
                serviceSampleSum[queueIdx][classId] += baseServiceTime;
                serviceSampleCount[queueIdx][classId]++;
            }

            return baseServiceTime;
        }

        /**
         * Initiates server setup (cold start) phase.
         * The server transitions to SETUP state and schedules a SetupCompletion event.
         * Kotlin line 12509.
         */
        private void startServerSetup(int queueIdx, int serverId, int classId) {
            // Update server state to SETUP
            serverState[queueIdx][serverId] = ServerState.SETUP;

            // Track which class triggered this setup
            serverLastClass[queueIdx][serverId] = classId;

            // Generate setup time from distribution
            double setupTime = generateSetupTime(queueIdx, classId);

            // Update setup statistics
            updateSetupStats(queueIdx, classId);
            currentServersInSetup[queueIdx][classId]++;
            lastSetupUpdateTime[queueIdx][classId] = ssjSim.time();

            // Schedule setup completion event
            SetupCompletion setupEvent = new SetupCompletion(queueIdx, serverId, classId);
            setupEvent.schedule(setupTime);

            // Log setup start
            int stationIdx = serviceStations.get(queueIdx);
            logEvent("SETUP_START", stationIdx, classId, serverId, 0);
        }

        /**
         * Cancels an ongoing delayoff phase and returns the server to ACTIVE state.
         * This happens when a job arrives while the server is in the delayoff phase.
         * Kotlin line 12573.
         */
        private void cancelDelayoff(int queueIdx, int serverId) {
            // Cancel the pending delayoff event
            Event event = pendingDelayoffEvents[queueIdx][serverId];
            if (event != null) {
                event.cancel();
                pendingDelayoffEvents[queueIdx][serverId] = null;
            }

            // Get the class that was in delayoff
            int classId = serverLastClass[queueIdx][serverId];
            if (classId >= 0) {
                // Update delayoff statistics
                updateDelayoffStats(queueIdx, classId);
                currentServersInDelayoff[queueIdx][classId]--;
            }

            // Update server state to ACTIVE (resume without setup)
            serverState[queueIdx][serverId] = ServerState.ACTIVE;

            // Log delayoff cancellation
            int stationIdx = serviceStations.get(queueIdx);
            logEvent("DELAYOFF_CANCEL", stationIdx, classId, serverId, 0);
        }

        /**
         * Starts service for a customer on a specific server.
         * Updates statistics and schedules a departure event.
         * Kotlin line 12604.
         */
        private void startService(int queueIdx, int serverId, Customer customer, int serverTypeId) {
            int classId = customer.classId;

            // Cancel any scheduled reneging event (customer is starting service, no longer waiting)
            cancelRenegingIfScheduled(queueIdx, customer);

            // Mark server as busy (with heterogeneous tracking)
            markServerBusy(queueIdx, serverId, serverTypeId);
            customersInService[queueIdx]++;

            // Store assigned server type on customer
            customer.assignedServerType = serverTypeId;

            // Update busy time statistics
            updateBusyStats(queueIdx, classId);
            currentBusyServers[queueIdx][classId]++;
            lastBusyUpdateTime[queueIdx][classId] = ssjSim.time();

            // Generate or use existing service time (heterogeneous-aware)
            double serviceTime;
            if (customer.serviceTime > 0) {
                serviceTime = customer.serviceTime;
            } else {
                serviceTime = generateHeteroServiceTime(queueIdx, classId, serverTypeId);
            }

            // Schedule departure event
            Departure departureEvent = new Departure(queueIdx, serverId, customer);
            departureEvent.schedule(serviceWallDelay(queueIdx, classId, serviceTime));

            // Track class-dependent departure for dynamic re-scaling
            if ((hasCd || hasLld) && sdDepartureEvents != null) {
                sdDepartureEvents[queueIdx].put(serverId, departureEvent);
                sdInServiceCustomers[queueIdx].put(serverId, customer);
            }

            // Track for signal-based removal (G-networks)
            if (hasRemovalSignals) {
                inServiceJobs.put(new IntPair(queueIdx, serverId),
                        new InServiceJob(customer, departureEvent));
            }
        }

        /** Departure event for a Queue (FCFS-family). Kotlin line 10378.
         *  Body translated in PART 8 (see {@link #departureActions(int, int, Customer)}). */
        private final class Departure extends SimEvent {
            final int queueIdx;
            final int serverId;
            final Customer customer;

            Departure(int queueIdx, int serverId, Customer customer) {
                this.queueIdx = queueIdx;
                this.serverId = serverId;
                this.customer = customer;
            }

            @Override
            public void actions() {
                departureActions(queueIdx, serverId, customer);
            }
        }

        /**
         * Firing of a Batch Markovian Service Process (BMSP) at a bulk-service
         * station. Serves min(k, N) jobs FCFS, where k is the batch size drawn
         * with the inter-firing time and N is the number present at the firing
         * instant (partial-batch truncation). If jobs remain the clock is
         * re-armed; otherwise the server idles until the next arrival.
         */
        private final class BatchServiceDeparture extends SimEvent {
            final int queueIdx;
            final int classId;

            BatchServiceDeparture(int queueIdx, int classId) {
                this.queueIdx = queueIdx;
                this.classId = classId;
            }

            @Override
            public void actions() {
                batchServiceDepartureActions(queueIdx, classId);
            }
        }

        /**
         * Enqueues an arrival at a bulk-service station and starts the BMSP
         * clock if the station was previously empty. The single server is busy
         * for the whole non-empty period, so utilization equals the fraction of
         * time the station is non-empty.
         */
        private void handleBatchServiceArrival(int queueIdx, Customer customer) {
            int classId = customer.classId;
            waitQueues[queueIdx].add(customer);
            // currentQueueLength[queueIdx][classId] was already incremented by
            // arriveAtQueue before this call; == number present for a
            // single-class station.
            if (currentQueueLength[queueIdx][classId] == 1) {
                markServerBusy(queueIdx, 0, -1);
                customersInService[queueIdx] = 1;
                updateBusyStats(queueIdx, classId);
                currentBusyServers[queueIdx][classId] = 1;
                lastBusyUpdateTime[queueIdx][classId] = ssjSim.time();
                scheduleNextBatchService(queueIdx, classId);
            }
        }

        /** Draws the next (inter-firing time, batch size) pair from the station's
         *  BMSP clock, retaining the modulating phase, and schedules the firing. */
        private void scheduleNextBatchService(int queueIdx, int classId) {
            jline.api.mam.BmapSample sample =
                    batchServiceSampler[queueIdx].next(batchServiceRng[queueIdx]);
            batchPendingSize[queueIdx] = sample.getBatchSize();
            new BatchServiceDeparture(queueIdx, classId)
                    .schedule(sample.getInterarrivalTime());
        }

        private void batchServiceDepartureActions(int queueIdx, int classId) {
            trackEvent();
            int k = batchPendingSize[queueIdx];
            int present = currentQueueLength[queueIdx][classId];
            int d = Math.min(k, present);

            // Phase A: remove the served jobs FCFS and record their completions,
            // so that the remaining count is settled before any routing (a job
            // routed back to this same station would otherwise re-arm the clock).
            List<Customer> served = new ArrayList<Customer>(d);
            for (int i = 0; i < d; i++) {
                Customer cust = waitQueues[queueIdx].poll();
                if (cust == null) {
                    break;
                }
                recordBatchMemberCompletion(queueIdx, cust);
                served.add(cust);
            }

            // Re-arm the clock while busy, or idle the server when empty.
            if (currentQueueLength[queueIdx][classId] >= 1) {
                scheduleNextBatchService(queueIdx, classId);
            } else {
                updateBusyStats(queueIdx, classId);
                currentBusyServers[queueIdx][classId] = 0;
                markServerIdle(queueIdx, 0);
                customersInService[queueIdx] = 0;
            }

            // Phase B: route each served job to its destination.
            for (int i = 0; i < served.size(); i++) {
                routeBatchMember(queueIdx, served.get(i));
            }
        }

        /** Records a single bulk-service completion (response time, throughput,
         *  queue-length decrement, region leave). Mirrors the accounting head of
         *  {@link #departureActions(int, int, Customer)} without per-job server
         *  bookkeeping (the bulk server is tracked at station granularity). */
        private void recordBatchMemberCompletion(int queueIdx, Customer customer) {
            int classId = customer.classId;
            double queueResponseTime = ssjSim.time() - customer.queueArrivalTime;
            responseTimeTally[queueIdx][classId].add(queueResponseTime);
            responseTimeSamples[queueIdx][classId].add(queueResponseTime);
            completedCustomers[queueIdx][classId]++;

            double queueTardiness = Math.max(0.0, ssjSim.time() - customer.absoluteDeadline);
            tardinessTally[queueIdx][classId].add(queueTardiness);

            checkEventCountStop();

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
            }

            logEvent("DEPARTURE", stationIdx, classId,
                    currentQueueLength[queueIdx][classId],
                    currentBusyServers[queueIdx][classId]);
        }

        /** Routes a completed bulk-service job to its destination (sink or the
         *  next station). Fork/join and signal destinations are unsupported for
         *  bulk service and rejected explicitly rather than mishandled. */
        private void routeBatchMember(int queueIdx, Customer customer) {
            int classId = customer.classId;
            int currentNode = serviceNodes.get(queueIdx);
            RoutingResult routingResult = selectDestination(currentNode, classId);
            int destNode = routingResult.destNode;
            int destClassId = routingResult.destClassId;

            regionExitCompleted(queueIdx, classId, destNode, destClassId);

            if (destNode >= 0 && sinkNodes.contains(Integer.valueOf(destNode))) {
                double systemResponseTime = ssjSim.time() - customer.systemArrivalTime;
                systemResponseTimeTally[classId].add(systemResponseTime);
                double systemTardiness = Math.max(0.0, ssjSim.time() - customer.absoluteDeadline);
                systemTardinessTally[classId].add(systemTardiness);
                systemCompletedCustomers[classId]++;
                return;
            }
            if (destNode >= 0 && (forkNodes.contains(Integer.valueOf(destNode))
                    || joinNodes.contains(Integer.valueOf(destNode)))) {
                throw new RuntimeException("LDES: batch (BMSP) service routing to a "
                        + "Fork/Join node is not supported.");
            }
            if (destNode >= 0
                    && ((hasNegativeSignals && isNegativeSignal[destClassId])
                        || (hasReplySignals && isReplySignal[destClassId])
                        || (hasCatastropheSignals && isCatastropheSignal[destClassId]))) {
                throw new RuntimeException("LDES: batch (BMSP) service routing to a "
                        + "signal class is not supported.");
            }
            if (destNode >= 0) {
                int nextQueueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
                if (nextQueueIdx >= 0) {
                    Customer nextCustomer = new Customer(
                            destClassId, classPrio[destClassId],
                            customer.systemArrivalTime, ssjSim.time(),
                            siroRng.nextDouble(), -1.0, customer.jobId,
                            customer.absoluteDeadline, -1, null);
                    arriveAtQueue(nextQueueIdx, nextCustomer);
                }
            }
        }

        /** DelayDeparture event for a Delay (infinite-server) node. Kotlin line 10958.
         *  Body translated in PART 8 (see {@link #delayDepartureActions(int, Customer)}). */
        private final class DelayDeparture extends SimEvent {
            final int queueIdx;
            final Customer customer;

            DelayDeparture(int queueIdx, Customer customer) {
                this.queueIdx = queueIdx;
                this.customer = customer;
            }

            @Override
            public void actions() {
                delayDepartureActions(queueIdx, customer);
            }
        }

    // ==================== PART 7 TRANSLATION (Kotlin lines 9001-10500) ====================

        /**
         * Forwards stub: helpers from later chunks.
         */
        /** Forward stub — Kotlin line 11457 (handleForkArrival already declared above). */
        /** Forward stub — Kotlin line 11606 (arriveAtQueueForked already declared above). */
        /** Forward stub — Kotlin line 11631 (handleJoinArrival already declared above). */
        /** Forward stub — Kotlin line 11904 (handleUnknownJoinArrival already declared above). */

        /** Map key (queueIdx, queueArrivalTime, randomRank) — used by
         *  forkedCustomerMap. The rank is what tells sibling tasks apart:
         *  a Fork routes all of them at the same instant, so two siblings
         *  entering the same station share (queueIdx, queueArrivalTime) and
         *  one would overwrite the other, losing its forked identity and
         *  stalling the Join forever.
         *  Forward declaration; the map itself is wired up in a later chunk. */
        private static final class IntDoubleKey {
            final int a;
            final double b;
            final double c;
            IntDoubleKey(int a, double b, double c) { this.a = a; this.b = b; this.c = c; }
            @Override public boolean equals(Object o) {
                if (this == o) return true;
                if (!(o instanceof IntDoubleKey)) return false;
                IntDoubleKey k = (IntDoubleKey) o;
                return a == k.a
                        && Double.doubleToLongBits(b) == Double.doubleToLongBits(k.b)
                        && Double.doubleToLongBits(c) == Double.doubleToLongBits(k.c);
            }
            @Override public int hashCode() {
                long bb = Double.doubleToLongBits(b);
                long cc = Double.doubleToLongBits(c);
                int h = 31 * a + (int) (bb ^ (bb >>> 32));
                return 31 * h + (int) (cc ^ (cc >>> 32));
            }
        }
        /** Map: (queueIdx, queueArrivalTime) -> ForkedJob.
         *  Populated by arriveAtQueueForked (translated in later chunk). */
        private final Map<IntDoubleKey, ForkedJob> forkedCustomerMap =
                new HashMap<IntDoubleKey, ForkedJob>();

        /** Body translated in PART 8 (see {@link #findFreeServerImpl(int)}). */
        private int findFreeServer(int queueIdx) {
            return findFreeServerImpl(queueIdx);
        }

        /**
         * Preempts a currently serving job, saving its remaining work and adding
         * it back to the wait queue.
         */
        private void preemptJob(int queueIdx, PreemptiveCustomer job) {
            // Cancel the scheduled departure
            if (job.scheduledDepartureEvent != null) {
                job.scheduledDepartureEvent.cancel();
            }

            // Compute remaining work and elapsed time in operational time, so a
            // job preempted under a time-varying rate resumes with the work it
            // actually accrued (int mu over the served interval), not the wall
            // duration. serviceWorkBetween is the wall elapsed off a schedule.
            double currentTime = ssjSim.time();
            double elapsed = serviceWorkBetween(queueIdx, job.classId,
                    job.serviceStartTime, currentTime);
            double newRemainingWork = job.remainingServiceWork - elapsed;

            // Update total elapsed service time on the job (for FB scheduling)
            job.elapsedServiceTime += elapsed;

            // Remove from in-service list
            preemptiveJobsInService[queueIdx].remove(job);

            // For LCFSPR (preemptive resume) and FCFSPR/FCFSPRPRIO (size-based preemptive),
            // save state for accurate resumption
            SchedStrategy strategy = schedStrategies[queueIdx];
            if (strategy == SchedStrategy.LCFSPR || strategy == SchedStrategy.LCFSPRPRIO ||
                    strategy == SchedStrategy.FCFSPR || strategy == SchedStrategy.FCFSPRPRIO ||
                    isPSJFScheduling(strategy) || isFBScheduling(strategy) || isLRPTScheduling(strategy) ||
                    isSRPTScheduling(strategy) || isFSPScheduling(strategy)) {
                PreemptionKey customerKey = new PreemptionKey(queueIdx,
                        job.systemArrivalTime, job.queueArrivalTime);
                // Capture distribution type and phase parameter for PH distribution tracking
                ProcessType distType = serviceProcessType[queueIdx][job.classId];
                Integer phaseParam = null;

                // For Erlang distributions, extract the shape parameter (number of phases)
                if (distType == ProcessType.ERLANG) {
                    RandomVariateGen gen = serviceGens[queueIdx][job.classId];
                    if (gen != null) {
                        try {
                            // Attempt to get the shape parameter from ErlangGen using reflection
                            // ErlangGen has a 'k' field representing the number of phases
                            java.lang.reflect.Field kField = gen.getClass().getDeclaredField("k");
                            kField.setAccessible(true);
                            Object val = kField.get(gen);
                            phaseParam = (val instanceof Number) ? ((Number) val).intValue() : null;
                        } catch (Exception e) {
                            // If we can't access the field, just use null
                            phaseParam = null;
                        }
                    }
                }

                // Save residual work, original total (for accurate statistics), elapsed time, and distribution tracking info
                preemptedJobHistory.put(customerKey, new PreemptionRecord(
                        newRemainingWork,
                        job.totalServiceRequirement,
                        job.elapsedServiceTime,
                        distType,
                        phaseParam));
            }

            // Add back to wait queue as a regular customer with remaining work preserved
            Customer preemptedCustomer = new Customer(
                    job.classId,
                    job.priority,
                    job.systemArrivalTime,
                    job.queueArrivalTime,
                    job.randomRank,
                    newRemainingWork  // Preserve remaining work for SRPT comparisons
            );
            waitQueues[queueIdx].add(preemptedCustomer);

            // Update statistics (job leaves service but stays in queue)
            updateBusyStats(queueIdx, job.classId);
            currentBusyServers[queueIdx][job.classId]--;
            customersInService[queueIdx]--;

            // DO NOT update queue length - job stays in system

            logEvent("PREEMPT", serviceStations.get(queueIdx), job.classId,
                    currentQueueLength[queueIdx][job.classId],
                    currentBusyServers[queueIdx][job.classId]);
        }

        /**
         * Starts preemptive service for a customer.
         * For new arrivals, generates fresh service time.
         * For resumed jobs, uses saved remaining work (LCFSPR) or resamples (LCFSPI).
         */
        private void startPreemptiveService(int queueIdx, Customer customer, int serverId) {
            startPreemptiveService(queueIdx, customer, serverId, -1);
        }

        private void startPreemptiveService(int queueIdx, Customer customer, int serverId, int serverTypeId) {
            int classId = customer.classId;
            double currentTime = ssjSim.time();
            SchedStrategy strategy = schedStrategies[queueIdx];

            // Mark server busy (with heterogeneous tracking)
            markServerBusy(queueIdx, serverId, serverTypeId);
            customersInService[queueIdx]++;

            // Store assigned server type on customer
            customer.assignedServerType = serverTypeId;

            // Update busy statistics
            updateBusyStats(queueIdx, classId);
            currentBusyServers[queueIdx][classId]++;

            // Determine service time based on whether this is a resumed job
            // Use stable key based on arrival times (matches key used in preemptJob)
            PreemptionKey customerKey = new PreemptionKey(queueIdx,
                    customer.systemArrivalTime, customer.queueArrivalTime);
            PreemptionRecord savedRecord = preemptedJobHistory.remove(customerKey);

            double serviceTime;
            double totalRequirement;
            double remainingWork;
            double elapsedTime;

            if (savedRecord != null && (strategy == SchedStrategy.LCFSPR
                    || strategy == SchedStrategy.LCFSPRPRIO
                    || strategy == SchedStrategy.FCFSPR
                    || strategy == SchedStrategy.FCFSPRPRIO
                    || isSizeBasedPreemptiveScheduling(strategy))) {
                // Resume with saved state (LCFSPR, FCFSPR, SRPT, PSJF, FB, LRPT) - use residual remaining work
                // All size-based preemptive policies and PR variants are preemptive-resume
                // Distribution type and phase parameter are tracked for PH distributions (no correction applied)
                serviceTime = savedRecord.remainingWork;
                totalRequirement = savedRecord.originalTotal;  // FIXED: use original total for statistics
                remainingWork = savedRecord.remainingWork;
                elapsedTime = savedRecord.elapsedTime;
            } else {
                // Fresh service time (new arrival or LCFSPI resample) - heterogeneous-aware
                double freshServiceTime;
                if (isSizeBasedPreemptiveScheduling(strategy)
                        && (numServerTypes[queueIdx] == 0 || serverTypeId < 0)
                        && customer.serviceTime > 0.0) {
                    // Size-based preemptive policies (SRPT, PSJF, FB, LRPT, FSP, EDF)
                    // sample the job size once at arrival and use it both for the
                    // preemption decision and for wait-queue ordering. The realized
                    // service requirement must be that same size: re-sampling a fresh
                    // draw here would decouple the served work from the scheduling
                    // priority, degrading SRPT toward random-order preemption and
                    // inflating response times. Homogeneous servers only; with
                    // heterogeneous servers the service time depends on the server
                    // type (unknown at arrival), so it is redrawn below.
                    freshServiceTime = customer.serviceTime;
                } else {
                    freshServiceTime = generateHeteroServiceTime(queueIdx, classId, serverTypeId);
                }
                serviceTime = freshServiceTime;
                totalRequirement = freshServiceTime;
                remainingWork = freshServiceTime;
                elapsedTime = 0.0;
            }

            // Create preemptive customer record
            PreemptiveCustomer preemptiveCustomer = new PreemptiveCustomer(
                    classId,
                    customer.priority,
                    customer.systemArrivalTime,
                    customer.queueArrivalTime,
                    customer.randomRank,
                    totalRequirement,
                    remainingWork,
                    elapsedTime,  // Preserve elapsed time from before preemption
                    currentTime,
                    null,
                    serverId,
                    serverTypeId,
                    customer.absoluteDeadline);

            // Schedule departure
            PreemptiveDeparture departureEvent = new PreemptiveDeparture(queueIdx, serverId, preemptiveCustomer);
            departureEvent.schedule(serviceWallDelay(queueIdx, classId, serviceTime));
            preemptiveCustomer.scheduledDepartureEvent = departureEvent;

            // Track in service
            preemptiveJobsInService[queueIdx].add(preemptiveCustomer);
        }

        /** PreemptiveDeparture event for preemptive LCFS family. Kotlin line 10751.
         *  Body translated in PART 8 (see {@link #preemptiveDepartureActions(int, int, PreemptiveCustomer)}). */
        private final class PreemptiveDeparture extends SimEvent {
            final int queueIdx;
            final int serverId;
            final PreemptiveCustomer customer;

            PreemptiveDeparture(int queueIdx, int serverId, PreemptiveCustomer customer) {
                this.queueIdx = queueIdx;
                this.serverId = serverId;
                this.customer = customer;
            }

            @Override
            public void actions() {
                preemptiveDepartureActions(queueIdx, serverId, customer);
            }
        }

        // =========================================================================
        // POLLING SCHEDULING IMPLEMENTATION
        // =========================================================================

        /**
         * Customer arrives at a polling station.
         * Jobs are added to per-class queues and served according to polling discipline.
         */
        private boolean arriveAtPollingQueue(int queueIdx, Customer customer) {
            int classId = customer.classId;

            // Update queue length statistics
            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]++;

            // Track max queue length reached
            int totalAtStation = getTotalCustomersAtStation(queueIdx);
            if (totalAtStation > maxQueueLengthReached) {
                maxQueueLengthReached = totalAtStation;
            }

            // Update region job counts if in a region
            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobEnter(regionIdx, classId);
                updateRegionArrivalTracking(regionIdx, classId);
            }

            // Add to per-class polling queue
            pollingQueues[queueIdx][classId].add(customer);

            // Log arrival event
            logEvent("ARRIVAL", serviceStations.get(queueIdx), classId,
                    currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId]);

            // If server is idle and not in switchover, try to start service
            if (customersInService[queueIdx] == 0 && !pollingInSwitchover[queueIdx]) {
                // Check if the current class has jobs or find next class with jobs
                pollingTryStartService(queueIdx);
            }

            return true;
        }

        /**
         * Try to start service at a polling station.
         * Called when server becomes idle or after switchover completes.
         */
        private void pollingTryStartService(int queueIdx) {
            PollingType type = pollingType[queueIdx];
            if (type == null) {
                type = PollingType.EXHAUSTIVE;
            }

            // For GATED and DECREMENTING: at the start of a visit to a new class, record the
            // number of waiting jobs found at the polling instant. GATED serves exactly this
            // many; DECREMENTING serves until the count drops to one less than this value.
            if ((type == PollingType.GATED || type == PollingType.DECREMENTING)
                    && pollingJobsServedInRound[queueIdx] == 0) {
                int currentClass0 = pollingCurrentClass[queueIdx];
                pollingGateSize[queueIdx] = pollingQueues[queueIdx][currentClass0].size();
            }

            // Try to get a job from the current class
            int currentClass = pollingCurrentClass[queueIdx];
            LinkedList<Customer> queue = pollingQueues[queueIdx][currentClass];

            if (pollingCanServeFromCurrentClass(queueIdx, type)) {
                // Serve the next job from current class
                Customer job = queue.poll();
                if (job != null) {
                    pollingStartService(queueIdx, job);
                }
            } else {
                // Need to switch to next class
                pollingInitiateSwitchover(queueIdx);
            }
        }

        /**
         * Check if we can serve another job from the current class based on polling type.
         */
        private boolean pollingCanServeFromCurrentClass(int queueIdx, PollingType type) {
            int currentClass = pollingCurrentClass[queueIdx];
            int queueSize = pollingQueues[queueIdx][currentClass].size();

            switch (type) {
                case EXHAUSTIVE:
                    return queueSize > 0;
                case GATED: {
                    int gateSize = pollingGateSize[queueIdx];
                    int served = pollingJobsServedInRound[queueIdx];
                    return queueSize > 0 && served < gateSize;
                }
                case KLIMITED: {
                    int k = pollingK[queueIdx];
                    int served = pollingJobsServedInRound[queueIdx];
                    return queueSize > 0 && served < k;
                }
                case DECREMENTING: {
                    // Semiexhaustive: continue serving while the number of waiting jobs is
                    // still at least the count found at the polling instant; stop once it
                    // has decreased to one less than that count (net decrease of one).
                    int startCount = pollingGateSize[queueIdx];
                    return queueSize > 0 && queueSize >= startCount;
                }
                default:
                    return false;
            }
        }

        /**
         * Start service for a job at a polling station.
         */
        private void pollingStartService(int queueIdx, Customer customer) {
            int classId = customer.classId;

            // Find compatible server for heterogeneous-aware polling (if applicable)
            ServerSelection serverSelection = findFreeServerForClass(queueIdx, classId);
            int serverId = (serverSelection.serverId >= 0) ? serverSelection.serverId : 0;
            int serverTypeId = serverSelection.serverTypeId;

            // Mark server busy with heterogeneous tracking
            markServerBusy(queueIdx, serverId, serverTypeId);
            customersInService[queueIdx]++;
            customer.assignedServerType = serverTypeId;

            // Flush accumulated busy time for classId before updating tracker
            updateBusyStats(queueIdx, classId);
            currentBusyServers[queueIdx][classId]++;

            // Generate and schedule service time (heterogeneous-aware)
            double serviceTime = generateHeteroServiceTime(queueIdx, classId, serverTypeId);
            new PollingDeparture(queueIdx, customer).schedule(serviceWallDelay(queueIdx, classId, serviceTime));
        }

        /**
         * Initiate switchover to the next class at a polling station.
         */
        private void pollingInitiateSwitchover(int queueIdx) {
            // Classical cyclic polling: the server advances to the NEXT buffer in the
            // cyclic order, one step at a time, and pays that step's switchover whether
            // or not the buffer it reaches turns out to hold work. Walking over an empty
            // buffer therefore costs a switchover, and skipping ahead to the next
            // non-empty buffer for a single switchover would be a different (and much
            // faster) discipline: it shortens the cycle and understates the waiting time
            // by tens of percent against the exact solution. This is the model assumed by
            // Takagi's formulas (jline.api.polling, MATLAB matlab/src/api/polling) and by
            // SolverCTMC, both of which this now agrees with.
            //
            // A server with an empty station keeps lapping too, which is what makes the
            // cycle time r/(1-rho) tend to r as rho -> 0. The single exception is a lap
            // made entirely of zero-time legs, which would spin forever at the same
            // simulated instant and is unobservable anyway: there the server parks and
            // is restarted by the next arrival (see arriveAtPollingQueue).
            for (int step = 0; step < numClasses; step++) {
                int fromClass = pollingCurrentClass[queueIdx];
                int nextClass = (fromClass + 1) % numClasses;

                // Queue.setSwitchover(class_i, distrib) is documented as the time to
                // switch FROM buffer i to the next one, which is also Takagi's r_i, so
                // the leg out of fromClass is charged the switchover stored against
                // fromClass -- not the one against the buffer being entered.
                RandomVariateGen switchoverGen = pollingSwitchoverGens[queueIdx][fromClass];
                double switchoverTime = (switchoverGen != null) ? switchoverGen.nextDouble() : 0.0;

                if (switchoverTime > 0) {
                    pollingInSwitchover[queueIdx] = true;
                    new PollingSwitchover(queueIdx, nextClass).schedule(switchoverTime);
                    return;
                }

                // Zero-time leg: the server arrives at nextClass at this same instant.
                pollingCurrentClass[queueIdx] = nextClass;
                pollingJobsServedInRound[queueIdx] = 0;
                if (!pollingQueues[queueIdx][nextClass].isEmpty()) {
                    pollingTryStartService(queueIdx);
                    return;
                }
            }
            // A full lap of zero-time legs found no work anywhere: park.
            pollingJobsServedInRound[queueIdx] = 0;
        }

        /**
         * Departure event for polling station.
         */
        private final class PollingDeparture extends SimEvent {
            private final int queueIdx;
            private final Customer customer;

            PollingDeparture(int queueIdx, Customer customer) {
                this.queueIdx = queueIdx;
                this.customer = customer;
            }

            @Override
            public void actions() {
                trackEvent();
                int classId = customer.classId;

                // Record queue response time
                double queueResponseTime = ssjSim.time() - customer.queueArrivalTime;
                responseTimeTally[queueIdx][classId].add(queueResponseTime);
                responseTimeSamples[queueIdx][classId].add(queueResponseTime);
                completedCustomers[queueIdx][classId]++;

                // Check event count for stopping/warmup/MSER sampling
                checkEventCountStop();

                // Update queue length statistics
                updateQueueStats(queueIdx, classId);
                currentQueueLength[queueIdx][classId]--;

                // Update region job counts if in a region and try releasing blocked customers
                int stationIdx = serviceStations.get(queueIdx);
                if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                    int regionIdx = fcRegionIndices.get(stationIdx);
                    updateRegionTimeWeightedStats(regionIdx);
                    regionJobLeave(regionIdx, classId);
                    // completion count and blocked-FIFO release are deferred to
                    // regionExitCompleted once the routing destination is known:
                    // intra-region hops must not release into the transient slot
                }

                // LQN phase-2: the spawn continuation takes over the freed slot
                maybeSpawnOnCompletion(queueIdx, classId, customer);

                // Update busy time before decrementing
                updateBusyStats(queueIdx, classId);
                currentBusyServers[queueIdx][classId]--;

                // Free server (with heterogeneous tracking)
                markServerIdle(queueIdx, 0);
                customersInService[queueIdx]--;

                // Increment jobs served in this round
                pollingJobsServedInRound[queueIdx]++;

                // Log departure event
                logEvent("DEPARTURE", stationIdx, classId,
                        currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId]);

                // Route customer to next destination
                int currentNode = serviceNodes.get(queueIdx).intValue();
                RoutingResult routingResult = selectDestination(currentNode, classId);
                int destNode = routingResult.destNode;
                int destClassId = routingResult.destClassId;
                // FCR: on a true region exit, count the completion and release
                // blocked customers; a no-op for intra-region hops
                regionExitCompleted(queueIdx, classId, destNode, destClassId);

                if (destNode >= 0 && sinkNodes.contains(Integer.valueOf(destNode))) {
                    // Customer leaves system
                    double systemResponseTime = ssjSim.time() - customer.systemArrivalTime;
                    systemResponseTimeTally[classId].add(systemResponseTime);
                    systemCompletedCustomers[classId]++;
                } else if (destNode >= 0 && forkNodes.contains(Integer.valueOf(destNode))) {
                    // Destination is a Fork node
                    long parentJobId = nextJobId++;
                    handleForkArrival(destNode, parentJobId, destClassId, customer.systemArrivalTime);
                } else if (destNode >= 0 && joinNodes.contains(Integer.valueOf(destNode))) {
                    // Destination is a Join node
                    IntDoubleKey forkedJobKey = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
                    ForkedJob forkedJob = forkedCustomerMap.remove(forkedJobKey);
                    if (forkedJob != null) {
                        ForkedJob updatedForkedJob = new ForkedJob(
                                forkedJob.forkJobId,
                                forkedJob.parentJobId,
                                destClassId,
                                forkedJob.priority,
                                forkedJob.systemArrivalTime,
                                ssjSim.time(),
                                forkedJob.randomRank);
                        handleJoinArrival(destNode, updatedForkedJob);
                    } else {
                        handleUnknownJoinArrival(destNode, destClassId, customer.systemArrivalTime);
                    }
                } else if (destNode >= 0) {
                    // Route to another queue (possibly with switched class)
                    int nextQueueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
                    if (nextQueueIdx >= 0) {
                        IntDoubleKey forkedJobKey = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
                        ForkedJob forkedJob = forkedCustomerMap.remove(forkedJobKey);
                        Customer nextCustomer = new Customer(
                                destClassId, classPrio[destClassId],
                                customer.systemArrivalTime, ssjSim.time(),
                                siroRng.nextDouble(),
                                -1.0, -1L, customer.absoluteDeadline, -1, null);
                        if (forkedJob != null) {
                            ForkedJob nextForkedJob = new ForkedJob(
                                    forkedJob.forkJobId,
                                    forkedJob.parentJobId,
                                    destClassId,
                                    forkedJob.priority,
                                    forkedJob.systemArrivalTime,
                                    ssjSim.time(),
                                    forkedJob.randomRank);
                            arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob);
                        } else {
                            arriveAtQueue(nextQueueIdx, nextCustomer);
                        }
                    }
                }

                // Try to serve next job or initiate switchover
                if (!pollingInSwitchover[queueIdx]) {
                    pollingTryStartService(queueIdx);
                }
            }
        }

        /**
         * Switchover event for polling station.
         * Completes transition to the next class.
         */
        private final class PollingSwitchover extends SimEvent {
            private final int queueIdx;
            private final int nextClass;

            PollingSwitchover(int queueIdx, int nextClass) {
                this.queueIdx = queueIdx;
                this.nextClass = nextClass;
            }

            @Override
            public void actions() {
                trackEvent();
                pollingInSwitchover[queueIdx] = false;
                pollingCurrentClass[queueIdx] = nextClass;
                pollingJobsServedInRound[queueIdx] = 0;

                // Try to start service in the new class
                pollingTryStartService(queueIdx);
            }
        }

        // =========================================================================
        // END OF POLLING SCHEDULING IMPLEMENTATION
        // =========================================================================

        /**
         * Returns the total number of customers at a station (waiting + in service).
         * Note: currentQueueLength already includes ALL jobs at the station (both waiting and in service).
         */
        private int getTotalCustomersAtStation(int queueIdx) {
            int total = 0;
            for (int k = 0; k < numClasses; k++) {
                total += currentQueueLength[queueIdx][k];
            }
            return total;
        }

        // ==================== Impatience Helper Functions ====================

        /**
         * Schedules a reneging event for an impatient customer joining a queue.
         * Called when a customer joins a queue with reneging configured.
         */
        private void scheduleRenegingEvent(int queueIdx, Customer customer) {
            int classId = customer.classId;
            if (!hasPatienceConfig[queueIdx][classId]) return;

            RandomVariateGen patienceGen = patienceGens[queueIdx][classId];
            if (patienceGen == null) return;
            double patienceTime = patienceGen.nextDouble();
            if (patienceTime <= 0) return;

            double currentTime = ssjSim.time();
            double deadline = currentTime + patienceTime;

            ImpatientKey key = new ImpatientKey(queueIdx,
                    Double.doubleToLongBits(customer.systemArrivalTime), classId);
            RenegingEvent renegingEvent = new RenegingEvent(queueIdx, customer, key);
            renegingEvent.schedule(patienceTime);

            ImpatientCustomer impatient = new ImpatientCustomer(customer, patienceTime, deadline, renegingEvent);
            waitingImpatientCustomers.put(key, impatient);
        }

        /** RenegingEvent — Kotlin line 11310. Body translated in PART 8
         *  (see {@link #renegingActions(int, Customer, ImpatientKey)}). */
        private final class RenegingEvent extends SimEvent {
            @Override
            protected double slotPhase() {
                return SLOT_PHASE_INTERNAL;
            }

            final int queueIdx;
            final Customer customer;
            final ImpatientKey key;

            RenegingEvent(int queueIdx, Customer customer, ImpatientKey key) {
                this.queueIdx = queueIdx;
                this.customer = customer;
                this.key = key;
            }

            @Override
            public void actions() {
                renegingActions(queueIdx, customer, key);
            }
        }

        /**
         * Cancels a scheduled reneging event when customer starts service.
         */
        private void cancelRenegingIfScheduled(int queueIdx, Customer customer) {
            int classId = customer.classId;
            if (!hasPatienceConfig[queueIdx][classId]) return;

            ImpatientKey key = new ImpatientKey(queueIdx,
                    Double.doubleToLongBits(customer.systemArrivalTime), classId);
            ImpatientCustomer impatient = waitingImpatientCustomers.remove(key);
            if (impatient != null && impatient.renegingEvent != null) {
                impatient.renegingEvent.cancel();
            }
        }

        /**
         * Schedules a retrial event for a customer entering the orbit.
         */
        private void scheduleRetrial(int queueIdx, Customer customer, int attempts, int maxAttempts) {
            int classId = customer.classId;
            RandomVariateGen retrialGen = retrialGens[queueIdx][classId];
            if (retrialGen == null) return;

            double retrialDelay = retrialGen.nextDouble();
            if (retrialDelay <= 0) return;

            OrbitJob orbitJob = new OrbitJob(
                    customer,
                    ssjSim.time(),
                    attempts,
                    maxAttempts,
                    queueIdx);

            RetrialEvent retrialEvent = new RetrialEvent(queueIdx, orbitJob);
            retrialEvent.schedule(retrialDelay);
            orbitJob.retrialEvent = retrialEvent;

            orbitJobs[queueIdx].add(orbitJob);

            // Update orbit statistics
            updateOrbitTimeStats(queueIdx, classId);
            currentOrbitSize[queueIdx][classId]++;
        }

        /** RetrialEvent — Kotlin line 11364. Body translated in PART 8
         *  (see {@link #retrialActions(int, OrbitJob)}). */
        private final class RetrialEvent extends SimEvent {
            @Override
            protected double slotPhase() {
                return SLOT_PHASE_INTERNAL;
            }

            final int queueIdx;
            final OrbitJob orbitJob;

            RetrialEvent(int queueIdx, OrbitJob orbitJob) {
                this.queueIdx = queueIdx;
                this.orbitJob = orbitJob;
            }

            @Override
            public void actions() {
                retrialActions(queueIdx, orbitJob);
            }
        }

        /**
         * Handle customer re-entering queue after successful retrial.
         * Similar to regular arrival but skips balking check (they already decided to join).
         */
        private void arriveAtQueueFromRetrial(int queueIdx, Customer customer) {
            int classId = customer.classId;

            // Update queue length statistics
            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]++;

            // Track successful retrial
            retriedCustomers[queueIdx][classId]++;

            // Try to start service or join queue
            if (isDelayNode.get(queueIdx)) {
                // Delay node: always start service immediately
                customersInService[queueIdx]++;
                updateBusyStats(queueIdx, classId);
                currentBusyServers[queueIdx][classId]++;
                lastBusyUpdateTime[queueIdx][classId] = ssjSim.time();

                double serviceTime = generateServiceTime(queueIdx, classId);
                DelayDeparture departureEvent = new DelayDeparture(queueIdx, customer);
                departureEvent.schedule(serviceWallDelay(queueIdx, classId, serviceTime));
            } else {
                int freeServer = findFreeServer(queueIdx);
                if (freeServer >= 0) {
                    // Start service immediately
                    serverBusy[queueIdx][freeServer] = true;
                    customersInService[queueIdx]++;
                    updateBusyStats(queueIdx, classId);
                    currentBusyServers[queueIdx][classId]++;
                    lastBusyUpdateTime[queueIdx][classId] = ssjSim.time();

                    double serviceTime = generateServiceTime(queueIdx, classId);
                    Departure departureEvent = new Departure(queueIdx, freeServer, customer);
                    departureEvent.schedule(serviceWallDelay(queueIdx, classId, serviceTime));
                } else {
                    // Join wait queue
                    waitQueues[queueIdx].add(customer);

                    // Schedule reneging if configured
                    scheduleRenegingEvent(queueIdx, customer);
                }
            }
        }

        /**
         * Determines if an arriving customer should balk based on configured thresholds.
         */
        private boolean shouldBalk(int queueIdx, int classId) {
            if (!hasBalkingConfig[queueIdx][classId]) return false;

            if (sn.balkingThresholds == null) return false;
            jline.lang.nodes.Station station = sn.stations.get(serviceStations.get(queueIdx));
            Map<jline.lang.JobClass, java.util.List<jline.lang.constant.BalkingThreshold>> classMap =
                    sn.balkingThresholds.get(station);
            if (classMap == null) return false;
            jline.lang.JobClass jobClass = sn.jobclasses.get(classId);
            java.util.List<jline.lang.constant.BalkingThreshold> thresholds = classMap.get(jobClass);
            if (thresholds == null) return false;

            int currentQueueLen = getTotalCustomersAtStation(queueIdx);

            for (jline.lang.constant.BalkingThreshold threshold : thresholds) {
                if (threshold.matches(currentQueueLen)) {
                    // Within range - balk with given probability
                    return routingRng.nextDouble() < threshold.getProbability();
                }
            }
            return false;
        }

        /**
         * Checks if BAS (Blocking After Service) is enabled for the given station and class.
         * JMT convention: blocking policy is specified on the LDESTINATION queue (the queue with limited capacity).
         * When a job completes service at an upstream queue and the destination is full, the server blocks.
         */
        private boolean hasBASBlocking(int destQueueIdx, int classId) {
            if (destQueueIdx < 0 || destQueueIdx >= stationDropRule.length) return false;
            if (classId < 0 || classId >= stationDropRule[destQueueIdx].length) return false;
            return stationDropRule[destQueueIdx][classId] == DropStrategy.BlockingAfterService.getID();
        }

        /**
         * Checks if BBS (Blocking Before Service) is enabled for the given station and class.
         * JMT convention: blocking policy is specified on the LDESTINATION queue (the queue with limited capacity).
         */
        private boolean hasBBSBlocking(int destQueueIdx, int classId) {
            if (destQueueIdx < 0 || destQueueIdx >= stationDropRule.length) return false;
            if (classId < 0 || classId >= stationDropRule[destQueueIdx].length) return false;
            return stationDropRule[destQueueIdx][classId] == DropStrategy.BlockingBeforeService.getID();
        }

        /**
         * Checks if RSRD (Re-Service on Rejection with Delay) is enabled for the given station and class.
         */
        private boolean hasRSRDBlocking(int destQueueIdx, int classId) {
            if (destQueueIdx < 0 || destQueueIdx >= stationDropRule.length) return false;
            if (classId < 0 || classId >= stationDropRule[destQueueIdx].length) return false;
            return stationDropRule[destQueueIdx][classId] == DropStrategy.ReServiceOnRejection.getID();
        }

        /**
         * Checks if a destination queue has capacity to accept a new customer.
         * Returns true if the destination can accept, false if full.
         */
        private boolean destinationHasCapacity(int destQueueIdx, int destClassId) {
            if (destQueueIdx < 0 || destQueueIdx >= numServiceNodes) return true;  // Unknown destination, allow

            // Check total station capacity
            int currentTotal = getTotalCustomersAtStation(destQueueIdx);
            if (currentTotal >= bufferCapacities[destQueueIdx]) {
                return false;
            }

            // Check per-class capacity constraint
            if (classCapacities[destQueueIdx][destClassId] < Integer.MAX_VALUE) {
                int currentClassCount = currentQueueLength[destQueueIdx][destClassId];
                if (currentClassCount >= classCapacities[destQueueIdx][destClassId]) {
                    return false;
                }
            }

            return true;
        }

        /**
         * Attempts to unblock BBS-blocked servers when destination queue gets space.
         * Finds oldest blocked server across ALL SOURCES (global FIFO by block start time).
         * Called when a customer departs from a service node, freeing capacity.
         *
         * BBS (Blocking Before Service): The server at the source is blocked until
         * the destination has capacity. When unblocked, the server can start serving
         * the next waiting customer.
         */
        private void tryUnblockBBSServers(int destQueueIdx) {
            // Get all sources with blocked servers waiting for this destination
            Set<Integer> sourcesWithBlocked = bbsDestinationToSources.get(destQueueIdx);
            if (sourcesWithBlocked == null) return;
            if (sourcesWithBlocked.isEmpty()) return;

            // Keep trying to unblock servers while destination has capacity
            while (true) {
                // Find oldest blocked server across ALL sources (global FIFO by block start time)
                BBSBlockedServer oldestBlocked = null;
                int oldestSourceIdx = -1;
                double oldestTimestamp = Double.POSITIVE_INFINITY;

                List<Integer> sourcesSnapshot = new ArrayList<Integer>(sourcesWithBlocked);
                for (Integer sourceIdx : sourcesSnapshot) {
                    List<BBSBlockedServer> blockedList = bbsBlockedServers.get(sourceIdx);
                    if (blockedList == null) continue;

                    // Find first blocked server for this destination from this source
                    for (BBSBlockedServer blocked : blockedList) {
                        if (blocked.destQueueIdx == destQueueIdx) {
                            if (blocked.blockStartTime < oldestTimestamp) {
                                oldestTimestamp = blocked.blockStartTime;
                                oldestBlocked = blocked;
                                oldestSourceIdx = sourceIdx.intValue();
                            }
                            break;  // Only check first match per source (FIFO per source)
                        }
                    }
                }

                // No more blocked servers waiting for this destination
                if (oldestBlocked == null) {
                    break;
                }

                // Check if destination has capacity for this specific class
                if (!destinationHasCapacity(destQueueIdx, oldestBlocked.destClassId)) {
                    break;  // No capacity for this class
                }

                // Unblock the oldest blocked server
                int sourceQueueIdx = oldestSourceIdx;
                Customer customer = oldestBlocked.customer;
                int destClassId = oldestBlocked.destClassId;
                int serverId = oldestBlocked.serverId;
                int sourceClassId = oldestBlocked.sourceClassId;

                // Remove from blocked servers list
                List<BBSBlockedServer> blockedList = bbsBlockedServers.get(sourceQueueIdx);
                blockedList.remove(oldestBlocked);

                // Clean up empty lists and indices
                boolean stillHasDestMatch = false;
                for (BBSBlockedServer it : blockedList) {
                    if (it.destQueueIdx == destQueueIdx) {
                        stillHasDestMatch = true;
                        break;
                    }
                }
                if (!stillHasDestMatch) {
                    sourcesWithBlocked.remove(Integer.valueOf(sourceQueueIdx));
                }
                if (blockedList.isEmpty()) {
                    bbsBlockedServers.remove(sourceQueueIdx);
                }

                int sourceStationIdx = serviceStations.get(sourceQueueIdx);
                logEvent("BBS_UNBLOCK", sourceStationIdx, sourceClassId, serverId, destQueueIdx);

                // BBS: blocked job was counted at DESTINATION, so decrement there
                // Update destination queue stats before decrementing blocked count (to properly time-weight)
                updateQueueStats(destQueueIdx, destClassId);
                bbsBlockedAtDest[destQueueIdx][destClassId]--;

                // Unblock server at source
                serverBlocked[sourceQueueIdx][serverId] = false;

                // Create customer for destination
                Customer nextCustomer = new Customer(
                        destClassId,
                        classPrio[destClassId],
                        customer.systemArrivalTime,
                        ssjSim.time(),  // Queue arrival time at destination is now
                        siroRng.nextDouble(),
                        -1.0,
                        customer.jobId,
                        customer.absoluteDeadline,
                        -1,
                        null);

                // Send job to destination queue
                arriveAtQueue(destQueueIdx, nextCustomer, true);

                // Free server and start next customer at source queue (with heterogeneous tracking)
                markServerIdle(sourceQueueIdx, serverId);
                customersInService[sourceQueueIdx]--;

                // Start service for next customer in source queue (if any) - heterogeneous-aware
                if (!waitQueues[sourceQueueIdx].isEmpty()) {
                    Customer nextWaiting = waitQueues[sourceQueueIdx].poll();
                    int nextClassId = nextWaiting.classId;

                    // Find compatible server for next customer
                    ServerSelection serverSelection = findFreeServerForClass(sourceQueueIdx, nextClassId);
                    if (serverSelection.serverId >= 0) {
                        markServerBusy(sourceQueueIdx, serverSelection.serverId, serverSelection.serverTypeId);
                        customersInService[sourceQueueIdx]++;
                        nextWaiting.assignedServerType = serverSelection.serverTypeId;

                        // Flush accumulated busy time for nextClassId before updating tracker
                        updateBusyStats(sourceQueueIdx, nextClassId);
                        currentBusyServers[sourceQueueIdx][nextClassId]++;

                        double serviceTime = (nextWaiting.serviceTime > 0)
                                ? nextWaiting.serviceTime
                                : generateHeteroServiceTime(sourceQueueIdx, nextClassId, serverSelection.serverTypeId);
                        Departure departureEvent = new Departure(sourceQueueIdx, serverSelection.serverId, nextWaiting);
                        departureEvent.schedule(serviceWallDelay(sourceQueueIdx, nextClassId, serviceTime));

                        // Track for signal-based removal (G-networks)
                        if (hasRemovalSignals) {
                            inServiceJobs.put(new IntPair(sourceQueueIdx, serverSelection.serverId),
                                    new InServiceJob(nextWaiting, departureEvent));
                        }

                        logEvent("BBS_NEXT_SERVICE", sourceStationIdx, nextClassId,
                                currentQueueLength[sourceQueueIdx][nextClassId],
                                currentBusyServers[sourceQueueIdx][nextClassId]);
                    } else {
                        // No compatible server - put customer back
                        waitQueues[sourceQueueIdx].add(nextWaiting);
                    }
                }
            }

            // Clean up reverse index if no more sources blocked for this destination
            if (sourcesWithBlocked.isEmpty()) {
                bbsDestinationToSources.remove(destQueueIdx);
            }
        }

        /**
         * Attempts to admit BAS waiting jobs when destination queue gets space.
         * Finds oldest waiting job across ALL SOURCES (global FIFO by arrival time).
         * Called when a customer departs from a service node, freeing capacity.
         *
         * BAS (Blocking After Service): Jobs wait at the SOURCE occupying the server.
         * The server is blocked until the destination has capacity. When unblocked,
         * the job moves to the destination and the server resumes processing.
         */
        private void tryAdmitBASWaitingJobs(int destQueueIdx) {
            // Get all sources with waiting jobs for this destination
            Set<Integer> sourcesWithWaiting = basDestinationToSources.get(destQueueIdx);
            if (sourcesWithWaiting == null) return;
            if (sourcesWithWaiting.isEmpty()) return;

            // Keep admitting jobs while destination has capacity
            while (true) {
                // Find oldest waiting job across ALL sources (global FIFO by arrival time)
                BASWaitingJob oldestWaiting = null;
                int oldestSourceIdx = -1;
                double oldestTimestamp = Double.POSITIVE_INFINITY;

                List<Integer> sourcesSnapshot = new ArrayList<Integer>(sourcesWithWaiting);
                for (Integer sourceIdx : sourcesSnapshot) {
                    List<BASWaitingJob> outgoingBuffer = basOutgoingBuffer.get(sourceIdx);
                    if (outgoingBuffer == null) continue;

                    // Find first waiting job for this destination from this source
                    for (BASWaitingJob waiting : outgoingBuffer) {
                        if (waiting.destQueueIdx == destQueueIdx) {
                            if (waiting.arrivalTime < oldestTimestamp) {
                                oldestTimestamp = waiting.arrivalTime;
                                oldestWaiting = waiting;
                                oldestSourceIdx = sourceIdx.intValue();
                            }
                            break;  // Only check first match per source (FIFO per source)
                        }
                    }
                }

                // No more waiting jobs for this destination
                if (oldestWaiting == null) {
                    break;
                }

                // Check if destination has capacity for this specific class
                if (!destinationHasCapacity(destQueueIdx, oldestWaiting.destClassId)) {
                    break;  // No capacity for this class
                }

                // Remove from source's outgoing buffer
                int sourceQueueIdx = oldestSourceIdx;
                List<BASWaitingJob> outgoingBuffer = basOutgoingBuffer.get(sourceQueueIdx);
                outgoingBuffer.remove(oldestWaiting);

                Customer customer = oldestWaiting.customer;
                int destClassId = oldestWaiting.destClassId;
                int serverId = oldestWaiting.serverId;
                int sourceClassId = oldestWaiting.sourceClassId;

                int sourceStationIdx = serviceStations.get(sourceQueueIdx);
                logEvent("BAS_UNBLOCK", sourceStationIdx, sourceClassId, serverId, destQueueIdx);
                basUnblockCount++;

                // BAS: waiting job was counted at DESTINATION (JMT convention), so decrement there
                // Update dest queue stats before decrementing (to properly time-weight)
                updateQueueStats(destQueueIdx, destClassId);
                basBlockedAtDest[destQueueIdx][destClassId]--;

                // Unblock server at source
                serverBlocked[sourceQueueIdx][serverId] = false;

                // Create customer for destination
                // Use the original buffer arrival time so response time includes waiting in BAS buffer
                // (job is conceptually at destination from qlen perspective, so waiting counts)
                Customer nextCustomer = new Customer(
                        destClassId,
                        classPrio[destClassId],
                        customer.systemArrivalTime,
                        oldestWaiting.arrivalTime,  // When job arrived at BAS buffer (counts toward dest response time)
                        siroRng.nextDouble(),
                        -1.0,
                        customer.jobId,
                        customer.absoluteDeadline,
                        -1,
                        null);

                // Send job to destination queue
                arriveAtQueue(destQueueIdx, nextCustomer, true);

                // Free server and start next customer at source queue (with heterogeneous tracking)
                markServerIdle(sourceQueueIdx, serverId);
                customersInService[sourceQueueIdx]--;

                // Start service for next customer in source queue (if any) - heterogeneous-aware
                if (!waitQueues[sourceQueueIdx].isEmpty()) {
                    Customer nextWaiting = waitQueues[sourceQueueIdx].poll();
                    int nextClassId = nextWaiting.classId;

                    // Find compatible server for next customer
                    ServerSelection serverSelection = findFreeServerForClass(sourceQueueIdx, nextClassId);
                    if (serverSelection.serverId >= 0) {
                        markServerBusy(sourceQueueIdx, serverSelection.serverId, serverSelection.serverTypeId);
                        customersInService[sourceQueueIdx]++;
                        nextWaiting.assignedServerType = serverSelection.serverTypeId;

                        // Flush accumulated busy time for nextClassId before updating tracker
                        updateBusyStats(sourceQueueIdx, nextClassId);
                        currentBusyServers[sourceQueueIdx][nextClassId]++;

                        double serviceTime = (nextWaiting.serviceTime > 0)
                                ? nextWaiting.serviceTime
                                : generateHeteroServiceTime(sourceQueueIdx, nextClassId, serverSelection.serverTypeId);
                        Departure departureEvent = new Departure(sourceQueueIdx, serverSelection.serverId, nextWaiting);
                        departureEvent.schedule(serviceWallDelay(sourceQueueIdx, nextClassId, serviceTime));

                        // Track for signal-based removal (G-networks)
                        if (hasRemovalSignals) {
                            inServiceJobs.put(new IntPair(sourceQueueIdx, serverSelection.serverId),
                                    new InServiceJob(nextWaiting, departureEvent));
                        }

                        logEvent("BAS_NEXT_SERVICE", sourceStationIdx, nextClassId,
                                currentQueueLength[sourceQueueIdx][nextClassId],
                                currentBusyServers[sourceQueueIdx][nextClassId]);
                    } else {
                        // No compatible server - put customer back
                        waitQueues[sourceQueueIdx].add(nextWaiting);
                    }
                }

                // Clean up empty lists and indices
                boolean stillHasDestMatch = false;
                for (BASWaitingJob it : outgoingBuffer) {
                    if (it.destQueueIdx == destQueueIdx) {
                        stillHasDestMatch = true;
                        break;
                    }
                }
                if (!stillHasDestMatch) {
                    sourcesWithWaiting.remove(Integer.valueOf(sourceQueueIdx));
                }
                if (outgoingBuffer.isEmpty()) {
                    basOutgoingBuffer.remove(sourceQueueIdx);
                }
            }

            // Clean up reverse index if no more sources with waiting jobs for this destination
            if (sourcesWithWaiting.isEmpty()) {
                basDestinationToSources.remove(destQueueIdx);
            }
        }

        /**
         * Returns the effective number of servers for PS queues, accounting for load dependence.
         * For load-dependent stations, the effective server count is the lldscaling value for the current population.
         * For regular stations, returns the static numServers value.
         *
         * @param queueIdx Service node index
         * @return Effective number of servers as a double
         */
        private double getEffectivePSServerCount(int queueIdx) {
            if (isLoadDependent[queueIdx]) {
                int totalJobs = getTotalCustomersAtStation(queueIdx);
                if (totalJobs > 0) {
                    double[] scalingArray = lldScaling[queueIdx];
                    if (scalingArray != null) {
                        // lldscaling[station, n-1] gives effective server count when n jobs are present
                        int scalingIdx = Math.min(totalJobs - 1, scalingArray.length - 1);
                        return scalingArray[scalingIdx];
                    }
                }
                // When 0 jobs, use 1 as effective server count (first job gets full service)
                return 1.0;
            }
            return (double) numServers[queueIdx];
        }

        /**
         * Returns the total number of customers in a finite capacity region (all classes).
         */
        private int getTotalCustomersInRegion(int regionIdx) {
            int total = 0;
            for (int k = 0; k < numClasses; k++) {
                total += currentJobsInRegion[regionIdx][k];
            }
            return total;
        }

        /**
         * Returns the memory currently occupied in a region: sum over classes of
         * the in-region job count times the per-class memory footprint (classSize).
         * Used to enforce the region global memory budget (JMT globalMemoryConstraint).
         */
        private double getWeightedMemoryInRegion(int regionIdx) {
            return currentMemInRegion[regionIdx];
        }

        /** Registers one classId job entering a region, keeping the occupied-
         *  memory counter in sync with the job counter. */
        private void regionJobEnter(int regionIdx, int classId) {
            currentJobsInRegion[regionIdx][classId]++;
            currentMemInRegion[regionIdx] += fcRegionClassSize[regionIdx][classId];
        }

        /** Registers one classId job leaving a region (inverse of regionJobEnter). */
        private void regionJobLeave(int regionIdx, int classId) {
            currentJobsInRegion[regionIdx][classId]--;
            currentMemInRegion[regionIdx] -= fcRegionClassSize[regionIdx][classId];
        }

        /**
         * True if admitting one classId job would exceed the region's global
         * memory budget (usedMemory + incomingSize > maxMemory, JMT
         * globalMemoryConstraint semantics). The epsilon absorbs floating-point
         * accumulation now that footprints may be fractional.
         */
        private boolean exceedsRegionMemBudget(int regionIdx, int classId) {
            double regionMaxMem = (regionIdx < fcRegionGlobalMaxMem.length)
                    ? fcRegionGlobalMaxMem[regionIdx] : -1.0;
            if (regionMaxMem < 0) {
                return false;
            }
            double incomingMem = (classId < fcRegionClassSize[regionIdx].length)
                    ? fcRegionClassSize[regionIdx][classId] : 1.0;
            return getWeightedMemoryInRegion(regionIdx) + incomingMem > regionMaxMem + 1e-9;
        }

        /**
         * True if destNode is a station inside region regionIdx, i.e. a hop to
         * it does not cross the region boundary. Sinks, fork/join nodes and
         * non-member stations are outside every region.
         */
        private boolean isSameRegionDestination(int regionIdx, int destNode) {
            if (regionIdx < 0 || destNode < 0) {
                return false;
            }
            int nextQueueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
            if (nextQueueIdx < 0) {
                return false;
            }
            int stationIdx = serviceStations.get(nextQueueIdx);
            return stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) == regionIdx;
        }

        /**
         * Completes the region bookkeeping of a departure once the routing
         * destination is known. Called after regionJobLeave: on a true region
         * exit it counts the completion and releases blocked FIFO customers;
         * on an intra-region hop it does nothing, so the freed slot stays
         * reserved for the moving job (its admission gate then self-passes).
         * Releasing during the transient dip of an intra-region hop admitted a
         * blocked head into the reserved slot and re-parked the moving job,
         * which distorted per-class-cap dynamics (JMT parity fix).
         */
        private void regionExitCompleted(int queueIdx, int classId, int destNode, int destClassId) {
            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx >= fcRegionIndices.size() || fcRegionIndices.get(stationIdx) < 0) {
                return;
            }
            int regionIdx = fcRegionIndices.get(stationIdx);
            if (destClassId == classId && isSameRegionDestination(regionIdx, destNode)) {
                // same-class intra-region hop: the slot stays reserved for the
                // moving job. A class-switching hop instead crosses the border
                // in JMT (the ClassSwitch node sits outside the region), so it
                // is treated as a true exit followed by a gated re-entry.
                return;
            }
            if (warmupDone) {
                regionCompletions[regionIdx][classId]++;
            }
            tryReleaseBlockedCustomers(regionIdx);
        }

        /**
         * Handles a customer refused entry to a finite capacity region: DROP
         * discards it, WAITQ parks it in the region's blocked queue. Returns
         * the arriveAtQueue return value (false = dropped, true = blocked and
         * consumed). kind names the violated constraint for event logging,
         * e.g. "REGION_MEMORY" -> DROP_REGION_MEMORY / BLOCK_REGION_MEMORY.
         */
        private boolean rejectAtRegion(int regionIdx, int classId, int queueIdx,
                int stationIdx, Customer customer, String kind) {
            boolean shouldDrop = (regionIdx >= 0
                    && regionIdx < fcRegionDropRule.length
                    && classId < fcRegionDropRule[regionIdx].length)
                    ? fcRegionDropRule[regionIdx][classId] : true;
            if (shouldDrop) {
                if (warmupDone) {
                    droppedByRegion[regionIdx][classId]++;
                }
                logEvent("DROP_" + kind, stationIdx, classId,
                        currentQueueLength[queueIdx][classId],
                        currentBusyServers[queueIdx][classId]);
                return false;
            } else {
                fcRegionBlockedQueue[regionIdx].add(new BlockedCustomer(customer, queueIdx));
                updateRegionTimeWeightedStats(regionIdx);
                blockedInRegion[regionIdx][classId]++;
                updateQueueStats(queueIdx, classId);
                fcrBlockedAtDest[queueIdx][classId]++;
                logEvent("BLOCK_" + kind, stationIdx, classId,
                        currentQueueLength[queueIdx][classId],
                        currentBusyServers[queueIdx][classId]);
                return true;
            }
        }

        /**
         * Returns the number of customers of a specific class in a finite capacity region.
         */
        private int getClassJobsInRegion(int regionIdx, int classId) {
            if (regionIdx >= 0 && regionIdx < numRegions) {
                return currentJobsInRegion[regionIdx][classId];
            }
            return 0;
        }

        /**
         * Updates time-weighted statistics for a finite capacity region.
         * Should be called before changing currentJobsInRegion or blockedInRegion.
         * JMT semantics: FCR QLen includes only jobs inside the region (blocked jobs are NOT included).
         * In JMT, BlockingRegion.increaseOccupation() is only called when a job successfully enters.
         */
        private void updateRegionTimeWeightedStats(int regionIdx) {
            if (regionIdx < 0 || regionIdx >= numRegions) return;

            double currentTime = ssjSim.time();
            double elapsed = currentTime - lastRegionUpdateTime[regionIdx];
            if (elapsed > 0 && warmupDone) {
                for (int k = 0; k < numClasses; k++) {
                    // Blocked jobs are NOT included (they have not entered the region yet).
                    // JMT "Number of Customers" counts raw jobs; "FCR Capacity" applies
                    // classWeight (Total Weight); "FCR Memory" applies classSize.
                    double jobsForClass = (double) currentJobsInRegion[regionIdx][k];
                    totalRegionJobTime[regionIdx][k] += jobsForClass * elapsed;
                    totalRegionWeightTime[regionIdx][k] += jobsForClass
                            * fcRegionClassWeights[regionIdx][k] * elapsed;
                    totalRegionMemTime[regionIdx][k] += jobsForClass
                            * fcRegionClassSize[regionIdx][k] * elapsed;
                }
            }
            lastRegionUpdateTime[regionIdx] = currentTime;
        }

        /**
         * Updates FCR arrival rate tracking for a region.
         * Should be called when a job successfully enters a region (after updating currentJobsInRegion).
         * Tracks inter-arrival times for computing arrival rate = 1 / mean(inter-arrival time).
         */
        private void updateRegionArrivalTracking(int regionIdx, int classId) {
            if (regionIdx < 0 || regionIdx >= numRegions || !warmupDone) return;

            double currentTime = ssjSim.time();
            if (regionArrivalCount[regionIdx][classId] > 0) {
                // Not the first arrival - record inter-arrival time
                double interArrivalTime = currentTime - lastRegionArrivalTime[regionIdx][classId];
                regionInterArrivalTimeSum[regionIdx][classId] += interArrivalTime;
            }
            regionArrivalCount[regionIdx][classId]++;
            lastRegionArrivalTime[regionIdx][classId] = currentTime;
        }

        /**
         * Attempts to release blocked customers from an FCR's waiting queue after a departure.
         * Called when a customer departs from a station in a finite capacity region.
         * Releases customers in FIFO order as long as capacity permits.
         */
        private void tryReleaseBlockedCustomers(int regionIdx) {
            if (regionIdx < 0 || regionIdx >= numRegions) return;

            LinkedList<BlockedCustomer> blockedQueue = fcRegionBlockedQueue[regionIdx];
            while (!blockedQueue.isEmpty()) {
                BlockedCustomer blockedCustomer = blockedQueue.peek();
                Customer customer = blockedCustomer.customer;
                int destQueueIdx = blockedCustomer.destQueueIdx;
                int classId = customer.classId;

                // Check if we can admit this customer now
                int currentRegionTotal = getTotalCustomersInRegion(regionIdx);
                int globalMax = (regionIdx >= 0 && regionIdx < fcRegionGlobalMax.size())
                        ? fcRegionGlobalMax.get(regionIdx) : Integer.MAX_VALUE;

                // Check global region capacity
                if (currentRegionTotal >= globalMax) {
                    // Still full, can't release more
                    break;
                }

                // Check global region memory budget (weighted by per-class size)
                if (exceedsRegionMemBudget(regionIdx, classId)) {
                    // Still over the memory budget, can't release more (FIFO order preserved)
                    break;
                }

                // Check per-class region capacity
                int regionClassMax = (int) fcRegionClassMax.get(regionIdx, classId);
                if (regionClassMax < Integer.MAX_VALUE) {
                    int currentRegionClass = getClassJobsInRegion(regionIdx, classId);
                    if (currentRegionClass >= regionClassMax) {
                        // Per-class still full for this customer's class
                        // To maintain FIFO ordering properly, we stop here - blocked customers must be released in order
                        break;
                    }
                }

                // Check linear admission constraints before releasing
                if (fcRegionLinConA[regionIdx].length > 0) {
                    boolean linConViolated = false;
                    for (int c = 0; c < fcRegionLinConA[regionIdx].length; c++) {
                        double lhs = fcRegionLinConA[regionIdx][c][classId];  // contribution of released job
                        for (int r = 0; r < numClasses; r++) {
                            lhs += fcRegionLinConA[regionIdx][c][r] * currentJobsInRegion[regionIdx][r];
                        }
                        if (lhs > fcRegionLinConb[regionIdx][c]) {
                            linConViolated = true;
                            break;
                        }
                    }
                    if (linConViolated) break;
                }

                // Remove from blocked queue and route to destination
                blockedQueue.poll();

                // Update region stats before decrementing blocked count
                updateRegionTimeWeightedStats(regionIdx);
                blockedInRegion[regionIdx][classId]--;

                // Update queue stats before decrementing FCR blocked count
                updateQueueStats(destQueueIdx, classId);
                fcrBlockedAtDest[destQueueIdx][classId]--;

                // The blocked customer is released as-is: its queueArrivalTime
                // keeps the original arrival instant (JMT semantics: blocking
                // time is included in response time) and its jobId, serviceTime,
                // deadline and forked identity survive the release. Rebuilding
                // it with the short constructor dropped the jobId, so a blocked
                // job of a class that expects a reply could never pair with the
                // reply of the call it issues, and a blocked forked sibling
                // stalled its Join forever.
                Customer releasedCustomer = customer;
                // ... but order in the destination buffer by the release instant:
                // customers admitted through the gate meanwhile are served first
                releasedCustomer.orderTime = ssjSim.time();

                logEvent("UNBLOCK_REGION", serviceStations.get(destQueueIdx), classId,
                        currentQueueLength[destQueueIdx][classId],
                        currentBusyServers[destQueueIdx][classId]);

                // Now route the released customer to its destination queue
                // This will do all the normal processing (including updating region counts)
                processReleasedCustomerArrival(destQueueIdx, releasedCustomer);
            }
        }

        /**
         * Processes a released blocked customer's arrival at a queue.
         * This is called after a customer is released from an FCR blocked queue.
         * It performs the same logic as arriveAtQueue but skips the FCR capacity check
         * since we already verified capacity before releasing.
         *
         * JMT semantics: Queue length and region job counts are incremented HERE when the customer
         * actually enters (blocked customers were only counted in FCR QLen, not queue QLen).
         */
        private void processReleasedCustomerArrival(int queueIdx, Customer customer) {
            int classId = customer.classId;
            int stationIdx = serviceStations.get(queueIdx);

            // Dispatch to PS handler for Processor Sharing scheduling
            // Pass fromBlocked=true to skip queue length increment
            if (isPSScheduling(schedStrategies[queueIdx])) {
                arriveAtPSQueueFromBlocked(queueIdx, customer);
                return;
            }

            // Dispatch to preemptive LCFS handler
            // Pass fromBlocked=true to skip queue length increment
            if (isPreemptiveScheduling[queueIdx]) {
                arriveAtPreemptiveLCFSQueueFromBlocked(queueIdx, customer);
                return;
            }

            // Dispatch to polling handler
            // Pass fromBlocked=true to skip queue length increment
            if (isPollingStation[queueIdx]) {
                arriveAtPollingQueueFromBlocked(queueIdx, customer);
                return;
            }

            // For SJF/LJF, we must generate service time upon arrival to sort the queue
            SchedStrategy strategy = schedStrategies[queueIdx];
            if (strategy == SchedStrategy.SJF || strategy == SchedStrategy.LJF) {
                customer.serviceTime = generateServiceTime(queueIdx, classId);
            }

            // JMT semantics: blocked jobs are NOT counted in queue QLen, only in FCR QLen
            // Now that customer is entering the queue, increment queue length
            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]++;

            // Track max queue length reached
            int totalAtStation = getTotalCustomersAtStation(queueIdx);
            if (totalAtStation > maxQueueLengthReached) {
                maxQueueLengthReached = totalAtStation;
            }

            // Update region job counts NOW that customer is actually entering the region
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobEnter(regionIdx, classId);
                updateRegionArrivalTracking(regionIdx, classId);
            }

            if (isDelayNode.get(queueIdx)) {
                // Delay node (infinite server): always start service immediately
                customersInService[queueIdx]++;

                // Update busy time statistics before incrementing
                updateBusyStats(queueIdx, classId);
                currentBusyServers[queueIdx][classId]++;
                // Reset busy time tracking for this job - prevents including warmup time in first job of class
                lastBusyUpdateTime[queueIdx][classId] = ssjSim.time();

                double serviceTime = generateServiceTime(queueIdx, classId);
                DelayDeparture departureEvent = new DelayDeparture(queueIdx, customer);
                departureEvent.schedule(serviceWallDelay(queueIdx, classId, serviceTime));
                // Track for signal-based removal (G-networks)
                if (hasRemovalSignals) {
                    long jobId = nextDelayJobId++;
                    delayJobs.put(jobId, new DelayJob(queueIdx, customer, departureEvent));
                }
            } else {
                // Queue node: check if server is available (heterogeneous-aware)
                ServerSelection serverSelection = findFreeServerForClass(queueIdx, classId);
                int freeServer = serverSelection.serverId;
                int serverTypeId = serverSelection.serverTypeId;
                if (freeServer >= 0) {
                    // Start service immediately with heterogeneous tracking
                    markServerBusy(queueIdx, freeServer, serverTypeId);
                    customersInService[queueIdx]++;
                    customer.assignedServerType = serverTypeId;

                    // Update busy time statistics before incrementing
                    updateBusyStats(queueIdx, classId);
                    currentBusyServers[queueIdx][classId]++;
                    // Reset busy time tracking for this job - prevents including warmup time in first job of class
                    lastBusyUpdateTime[queueIdx][classId] = ssjSim.time();

                    double serviceTime = (customer.serviceTime > 0)
                            ? customer.serviceTime
                            : generateHeteroServiceTime(queueIdx, classId, serverTypeId);
                    Departure departureEvent = new Departure(queueIdx, freeServer, customer);
                    departureEvent.schedule(serviceWallDelay(queueIdx, classId, serviceTime));
                    // Track for signal-based removal (G-networks)
                    if (hasRemovalSignals) {
                        inServiceJobs.put(new IntPair(queueIdx, freeServer),
                                new InServiceJob(customer, departureEvent));
                    }
                } else {
                    // Join queue (priority queue sorts by priority, then FCFS)
                    waitQueues[queueIdx].add(customer);
                }
            }

            // Log arrival event
            logEvent("RELEASED_ARRIVAL", stationIdx, classId,
                    currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId]);
        }

        /**
         * Arrival at PS queue from blocked state.
         * JMT semantics: Queue length and region counts are incremented HERE when customer enters.
         */
        private void arriveAtPSQueueFromBlocked(int queueIdx, Customer customer) {
            int classId = customer.classId;
            double currentTime = ssjSim.time();

            // Update remaining work for all current jobs based on elapsed time
            updatePSRemainingWork(queueIdx, currentTime);

            // JMT semantics: blocked jobs are NOT counted in queue QLen, only in FCR QLen
            // Now that customer is entering the queue, increment queue length
            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]++;

            // Track total customers in service
            customersInService[queueIdx]++;

            // Update region job counts NOW that customer is actually entering the region
            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobEnter(regionIdx, classId);
                updateRegionArrivalTracking(regionIdx, classId);
            }

            // LPS: if at capacity, job waits in FCFS queue instead of entering PS service
            int lpsLimit = lpsLimits[queueIdx];
            if (lpsLimit > 0 && psJobsInService[queueIdx].size() >= lpsLimit) {
                customer.serviceTime = generateServiceTime(queueIdx, classId);
                waitQueues[queueIdx].add(customer);
                logEvent("LPS_WAIT", serviceStations.get(queueIdx), classId,
                        currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId]);
                return;
            }

            // Generate service requirement for new customer
            double serviceRequirement = generateServiceTime(queueIdx, classId);

            // Create PS customer
            PSCustomer psCustomer = new PSCustomer(
                    classId,
                    customer.priority,
                    customer.systemArrivalTime,
                    customer.queueArrivalTime,
                    serviceRequirement,
                    serviceRequirement,
                    null);

            // Update busy time statistics for PS BEFORE adding the new job
            updatePSBusyStats(queueIdx);

            // Add to jobs in service (in PS, all jobs are always "in service")
            psJobsInService[queueIdx].add(psCustomer);

            // Reschedule all departures with new rates
            rescheduleAllPSDepartures(queueIdx);

            logEvent("PS_RELEASED_ARRIVAL", serviceStations.get(queueIdx), classId,
                    currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId]);
        }

        /**
         * Arrival at preemptive LCFS queue from blocked state.
         * JMT semantics: Queue length and region counts are incremented HERE when customer enters.
         */
        private void arriveAtPreemptiveLCFSQueueFromBlocked(int queueIdx, Customer customer) {
            int classId = customer.classId;
            SchedStrategy strategy = schedStrategies[queueIdx];

            // JMT semantics: blocked jobs are NOT counted in queue QLen, only in FCR QLen
            // Now that customer is entering the queue, increment queue length
            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]++;

            // Update region job counts NOW that customer is actually entering the region
            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobEnter(regionIdx, classId);
                updateRegionArrivalTracking(regionIdx, classId);
            }

            // Track max queue length reached
            int totalAtStation = getTotalCustomersAtStation(queueIdx);
            if (totalAtStation > maxQueueLengthReached) {
                maxQueueLengthReached = totalAtStation;
            }

            // Check if we should preempt current job
            boolean shouldPreempt = shouldPreemptForLCFS(queueIdx, customer);

            if (shouldPreempt) {
                // Preempt current job
                PreemptiveCustomer victimJob = findJobToPreempt(queueIdx, customer, strategy);
                if (victimJob != null) {
                    preemptJob(queueIdx, victimJob);
                    // Inherit victim's server type
                    startPreemptiveService(queueIdx, customer, victimJob.serverId, victimJob.assignedServerType);
                } else {
                    waitQueues[queueIdx].add(customer);
                }
            } else {
                // Find free server (heterogeneous-aware)
                ServerSelection serverSelection = findFreeServerForClass(queueIdx, classId);
                if (serverSelection.serverId >= 0) {
                    startPreemptiveService(queueIdx, customer, serverSelection.serverId, serverSelection.serverTypeId);
                } else {
                    waitQueues[queueIdx].add(customer);
                }
            }

            logEvent("LCFS_RELEASED_ARRIVAL", serviceStations.get(queueIdx), classId,
                    currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId]);
        }

        /**
         * Arrival at polling queue from blocked state.
         * JMT semantics: Queue length and region counts are incremented HERE when customer enters.
         */
        private void arriveAtPollingQueueFromBlocked(int queueIdx, Customer customer) {
            int classId = customer.classId;

            // JMT semantics: blocked jobs are NOT counted in queue QLen, only in FCR QLen
            // Now that customer is entering the queue, increment queue length
            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]++;

            // Track max queue length reached
            int totalAtStation = getTotalCustomersAtStation(queueIdx);
            if (totalAtStation > maxQueueLengthReached) {
                maxQueueLengthReached = totalAtStation;
            }

            // Update region job counts NOW that customer is actually entering the region
            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobEnter(regionIdx, classId);
                updateRegionArrivalTracking(regionIdx, classId);
            }

            // Add to per-class polling queue
            pollingQueues[queueIdx][classId].add(customer);

            logEvent("POLL_RELEASED_ARRIVAL", serviceStations.get(queueIdx), classId,
                    currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId]);

            // If server is idle and not in switchover, try to start service
            if (customersInService[queueIdx] == 0 && !pollingInSwitchover[queueIdx]) {
                pollingTryStartService(queueIdx);
            }
        }

    // ==================== PART 8 TRANSLATION (Kotlin lines 10501-12000) ====================

        // ---------- Forward stubs to be filled in subsequent chunks (Kotlin lines 12000+) ----------

        /** Checks if setup/delayoff is enabled for the given queue and class. Kotlin line 12496. */
        private boolean isSetupDelayoffEnabled(int queueIdx, int classId) {
            return hasSetupDelayoff[queueIdx]
                    && setupGens[queueIdx][classId] != null
                    && delayoffGens[queueIdx][classId] != null;
        }

        /**
         * Initiates server delayoff (teardown) phase.
         * The server transitions to DELAYOFF state and schedules a DelayoffCompletion event.
         * Kotlin line 12540.
         */
        private void startServerDelayoff(int queueIdx, int serverId, int classId) {
            // Update server state to DELAYOFF
            serverState[queueIdx][serverId] = ServerState.DELAYOFF;

            // Track which class initiated this delayoff
            serverLastClass[queueIdx][serverId] = classId;

            // Generate delayoff time from distribution
            double delayoffTime = generateDelayoffTime(queueIdx, classId);

            // Update delayoff statistics
            updateDelayoffStats(queueIdx, classId);
            currentServersInDelayoff[queueIdx][classId]++;
            lastDelayoffUpdateTime[queueIdx][classId] = ssjSim.time();

            // Schedule delayoff completion event
            DelayoffCompletion delayoffEvent = new DelayoffCompletion(queueIdx, serverId);
            delayoffEvent.schedule(delayoffTime);

            // Store event reference for potential cancellation
            pendingDelayoffEvents[queueIdx][serverId] = delayoffEvent;

            // Log delayoff start
            int stationIdx = serviceStations.get(queueIdx);
            logEvent("DELAYOFF_START", stationIdx, classId, serverId, 0);
        }

        // ---------- Departure event body (Kotlin lines 10383-10745) ----------

        private void departureActions(int queueIdx, int serverId, Customer customer) {
            trackEvent();
            int classId = customer.classId;

            // A job that self-loops back to this same LCFS-family station must not
            // be re-selected by the very server it just freed: under last-come
            // first-served the freed server serves the newest job WAITING at the
            // completion instant, and the completing job re-joins as newest only
            // afterwards. Its re-arrival is therefore deferred until after the
            // server has picked the next job (below). Under FCFS this is moot
            // (the re-looped job sorts to the back), so only LCFS needs it.
            Customer deferredSelfLoopArrival = null;

            // Record queue response time
            double queueResponseTime = ssjSim.time() - customer.queueArrivalTime;
            responseTimeTally[queueIdx][classId].add(queueResponseTime);
            responseTimeSamples[queueIdx][classId].add(queueResponseTime);
            completedCustomers[queueIdx][classId]++;

            // Record queue tardiness (relative to deadline)
            double queueTardiness = Math.max(0.0, ssjSim.time() - customer.absoluteDeadline);
            tardinessTally[queueIdx][classId].add(queueTardiness);

            // Check event count for stopping/warmup/MSER sampling
            checkEventCountStop();

            // Update queue length statistics
            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            // Clean up and re-scale class-dependence tracking
            if ((hasCd || hasLld) && sdDepartureEvents != null) {
                sdDepartureEvents[queueIdx].remove(Integer.valueOf(serverId));
                sdInServiceCustomers[queueIdx].remove(Integer.valueOf(serverId));
                rescaleStateDepInServiceJobs(queueIdx);
            }

            // Update region job counts and try releasing blocked customers
            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
                // completion count and blocked-FIFO release are deferred to
                // regionExitCompleted once the routing destination is known:
                // intra-region hops must not release into the transient slot
            }

            // LQN phase-2: the spawn continuation takes over the freed slot
            maybeSpawnOnCompletion(queueIdx, classId, customer);

            // Try to admit blocked jobs now that this queue has space
            tryAdmitBASWaitingJobs(queueIdx);
            tryUnblockBBSServers(queueIdx);

            // Update busy time before decrementing
            updateBusyStats(queueIdx, classId);
            currentBusyServers[queueIdx][classId]--;

            // Clean up in-service tracking (G-networks)
            if (hasRemovalSignals) {
                inServiceJobs.remove(new IntPair(queueIdx, serverId));
            }

            // Log departure event
            logEvent("DEPARTURE", stationIdx, classId,
                    currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId]);

            // Check if this class expects a reply signal (synchronous call semantics)
            boolean expectsReply = synchCallReplyClass[classId] >= 0;
            long jobId = customer.jobId;

            // Route customer to next destination
            int currentNode = serviceNodes.get(queueIdx);
            RoutingResult routingResult = selectDestination(currentNode, classId);
            int destNode = routingResult.destNode;
            int destClassId = routingResult.destClassId;
            // A job that enters a blocking class by a class switch carries no
            // job id of its own: the outgoing customer inherits the departing
            // job's id, which is -1 whenever the departing class expects no
            // reply. Its server would then never block, since the pending
            // reply is keyed by job id. Mint one here so that the reply
            // arriving later can be paired with this call. Kept separate from
            // jobId, which keys the block registered by the departing class.
            long outJobId = mintCallJobId(destClassId, jobId, classId);
            // FCR: on a true region exit, count the completion and release
            // blocked customers; a no-op for intra-region hops
            regionExitCompleted(queueIdx, classId, destNode, destClassId);

            if (destNode >= 0 && sinkNodes.contains(Integer.valueOf(destNode))) {
                // Customer leaves system (count in original class for consistency)
                double systemResponseTime = ssjSim.time() - customer.systemArrivalTime;
                systemResponseTimeTally[classId].add(systemResponseTime);

                double systemTardiness = Math.max(0.0, ssjSim.time() - customer.absoluteDeadline);
                systemTardinessTally[classId].add(systemTardiness);

                systemCompletedCustomers[classId]++;
            } else if (destNode >= 0 && forkNodes.contains(Integer.valueOf(destNode))) {
                long parentJobId = nextJobId++;
                handleForkArrival(destNode, parentJobId, destClassId, customer.systemArrivalTime);
            } else if (destNode >= 0 && joinNodes.contains(Integer.valueOf(destNode))) {
                IntDoubleKey forkedJobKey = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
                ForkedJob forkedJob = forkedCustomerMap.remove(forkedJobKey);
                if (forkedJob != null) {
                    ForkedJob updatedForkedJob = new ForkedJob(
                            forkedJob.forkJobId,
                            forkedJob.parentJobId,
                            destClassId,
                            forkedJob.priority,
                            forkedJob.systemArrivalTime,
                            ssjSim.time(),
                            forkedJob.randomRank);
                    handleJoinArrival(destNode, updatedForkedJob);
                } else {
                    handleUnknownJoinArrival(destNode, destClassId, customer.systemArrivalTime);
                }
            } else if (destNode >= 0) {
                int nextQueueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
                if (nextQueueIdx >= 0) {
                    if ((hasNegativeSignals && isNegativeSignal[destClassId]) ||
                            (hasCatastropheSignals && isCatastropheSignal[destClassId])) {
                        handleNegativeSignalArrival(nextQueueIdx, destClassId, customer.systemArrivalTime);
                    } else if (hasReplySignals && isReplySignal[destClassId]) {
                        IntDoubleKey replyFjKey = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
                        ForkedJob replyFj = forkedCustomerMap.remove(replyFjKey);
                        handleReplySignalArrival(nextQueueIdx, destClassId, jobId, customer.systemArrivalTime, replyFj);
                    } else {
                        // Check blocking policies at BOTH source and destination queues
                        boolean hasBBS = hasBBSBlocking(queueIdx, classId)
                                || hasBBSBlocking(nextQueueIdx, destClassId);
                        boolean hasBAS = hasBASBlocking(queueIdx, classId)
                                || hasBASBlocking(nextQueueIdx, destClassId);
                        boolean destHasCap = destinationHasCapacity(nextQueueIdx, destClassId);

                        if (hasBBS && !destHasCap) {
                            // BBS (Blocking Before Service)
                            BBSBlockedServer blockedServer = new BBSBlockedServer(
                                    customer, nextQueueIdx, destClassId,
                                    queueIdx, serverId, classId,
                                    ssjSim.time());

                            List<BBSBlockedServer> blockedList = bbsBlockedServers.get(Integer.valueOf(queueIdx));
                            if (blockedList == null) {
                                blockedList = new ArrayList<BBSBlockedServer>();
                                bbsBlockedServers.put(Integer.valueOf(queueIdx), blockedList);
                            }
                            blockedList.add(blockedServer);

                            Set<Integer> sourcesForDest = bbsDestinationToSources.get(Integer.valueOf(nextQueueIdx));
                            if (sourcesForDest == null) {
                                sourcesForDest = new HashSet<Integer>();
                                bbsDestinationToSources.put(Integer.valueOf(nextQueueIdx), sourcesForDest);
                            }
                            sourcesForDest.add(Integer.valueOf(queueIdx));

                            serverBlocked[queueIdx][serverId] = true;

                            updateQueueStats(nextQueueIdx, destClassId);
                            bbsBlockedAtDest[nextQueueIdx][destClassId]++;
                            updateQueueStats(nextQueueIdx, destClassId);

                            logEvent("BBS_BLOCK", stationIdx, classId, serverId, nextQueueIdx);
                            return;
                        } else if (hasBAS && !destHasCap) {
                            // BAS (Blocking After Service)
                            BASWaitingJob waitingJob = new BASWaitingJob(
                                    customer, nextQueueIdx, destClassId,
                                    queueIdx, serverId, classId,
                                    ssjSim.time());

                            List<BASWaitingJob> outgoingBuffer = basOutgoingBuffer.get(Integer.valueOf(queueIdx));
                            if (outgoingBuffer == null) {
                                outgoingBuffer = new ArrayList<BASWaitingJob>();
                                basOutgoingBuffer.put(Integer.valueOf(queueIdx), outgoingBuffer);
                            }
                            outgoingBuffer.add(waitingJob);

                            Set<Integer> sourcesForDest = basDestinationToSources.get(Integer.valueOf(nextQueueIdx));
                            if (sourcesForDest == null) {
                                sourcesForDest = new HashSet<Integer>();
                                basDestinationToSources.put(Integer.valueOf(nextQueueIdx), sourcesForDest);
                            }
                            sourcesForDest.add(Integer.valueOf(queueIdx));

                            serverBlocked[queueIdx][serverId] = true;

                            updateQueueStats(nextQueueIdx, destClassId);
                            basBlockedAtDest[nextQueueIdx][destClassId]++;

                            logEvent("BAS_BLOCK", stationIdx, classId, serverId, nextQueueIdx);
                            basBlockCount++;
                            return;
                        } else if (!destHasCap && destClassId >= 0 && destClassId < isClosedClass.length
                                && isClosedClass[destClassId]) {
                            // A CLOSED job whose destination is full and that has NO explicit
                            // blocking rule (BAS/BBS handled above). It cannot be dropped --
                            // population conservation is an invariant -- so it is not routed;
                            // instead it REPETITIVELY re-serves at its current server, which
                            // is the default closed-class blocking the analytical engines use
                            // (CTMC/SSA disable the departure, so the exponential job simply
                            // re-attempts at rate mu -- repetitive service). Without this the
                            // job fell through to arriveAtQueue below and was DROPPED, draining
                            // the closed population (BUG-81: 3-queue N=6 classCap=2 drained
                            // 6 -> 2). The job stays as classId at queueIdx (it never switched
                            // class or moved), and a fresh service time is drawn -- exact for
                            // exponential service by memorylessness.
                            // Undo the completion bookkeeping that departureActions
                            // recorded at its top: this service completion did NOT release
                            // the job, it only re-serves it, so it must not count as a
                            // throughput departure (else Tput is inflated by every
                            // re-attempt). The queue-length decrement at the top is undone
                            // by the re-increment below, and the top's regionJobLeave is
                            // undone by the regionJobEnter below, so the state is unchanged.
                            completedCustomers[queueIdx][classId]--;

                            updateBusyStats(queueIdx, classId);
                            currentBusyServers[queueIdx][classId]++;
                            updateQueueStats(queueIdx, classId);
                            currentQueueLength[queueIdx][classId]++;

                            if (stationIdx < fcRegionIndices.size()
                                    && fcRegionIndices.get(stationIdx) >= 0) {
                                int regionIdxRs = fcRegionIndices.get(stationIdx);
                                updateRegionTimeWeightedStats(regionIdxRs);
                                regionJobEnter(regionIdxRs, classId);
                            }

                            // Preserve the ORIGINAL queue-arrival time so that when the
                            // job finally leaves, its residence time spans the whole stay
                            // including every re-service, not just the last segment.
                            Customer rsCustomer = new Customer(
                                    classId, classPrio[classId],
                                    customer.systemArrivalTime, customer.queueArrivalTime,
                                    siroRng.nextDouble(), -1.0, jobId,
                                    customer.absoluteDeadline, -1, null);
                            double rsServiceTime = generateHeteroServiceTime(
                                    queueIdx, classId, customer.assignedServerType);
                            Departure rsDeparture = new Departure(queueIdx, serverId, rsCustomer);
                            rsDeparture.schedule(rsServiceTime);

                            if ((hasCd || hasLld) && sdDepartureEvents != null) {
                                sdDepartureEvents[queueIdx].put(Integer.valueOf(serverId), rsDeparture);
                                sdInServiceCustomers[queueIdx].put(Integer.valueOf(serverId), rsCustomer);
                            }
                            if (hasRemovalSignals) {
                                inServiceJobs.put(new IntPair(queueIdx, serverId),
                                        new InServiceJob(rsCustomer, rsDeparture));
                            }

                            logEvent("RS_BLOCK", stationIdx, classId, serverId, nextQueueIdx);
                            return;
                        } else {
                            // Check for immediate feedback (self-loop staying in service)
                            boolean isImmediateFeedback = nextQueueIdx == queueIdx
                                    && sn.immfeed != null
                                    && stationIdx >= 0 && stationIdx < sn.immfeed.getNumRows()
                                    && destClassId >= 0 && destClassId < sn.immfeed.getNumCols()
                                    && sn.immfeed.get(stationIdx, destClassId) > 0.0;

                            if (isImmediateFeedback) {
                                // The completing job is fed back into the same server as
                                // destClassId, holding the server instead of re-queueing.
                                // departureActions has already removed it from the (classId)
                                // busy-server AND queue-length counts, so re-add it under
                                // destClassId for both. Unconditional: a same-class self-loop
                                // is a net no-op, a class-switching one moves the held job
                                // between classes. (Previously this double-decremented
                                // busy[classId] -> Util[class_a]=0, and never restored the
                                // queue length -> QLen[destClass]=0.)
                                updateBusyStats(queueIdx, destClassId);
                                currentBusyServers[queueIdx][destClassId]++;
                                updateQueueStats(queueIdx, destClassId);
                                currentQueueLength[queueIdx][destClassId]++;

                                // FCR: restore the region occupancy removed by the
                                // departure bookkeeping; the fed-back job never
                                // leaves the station (regionExitCompleted was a
                                // no-op for this self-loop)
                                if (stationIdx < fcRegionIndices.size()
                                        && fcRegionIndices.get(stationIdx) >= 0) {
                                    int regionIdxFb = fcRegionIndices.get(stationIdx);
                                    updateRegionTimeWeightedStats(regionIdxFb);
                                    regionJobEnter(regionIdxFb, destClassId);
                                }

                                Customer nextCustomer = new Customer(
                                        destClassId, classPrio[destClassId],
                                        customer.systemArrivalTime, ssjSim.time(),
                                        siroRng.nextDouble(), -1.0, outJobId,
                                        customer.absoluteDeadline, -1, null);

                                double serviceTime = generateHeteroServiceTime(
                                        queueIdx, destClassId, customer.assignedServerType);
                                Departure departureEvent = new Departure(queueIdx, serverId, nextCustomer);
                                departureEvent.schedule(serviceWallDelay(queueIdx, destClassId, serviceTime));

                                if ((hasCd || hasLld) && sdDepartureEvents != null) {
                                    sdDepartureEvents[queueIdx].put(Integer.valueOf(serverId), departureEvent);
                                    sdInServiceCustomers[queueIdx].put(Integer.valueOf(serverId), nextCustomer);
                                }

                                if (hasRemovalSignals) {
                                    inServiceJobs.put(new IntPair(queueIdx, serverId),
                                            new InServiceJob(nextCustomer, departureEvent));
                                }

                                logEvent("IMMFEED", stationIdx, destClassId,
                                        currentQueueLength[queueIdx][destClassId],
                                        currentBusyServers[queueIdx][destClassId]);
                                return;
                            }

                            // Normal routing: destination has capacity or no blocking policy
                            IntDoubleKey forkedJobKey = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
                            ForkedJob forkedJob = forkedCustomerMap.remove(forkedJobKey);
                            Customer nextCustomer = new Customer(
                                    destClassId, classPrio[destClassId],
                                    customer.systemArrivalTime, ssjSim.time(),
                                    siroRng.nextDouble(), -1.0, outJobId,
                                    customer.absoluteDeadline, -1, null);
                            if (forkedJob != null) {
                                ForkedJob nextForkedJob = new ForkedJob(
                                        forkedJob.forkJobId,
                                        forkedJob.parentJobId,
                                        destClassId,
                                        forkedJob.priority,
                                        forkedJob.systemArrivalTime,
                                        ssjSim.time(),
                                        forkedJob.randomRank);
                                arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob);
                            } else if (nextQueueIdx == queueIdx
                                    && isLCFSFamily(schedStrategies[queueIdx])) {
                                // Self-loop at an LCFS station: defer the re-arrival
                                // until the freed server has selected the next job.
                                deferredSelfLoopArrival = nextCustomer;
                            } else {
                                arriveAtQueue(nextQueueIdx, nextCustomer);
                            }
                        }
                    }
                }
            }

            // Handle synchronous call blocking
            boolean classSwitchedToReply = hasReplySignals && destClassId >= 0 && isReplySignal[destClassId];
            if (expectsReply && jobId >= 0 && !pendingReplyMap.containsKey(Long.valueOf(jobId)) && !classSwitchedToReply) {
                PendingReply pendingReply = new PendingReply(
                        jobId, classId, queueIdx, serverId, ssjSim.time());
                pendingReplyMap.put(Long.valueOf(jobId), pendingReply);

                serverBlocked[queueIdx][serverId] = true;

                updateQueueStats(queueIdx, classId);
                currentBlockedServers[queueIdx][classId]++;

                logEvent("SYNCH_BLOCK", stationIdx, classId, serverId, 0);
                if (deferredSelfLoopArrival != null) {
                    arriveAtQueue(queueIdx, deferredSelfLoopArrival);
                }
                return;
            }

            // Free server (normal case - no reply expected, with heterogeneous tracking)
            markServerIdle(queueIdx, serverId);
            customersInService[queueIdx]--;

            // Start service for next customer in queue (if any)
            if (!waitQueues[queueIdx].isEmpty()) {
                Customer nextCustomer = waitQueues[queueIdx].poll();
                int nextClassId = nextCustomer.classId;

                ServerSelection serverSelection = findFreeServerForClass(queueIdx, nextClassId);
                if (serverSelection.serverId >= 0) {
                    markServerBusy(queueIdx, serverSelection.serverId, serverSelection.serverTypeId);
                    customersInService[queueIdx]++;
                    nextCustomer.assignedServerType = serverSelection.serverTypeId;

                    updateBusyStats(queueIdx, nextClassId);
                    currentBusyServers[queueIdx][nextClassId]++;

                    double serviceTime;
                    if (nextCustomer.serviceTime > 0) {
                        serviceTime = nextCustomer.serviceTime;
                    } else {
                        serviceTime = generateHeteroServiceTime(queueIdx, nextClassId, serverSelection.serverTypeId);
                    }
                    Departure departureEvent = new Departure(queueIdx, serverSelection.serverId, nextCustomer);
                    departureEvent.schedule(serviceWallDelay(queueIdx, nextClassId, serviceTime));

                    if ((hasCd || hasLld) && sdDepartureEvents != null) {
                        sdDepartureEvents[queueIdx].put(Integer.valueOf(serverSelection.serverId), departureEvent);
                        sdInServiceCustomers[queueIdx].put(Integer.valueOf(serverSelection.serverId), nextCustomer);
                    }

                    if (hasRemovalSignals) {
                        inServiceJobs.put(new IntPair(queueIdx, serverSelection.serverId),
                                new InServiceJob(nextCustomer, departureEvent));
                    }
                } else {
                    waitQueues[queueIdx].add(nextCustomer);
                }
            } else {
                // Queue is empty - initiate delayoff if enabled
                if (isSetupDelayoffEnabled(queueIdx, classId)) {
                    startServerDelayoff(queueIdx, serverId, classId);
                }
            }

            // Now re-admit the deferred self-loop job (LCFS): it joins as the
            // newest waiter, after the freed server has already taken the job
            // that was newest at the completion instant.
            if (deferredSelfLoopArrival != null) {
                arriveAtQueue(queueIdx, deferredSelfLoopArrival);
            }
        }

        // ---------- PreemptiveDeparture event body (Kotlin lines 10756-10882) ----------

        private void preemptiveDepartureActions(int queueIdx, int serverId, PreemptiveCustomer customer) {
            trackEvent();
            int classId = customer.classId;
            double currentTime = ssjSim.time();

            // Remove from in-service tracking - if not present, job was preempted
            boolean wasInService = preemptiveJobsInService[queueIdx].remove(customer);
            if (!wasInService) {
                return;
            }

            // Record queue response time
            double queueResponseTime = currentTime - customer.queueArrivalTime;
            if (warmupDone) {
                responseTimeTally[queueIdx][classId].add(queueResponseTime);
                responseTimeSamples[queueIdx][classId].add(queueResponseTime);
                completedCustomers[queueIdx][classId]++;
            }

            checkEventCountStop();

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
                // completion count and blocked-FIFO release are deferred to
                // regionExitCompleted once the routing destination is known:
                // intra-region hops must not release into the transient slot
            }

            // LQN phase-2: the spawn continuation takes over the freed slot
            maybeSpawnOnCompletion(queueIdx, classId, customer.queueArrivalTime, customer.randomRank);

            updateBusyStats(queueIdx, classId);
            currentBusyServers[queueIdx][classId]--;

            markServerIdle(queueIdx, serverId);
            customersInService[queueIdx]--;

            logEvent("DEPARTURE", stationIdx, classId,
                    currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId]);

            int currentNode = serviceNodes.get(queueIdx);
            RoutingResult routingResult = selectDestination(currentNode, classId);
            int destNode = routingResult.destNode;
            int destClassId = routingResult.destClassId;
            // FCR: on a true region exit, count the completion and release
            // blocked customers; a no-op for intra-region hops
            regionExitCompleted(queueIdx, classId, destNode, destClassId);

            if (destNode >= 0 && sinkNodes.contains(Integer.valueOf(destNode))) {
                double systemResponseTime = currentTime - customer.systemArrivalTime;
                if (warmupDone) {
                    systemResponseTimeTally[classId].add(systemResponseTime);
                    systemCompletedCustomers[classId]++;
                }
            } else if (destNode >= 0 && forkNodes.contains(Integer.valueOf(destNode))) {
                long parentJobId = nextJobId++;
                handleForkArrival(destNode, parentJobId, destClassId, customer.systemArrivalTime);
            } else if (destNode >= 0 && joinNodes.contains(Integer.valueOf(destNode))) {
                IntDoubleKey forkedJobKey = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
                ForkedJob forkedJob = forkedCustomerMap.remove(forkedJobKey);
                if (forkedJob != null) {
                    ForkedJob updatedForkedJob = new ForkedJob(
                            forkedJob.forkJobId,
                            forkedJob.parentJobId,
                            destClassId,
                            forkedJob.priority,
                            forkedJob.systemArrivalTime,
                            currentTime,
                            forkedJob.randomRank);
                    handleJoinArrival(destNode, updatedForkedJob);
                } else {
                    handleUnknownJoinArrival(destNode, destClassId, customer.systemArrivalTime);
                }
            } else if (destNode >= 0) {
                int nextQueueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
                if (nextQueueIdx >= 0) {
                    IntDoubleKey forkedJobKey = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
                    ForkedJob forkedJob = forkedCustomerMap.remove(forkedJobKey);
                    Customer nextCustomer = new Customer(
                            destClassId, classPrio[destClassId],
                            customer.systemArrivalTime, currentTime, customer.randomRank);
                    if (forkedJob != null) {
                        ForkedJob nextForkedJob = new ForkedJob(
                                forkedJob.forkJobId,
                                forkedJob.parentJobId,
                                destClassId,
                                forkedJob.priority,
                                forkedJob.systemArrivalTime,
                                currentTime,
                                forkedJob.randomRank);
                        arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob);
                    } else {
                        arriveAtQueue(nextQueueIdx, nextCustomer);
                    }
                }
            }

            // Start service for next customer in queue (if any) - heterogeneous-aware
            if (!waitQueues[queueIdx].isEmpty()) {
                Customer nextCustomer = waitQueues[queueIdx].poll();
                int nextClassId = nextCustomer.classId;
                ServerSelection serverSelection = findFreeServerForClass(queueIdx, nextClassId);
                if (serverSelection.serverId >= 0) {
                    startPreemptiveService(queueIdx, nextCustomer, serverSelection.serverId, serverSelection.serverTypeId);
                } else {
                    waitQueues[queueIdx].add(nextCustomer);
                }
            }
        }

        // ---------- SetupCompletion event (Kotlin lines 10889-10917) ----------

        /** Setup completion event when a server finishes its cold start phase. */
        private final class SetupCompletion extends SimEvent {
            final int queueIdx;
            final int serverId;
            final int classId;

            SetupCompletion(int queueIdx, int serverId, int classId) {
                this.queueIdx = queueIdx;
                this.serverId = serverId;
                this.classId = classId;
            }

            @Override
            public void actions() {
                trackEvent();
                updateSetupStats(queueIdx, classId);
                currentServersInSetup[queueIdx][classId]--;

                serverState[queueIdx][serverId] = ServerState.ACTIVE;

                int stationIdx = serviceStations.get(queueIdx);
                logEvent("SETUP_COMPLETE", stationIdx, classId, serverId, 0);

                if (!waitQueues[queueIdx].isEmpty()) {
                    Customer nextCustomer = waitQueues[queueIdx].poll();
                    startService(queueIdx, serverId, nextCustomer, -1);
                }
            }
        }

        // ---------- DelayoffCompletion event (Kotlin lines 10923-10952) ----------

        /** Delayoff completion event when a server finishes its teardown phase. */
        private final class DelayoffCompletion extends SimEvent {
            final int queueIdx;
            final int serverId;

            DelayoffCompletion(int queueIdx, int serverId) {
                this.queueIdx = queueIdx;
                this.serverId = serverId;
            }

            @Override
            public void actions() {
                trackEvent();
                if (serverState[queueIdx][serverId] != ServerState.DELAYOFF) {
                    return;
                }

                int classId = serverLastClass[queueIdx][serverId];
                if (classId >= 0) {
                    updateDelayoffStats(queueIdx, classId);
                    currentServersInDelayoff[queueIdx][classId]--;
                }

                serverState[queueIdx][serverId] = ServerState.OFF;

                pendingDelayoffEvents[queueIdx][serverId] = null;

                int stationIdx = serviceStations.get(queueIdx);
                logEvent("DELAYOFF_COMPLETE", stationIdx, classId, serverId, 0);
            }
        }

        // ---------- DelayDeparture event body (Kotlin lines 10963-11086) ----------

        private void delayDepartureActions(int queueIdx, Customer customer) {
            trackEvent();
            int classId = customer.classId;

            double responseTime = ssjSim.time() - customer.queueArrivalTime;
            responseTimeTally[queueIdx][classId].add(responseTime);
            responseTimeSamples[queueIdx][classId].add(responseTime);
            completedCustomers[queueIdx][classId]++;

            checkEventCountStop();

            updateQueueStats(queueIdx, classId);
            currentQueueLength[queueIdx][classId]--;

            int stationIdx = serviceStations.get(queueIdx);
            if (stationIdx < fcRegionIndices.size() && fcRegionIndices.get(stationIdx) >= 0) {
                int regionIdx = fcRegionIndices.get(stationIdx);
                updateRegionTimeWeightedStats(regionIdx);
                regionJobLeave(regionIdx, classId);
                // completion count and blocked-FIFO release are deferred to
                // regionExitCompleted once the routing destination is known:
                // intra-region hops must not release into the transient slot
            }

            // LQN phase-2: the spawn continuation takes over the freed slot
            maybeSpawnOnCompletion(queueIdx, classId, customer);

            updateBusyStats(queueIdx, classId);
            currentBusyServers[queueIdx][classId]--;

            customersInService[queueIdx]--;

            // Clean up delay job tracking (G-networks)
            if (hasRemovalSignals) {
                Long jobIdToRemove = null;
                for (Map.Entry<Long, DelayJob> entry : delayJobs.entrySet()) {
                    DelayJob dj = entry.getValue();
                    if (dj.queueIdx == queueIdx && dj.customer == customer) {
                        jobIdToRemove = entry.getKey();
                        break;
                    }
                }
                if (jobIdToRemove != null) {
                    delayJobs.remove(jobIdToRemove);
                }
            }

            logEvent("DELAY_DEPARTURE", stationIdx, classId,
                    currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId]);

            int currentNode = serviceNodes.get(queueIdx);
            RoutingResult routingResult = selectDestination(currentNode, classId);
            int destNode = routingResult.destNode;
            int destClassId = routingResult.destClassId;
            // FCR: on a true region exit, count the completion and release
            // blocked customers; a no-op for intra-region hops
            regionExitCompleted(queueIdx, classId, destNode, destClassId);

            if (destNode >= 0 && sinkNodes.contains(Integer.valueOf(destNode))) {
                double systemResponseTime = ssjSim.time() - customer.systemArrivalTime;
                systemResponseTimeTally[classId].add(systemResponseTime);

                double systemTardiness = Math.max(0.0, ssjSim.time() - customer.absoluteDeadline);
                systemTardinessTally[classId].add(systemTardiness);

                systemCompletedCustomers[classId]++;
            } else if (destNode >= 0 && forkNodes.contains(Integer.valueOf(destNode))) {
                long parentJobId = nextJobId++;
                handleForkArrival(destNode, parentJobId, destClassId, customer.systemArrivalTime);
            } else if (destNode >= 0 && joinNodes.contains(Integer.valueOf(destNode))) {
                IntDoubleKey forkedJobKey = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
                ForkedJob forkedJob = forkedCustomerMap.remove(forkedJobKey);
                if (forkedJob != null) {
                    ForkedJob updatedForkedJob = new ForkedJob(
                            forkedJob.forkJobId,
                            forkedJob.parentJobId,
                            destClassId,
                            forkedJob.priority,
                            forkedJob.systemArrivalTime,
                            ssjSim.time(),
                            forkedJob.randomRank);
                    handleJoinArrival(destNode, updatedForkedJob);
                } else {
                    handleUnknownJoinArrival(destNode, destClassId, customer.systemArrivalTime);
                }
            } else if (destNode >= 0) {
                int nextQueueIdx = serviceNodes.indexOf(Integer.valueOf(destNode));
                if (nextQueueIdx >= 0) {
                    if ((hasNegativeSignals && isNegativeSignal[destClassId]) ||
                            (hasCatastropheSignals && isCatastropheSignal[destClassId])) {
                        handleNegativeSignalArrival(nextQueueIdx, destClassId, customer.systemArrivalTime);
                    } else if (hasReplySignals && isReplySignal[destClassId]) {
                        IntDoubleKey replyFjKey = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
                        ForkedJob replyFj = forkedCustomerMap.remove(replyFjKey);
                        handleReplySignalArrival(nextQueueIdx, destClassId, customer.jobId, customer.systemArrivalTime, replyFj);
                    } else {
                        IntDoubleKey forkedJobKey = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
                        ForkedJob forkedJob = forkedCustomerMap.remove(forkedJobKey);
                        // See the departure path of a queueing station: a job
                        // that enters a blocking class by a class switch needs
                        // a job id of its own, or the pending reply cannot be
                        // keyed and its server never blocks.
                        long outJobId = mintCallJobId(destClassId, customer.jobId, customer.classId);
                        Customer nextCustomer = new Customer(
                                destClassId, classPrio[destClassId],
                                customer.systemArrivalTime, ssjSim.time(),
                                siroRng.nextDouble(), -1.0, outJobId,
                                customer.absoluteDeadline, -1, null);
                        if (forkedJob != null) {
                            ForkedJob nextForkedJob = new ForkedJob(
                                    forkedJob.forkJobId,
                                    forkedJob.parentJobId,
                                    destClassId,
                                    forkedJob.priority,
                                    forkedJob.systemArrivalTime,
                                    ssjSim.time(),
                                    forkedJob.randomRank);
                            arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob);
                        } else {
                            arriveAtQueue(nextQueueIdx, nextCustomer);
                        }
                    }
                }
            }
        }

        // ---------- EndOfWarmup event (Kotlin lines 11092-11101) ----------

        /** End of warmup event - resets statistics and schedules end of simulation. */
        private final class EndOfWarmup extends SimEvent {
            @Override
            protected double slotPhase() {
                return SLOT_PHASE_BOOKKEEPING;
            }

            final double remainingSimTime;

            EndOfWarmup(double remainingSimTime) {
                this.remainingSimTime = remainingSimTime;
            }

            @Override
            public void actions() {
                resetStatistics();
                warmupDone = true;
                new EndOfSimulation().schedule(remainingSimTime);
            }
        }

        // ---------- resetStatistics() body (Kotlin lines 11106-11144) ----------

        private void resetStatisticsImpl() {
            warmupEndTime = ssjSim.time();
            // Flush any pending busy time before resetting
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    updateBusyStats(qIdx, k);
                }
            }
            // Discard warmup cache statistics; cache contents themselves are kept.
            if (cacheStates != null) {
                for (int cacheNodeIdx : cacheNodes) {
                    CacheStateInfo cs = cacheStates[cacheNodeIdx];
                    if (cs == null) continue;
                    java.util.Arrays.fill(cs.totalHits, 0L);
                    java.util.Arrays.fill(cs.totalMisses, 0L);
                    for (int k = 0; k < cs.hitsPerList.length; k++) {
                        java.util.Arrays.fill(cs.hitsPerList[k], 0L);
                    }
                    for (int it = 0; it < cs.itemLevelTime.length; it++) {
                        java.util.Arrays.fill(cs.itemLevelTime[it], 0.0);
                    }
                    cs.lastContentUpdateTime = ssjSim.time();
                    cs.occupancyStartTime = ssjSim.time();
                    if (cs.totalDelayedHits != null) {
                        java.util.Arrays.fill(cs.totalDelayedHits, 0L);
                    }
                    cs.totalDelayedHitWait = 0.0;
                    cs.totalFetchTime = 0.0;
                    cs.completedFetches = 0L;
                }
            }
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    responseTimeTally[qIdx][k].init();
                    completedCustomers[qIdx][k] = 0;
                    totalQueueTime[qIdx][k] = 0.0;
                    lastQueueUpdateTime[qIdx][k] = ssjSim.time();
                    totalBusyTime[qIdx][k] = 0.0;
                    lastBusyUpdateTime[qIdx][k] = ssjSim.time();
                    totalBlockingTime[qIdx][k] = 0.0;
                }
            }
            for (int k = 0; k < numClasses; k++) {
                systemResponseTimeTally[k].init();
                systemCompletedCustomers[k] = 0;
            }
            for (int joinListIdx = 0; joinListIdx < joinNodes.size(); joinListIdx++) {
                for (int k = 0; k < numClasses; k++) {
                    joinResponseTimeTally[joinListIdx][k].init();
                    joinCompletions[joinListIdx][k] = 0;
                    totalJoinQueueTime[joinListIdx][k] = 0.0;
                    lastJoinUpdateTime[joinListIdx][k] = ssjSim.time();
                    arrivedAtJoin[joinListIdx][k] = 0;
                    droppedByJoin[joinListIdx][k] = 0;
                }
            }
        }

        // ---------- updateBusyStats() body (Kotlin lines 11149-11157) ----------

        private void updateBusyStatsImpl(int queueIdx, int classId) {
            double currentTime = ssjSim.time();
            double elapsed = currentTime - lastBusyUpdateTime[queueIdx][classId];
            if (elapsed > 0) {
                int busy = currentBusyServers[queueIdx][classId];
                totalBusyTime[queueIdx][classId] += busy * elapsed;
            }
            lastBusyUpdateTime[queueIdx][classId] = currentTime;
        }

        // ---------- updateSetupStats() body (Kotlin lines 11164-11173) ----------

        private void updateSetupStats(int queueIdx, int classId) {
            if (!warmupDone) return;

            double currentTime = ssjSim.time();
            double elapsed = currentTime - lastSetupUpdateTime[queueIdx][classId];
            if (elapsed > 0) {
                totalSetupTime[queueIdx][classId] += currentServersInSetup[queueIdx][classId] * elapsed;
                lastSetupUpdateTime[queueIdx][classId] = currentTime;
            }
        }

        // ---------- updateDelayoffStats() body (Kotlin lines 11180-11189) ----------

        private void updateDelayoffStats(int queueIdx, int classId) {
            if (!warmupDone) return;

            double currentTime = ssjSim.time();
            double elapsed = currentTime - lastDelayoffUpdateTime[queueIdx][classId];
            if (elapsed > 0) {
                totalDelayoffTime[queueIdx][classId] += currentServersInDelayoff[queueIdx][classId] * elapsed;
                lastDelayoffUpdateTime[queueIdx][classId] = currentTime;
            }
        }

        // ---------- updatePSBusyStats() body (Kotlin lines 11195-11223) ----------

        private void updatePSBusyStatsImpl(int queueIdx) {
            double currentTime = ssjSim.time();
            double elapsed = currentTime - psLastBusyUpdateTime[queueIdx];
            if (elapsed <= 0) {
                psLastBusyUpdateTime[queueIdx] = currentTime;
                return;
            }

            List<PSCustomer> jobs = psJobsInService[queueIdx];
            if (!jobs.isEmpty()) {
                SchedStrategy strategy = schedStrategies[queueIdx];
                double c = getEffectivePSServerCount(queueIdx);

                double[] rates = calculatePSRates(queueIdx, strategy, jobs, c);

                for (int idx = 0; idx < jobs.size(); idx++) {
                    PSCustomer job = jobs.get(idx);
                    double rate = rates[idx];
                    if (rate > 0) {
                        totalBusyTime[queueIdx][job.classId] += rate * elapsed;
                    }
                }
            }

            psLastBusyUpdateTime[queueIdx] = currentTime;
        }

        // ---------- getActualSimTime() (Kotlin lines 11228-11230) ----------

        private double getActualSimTime() {
            return ssjSim.time() - warmupEndTime;
        }

        // ---------- EndOfSimulation actions (Kotlin lines 11236-11271) ----------

        private void endOfSimulationActions() {
            // Final update of queue and busy statistics
            for (int qIdx = 0; qIdx < numServiceNodes; qIdx++) {
                if (isPSScheduling(schedStrategies[qIdx])) {
                    updatePSBusyStats(qIdx);
                    for (int k = 0; k < numClasses; k++) {
                        updateQueueStats(qIdx, k);
                    }
                } else {
                    for (int k = 0; k < numClasses; k++) {
                        updateQueueStats(qIdx, k);
                        updateBusyStats(qIdx, k);
                    }
                }
            }

            // Apply MSER-5 truncation to determine warmup period
            if (mserEnabled) {
                applyMSER5Truncation();
            }

            closeTracing();
            closeLoggers();

            if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                System.out.printf("\b\b\b\b\b\b\b %6d", options.samples);
                System.out.println();
            }

            ssjSim.stop();
        }

        // ---------- ProgressEvent (Kotlin lines 11277-11301) ----------

        /** Progress reporting event. */
        private final class ProgressEvent extends SimEvent {
            @Override
            protected double slotPhase() {
                return SLOT_PHASE_BOOKKEEPING;
            }

            boolean printedStart = false;

            @Override
            public void actions() {
                double t = ssjSim.time();
                int samples = (int) t;
                int totalSamples = options.samples;

                if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                    if (!printedStart) {
                        System.out.printf("LDES samples: %6d ", samples);
                        System.out.flush();
                        printedStart = true;
                    } else {
                        System.out.printf("\b\b\b\b\b\b\b %6d", samples);
                        System.out.flush();
                    }
                }

                double updateInterval = totalSamples / 50.0;
                if (updateInterval > 0 && t + updateInterval <= totalSamples) {
                    schedule(updateInterval);
                }
            }
        }

        // ---------- RenegingEvent body (Kotlin lines 11315-11357) ----------

        private void renegingActions(int queueIdx, Customer customer, ImpatientKey impatientKey) {
            trackEvent();
            ImpatientCustomer impatient = waitingImpatientCustomers.remove(impatientKey);
            if (impatient == null) return;

            int classId = customer.classId;
            double currentTime = ssjSim.time();

            // Remove customer from wait queue
            boolean removed = false;
            Iterator<Customer> it = waitQueues[queueIdx].iterator();
            while (it.hasNext()) {
                Customer c = it.next();
                if (c.classId == classId && c.systemArrivalTime == customer.systemArrivalTime) {
                    it.remove();
                    removed = true;
                    break;
                }
            }

            if (removed) {
                updateQueueStats(queueIdx, classId);
                currentQueueLength[queueIdx][classId]--;

                if (warmupDone) {
                    renegedCustomers[queueIdx][classId]++;
                    double waitTime = currentTime - customer.queueArrivalTime;
                    totalRenegingWaitTime[queueIdx][classId] += waitTime;
                }

                int stationIdx = serviceStations.get(queueIdx);
                if (!fcRegionIndices.isEmpty() && stationIdx < fcRegionIndices.size()
                        && fcRegionIndices.get(stationIdx) >= 0) {
                    int regionIdx = fcRegionIndices.get(stationIdx);
                    updateRegionTimeWeightedStats(regionIdx);
                    regionJobLeave(regionIdx, classId);
                    // a reneged job exits the region for real: release blocked customers
                    tryReleaseBlockedCustomers(regionIdx);
                }

                if (options.verbose == VerboseLevel.DEBUG) {
                    System.out.println("RENEGE: Queue " + queueIdx + " class " + classId
                            + " at time " + currentTime + ", queue length now "
                            + currentQueueLength[queueIdx][classId]);
                }
            }
        }

        // ---------- RetrialEvent body (Kotlin lines 11368-11428) ----------

        private void retrialActions(int queueIdx, OrbitJob orbitJob) {
            trackEvent();
            int classId = orbitJob.customer.classId;
            double currentTime = ssjSim.time();

            // Remove from orbit tracking
            Iterator<OrbitJob> it = orbitJobs[queueIdx].iterator();
            while (it.hasNext()) {
                if (it.next() == orbitJob) {
                    it.remove();
                    break;
                }
            }

            updateOrbitTimeStats(queueIdx, classId);
            currentOrbitSize[queueIdx][classId]--;

            int currentTotal = getTotalCustomersAtStation(queueIdx);
            int capacity = bufferCapacities[queueIdx];

            if (currentTotal < capacity) {
                Customer newCustomer = new Customer(
                        orbitJob.customer.classId,
                        orbitJob.customer.priority,
                        orbitJob.customer.systemArrivalTime,
                        currentTime,
                        siroRng.nextDouble(),
                        -1.0,
                        orbitJob.customer.jobId,
                        orbitJob.customer.absoluteDeadline,
                        -1,
                        null);

                if (warmupDone) {
                    retriedCustomers[queueIdx][classId]++;
                }

                if (options.verbose == VerboseLevel.DEBUG) {
                    System.out.println("RETRIAL_SUCCESS: Queue " + queueIdx + " class "
                            + classId + " at time " + currentTime);
                }

                arriveAtQueueFromRetrial(queueIdx, newCustomer);
            } else {
                int newAttempts = orbitJob.retrialAttempts + 1;
                if (orbitJob.maxAttempts >= 0 && newAttempts >= orbitJob.maxAttempts) {
                    if (warmupDone) {
                        maxRetriesExceeded[queueIdx][classId]++;
                    }
                    if (options.verbose == VerboseLevel.DEBUG) {
                        System.out.println("RETRIAL_DROPPED: Queue " + queueIdx + " class "
                                + classId + " at time " + currentTime + " after "
                                + newAttempts + " attempts");
                    }
                } else {
                    scheduleRetrial(queueIdx, orbitJob.customer, newAttempts, orbitJob.maxAttempts);
                    if (options.verbose == VerboseLevel.DEBUG) {
                        System.out.println("RETRIAL_RESCHEDULE: Queue " + queueIdx + " class "
                                + classId + " attempt " + newAttempts + " at time " + currentTime);
                    }
                }
            }
        }

        // ---------- selectDestination (Kotlin lines 11437-11446) ----------

        private RoutingResult selectDestination(int fromNode, int classId) {
            if (!loggerNodes.isEmpty() || !routerNodes.isEmpty()
                    || !classSwitchNodes.isEmpty() || !cacheNodes.isEmpty()) {
                long jobId = nextJobId++;
                return routeThroughPassthroughNodes(fromNode, classId, jobId);
            }
            return selectDestinationWithClassSwitch(fromNode, classId);
        }

        // ---------- handleForkArrival (Kotlin lines 11457-11519) ----------

        private void handleForkArrival(int forkNodeIdx, long parentJobId, int classId,
                                       double systemArrivalTime) {
            int forkListIdx = forkNodes.indexOf(Integer.valueOf(forkNodeIdx));
            if (forkListIdx < 0) return;

            int fanOut = forkFanOut[forkListIdx];
            double currentTime = ssjSim.time();

            int R = numClasses;
            int I = numNodes;
            List<RoutingResult> outputDestinations = new ArrayList<RoutingResult>();

            for (int toNode = 0; toNode < I; toNode++) {
                for (int toClass = 0; toClass < R; toClass++) {
                    double prob = sn.rtnodes.get(forkNodeIdx * R + classId, toNode * R + toClass);
                    if (prob > 0) {
                        outputDestinations.add(new RoutingResult(toNode, toClass));
                    }
                }
            }

            if (outputDestinations.isEmpty()) {
                return;
            }

            int totalTasks = outputDestinations.size() * fanOut;

            ForkJobInfo forkInfo = new ForkJobInfo(
                    parentJobId, classId, systemArrivalTime, forkNodeIdx, totalTasks);
            forkJobInfoMap.put(Long.valueOf(parentJobId), forkInfo);

            for (RoutingResult dest : outputDestinations) {
                for (int taskIdx = 0; taskIdx < fanOut; taskIdx++) {
                    long forkedJobId = nextForkedJobId++;
                    forkedJobParentMap.put(Long.valueOf(forkedJobId), Long.valueOf(parentJobId));

                    ForkedJob forkedJob = new ForkedJob(
                            forkedJobId,
                            parentJobId,
                            dest.destClassId,
                            classPrio[dest.destClassId],
                            systemArrivalTime,
                            currentTime,
                            siroRng.nextDouble());

                    routeForkedJob(forkedJob, dest.destNode);
                }
            }
        }

        // ---------- routeForkedJob (Kotlin lines 11525-11600) ----------

        private void routeForkedJob(ForkedJob forkedJob, int destNode) {
            int currentNode = destNode;
            int currentClass = forkedJob.classId;
            int maxIterations = 100;

            while (maxIterations > 0) {
                maxIterations--;

                if (sinkNodes.contains(Integer.valueOf(currentNode))) {
                    return;
                }

                if (joinNodes.contains(Integer.valueOf(currentNode))) {
                    handleJoinArrival(currentNode, forkedJob);
                    return;
                }

                if (forkNodes.contains(Integer.valueOf(currentNode))) {
                    handleForkArrival(currentNode, forkedJob.forkJobId, currentClass,
                            forkedJob.systemArrivalTime);
                    return;
                }

                int queueIdx = serviceNodes.indexOf(Integer.valueOf(currentNode));
                if (queueIdx >= 0) {
                    // A sibling entering a class that issues a synchronous call
                    // needs a job id of its own, or the blocking branch at its
                    // departure can never engage (fork replication creates the
                    // siblings without one).
                    long siblingJobId = mintCallJobId(currentClass, -1L, -1);
                    Customer customer = new Customer(
                            currentClass,
                            forkedJob.priority,
                            forkedJob.systemArrivalTime,
                            ssjSim.time(),
                            forkedJob.randomRank,
                            -1.0, siblingJobId,
                            Double.POSITIVE_INFINITY, -1, null);
                    arriveAtQueueForked(queueIdx, customer, forkedJob);
                    return;
                }

                if (loggerNodes.contains(Integer.valueOf(currentNode))) {
                    logJobPassage(currentNode, currentClass, forkedJob.forkJobId);
                    RoutingResult result = selectDestinationWithClassSwitch(currentNode, currentClass);
                    if (result.destNode < 0) return;
                    currentNode = result.destNode;
                    currentClass = result.destClassId;
                    continue;
                }

                if (routerNodes.contains(Integer.valueOf(currentNode))) {
                    RoutingResult result = selectDestinationWithClassSwitch(currentNode, currentClass);
                    if (result.destNode < 0) return;
                    currentNode = result.destNode;
                    currentClass = result.destClassId;
                    continue;
                }

                if (classSwitchNodes.contains(Integer.valueOf(currentNode))) {
                    RoutingResult result = selectDestinationWithClassSwitch(currentNode, currentClass);
                    if (result.destNode < 0) return;
                    currentNode = result.destNode;
                    currentClass = result.destClassId;
                    continue;
                }

                return;
            }
        }

        // ---------- arriveAtQueueForked (Kotlin lines 11606-11619) ----------

        private void arriveAtQueueForked(int queueIdx, Customer customer, ForkedJob forkedJob) {
            if (isPSScheduling(schedStrategies[queueIdx])) {
                arriveAtPSQueue(queueIdx, customer, forkedJob);
                return;
            }

            IntDoubleKey key = new IntDoubleKey(queueIdx, customer.queueArrivalTime, customer.randomRank);
            forkedCustomerMap.put(key, forkedJob);

            arriveAtQueue(queueIdx, customer);
        }

        // ---------- handleJoinArrival (Kotlin lines 11631-11778) ----------

        private void handleJoinArrival(int joinNodeIdx, ForkedJob forkedJob) {
            long parentJobId = forkedJob.parentJobId;
            ForkJobInfo forkInfo = forkJobInfoMap.get(Long.valueOf(parentJobId));
            if (forkInfo == null) {
                // The parent already synchronized (quorum/PARTIAL join fired and its
                // record was removed); this late sibling is discarded. Count it as a
                // join loss, post-warmup, mirroring arrivedAtJoin bookkeeping.
                if (warmupDone) {
                    int jlIdx = joinNodes.indexOf(Integer.valueOf(joinNodeIdx));
                    if (jlIdx >= 0) {
                        droppedByJoin[jlIdx][forkedJob.classId]++;
                    }
                }
                return;
            }

            int joinListIdx = joinNodes.indexOf(Integer.valueOf(joinNodeIdx));
            if (joinListIdx < 0) return;

            int classId = forkedJob.classId;
            double currentTime = ssjSim.time();

            if (forkInfo.firstJoinArrivalTime < 0) {
                forkInfo.firstJoinArrivalTime = currentTime;
            }
            forkInfo.forkedJobJoinArrivalTimes.add(Double.valueOf(currentTime));
            forkInfo.forkedJobJoinClasses.add(Integer.valueOf(classId));

            if (warmupDone) {
                double elapsed = currentTime - lastJoinUpdateTime[joinListIdx][classId];
                if (elapsed > 0) {
                    totalJoinQueueTime[joinListIdx][classId] +=
                            currentJoinQueueLength[joinListIdx][classId] * elapsed;
                }
                lastJoinUpdateTime[joinListIdx][classId] = currentTime;
                currentJoinQueueLength[joinListIdx][classId]++;
                arrivedAtJoin[joinListIdx][classId]++;
            }

            forkInfo.completedTasks++;

            JoinStrategy strategy = joinStrategies[joinListIdx][classId];
            int required = joinRequired[joinListIdx][classId];

            boolean syncComplete;
            if (strategy == JoinStrategy.STD) {
                syncComplete = forkInfo.completedTasks >= forkInfo.totalTasks;
            } else if (strategy == JoinStrategy.PARTIAL || strategy == JoinStrategy.Quorum) {
                if (required > 0) {
                    syncComplete = forkInfo.completedTasks >= required;
                } else {
                    syncComplete = forkInfo.completedTasks >= forkInfo.totalTasks;
                }
            } else if (strategy == JoinStrategy.Guard) {
                syncComplete = forkInfo.completedTasks >= forkInfo.totalTasks;
            } else {
                syncComplete = forkInfo.completedTasks >= forkInfo.totalTasks;
            }

            if (syncComplete) {
                if (warmupDone) {
                    for (int i = 0; i < forkInfo.forkedJobJoinArrivalTimes.size(); i++) {
                        double arrivalTime = forkInfo.forkedJobJoinArrivalTimes.get(i);
                        int forkedClassId = forkInfo.forkedJobJoinClasses.get(i);
                        double joinResponseTime = currentTime - arrivalTime;
                        joinResponseTimeTally[joinListIdx][forkedClassId].add(joinResponseTime);
                    }

                    Map<Integer, Integer> classCountMap = new HashMap<Integer, Integer>();
                    for (Integer forkedClassId : forkInfo.forkedJobJoinClasses) {
                        Integer prev = classCountMap.get(forkedClassId);
                        classCountMap.put(forkedClassId, Integer.valueOf((prev == null ? 0 : prev.intValue()) + 1));
                    }

                    for (Map.Entry<Integer, Integer> entry : classCountMap.entrySet()) {
                        int forkedClassId = entry.getKey().intValue();
                        int count = entry.getValue().intValue();
                        double elapsed2 = currentTime - lastJoinUpdateTime[joinListIdx][forkedClassId];
                        if (elapsed2 > 0) {
                            totalJoinQueueTime[joinListIdx][forkedClassId] +=
                                    currentJoinQueueLength[joinListIdx][forkedClassId] * elapsed2;
                        }
                        lastJoinUpdateTime[joinListIdx][forkedClassId] = currentTime;

                        currentJoinQueueLength[joinListIdx][forkedClassId] -= count;
                        if (currentJoinQueueLength[joinListIdx][forkedClassId] < 0) {
                            currentJoinQueueLength[joinListIdx][forkedClassId] = 0;
                        }
                    }

                    joinCompletions[joinListIdx][forkInfo.parentClassId]++;
                }

                forkJobInfoMap.remove(Long.valueOf(parentJobId));

                RoutingResult routingResult = selectDestination(joinNodeIdx, forkInfo.parentClassId);
                int destNode = routingResult.destNode;
                int destClassId = routingResult.destClassId;

                if (destNode >= 0 && !sinkNodes.contains(Integer.valueOf(destNode))
                        && destClassId == forkInfo.parentClassId) {
                    boolean isQueueOrDelay = serviceNodes.contains(Integer.valueOf(destNode));
                    if (isQueueOrDelay) {
                        int stationIdx = (int) sn.nodeToStation.get(destNode);
                        if (sn.rates != null) {
                            double rateVal = sn.rates.get(stationIdx, destClassId);
                            if (Double.isNaN(rateVal)) {
                                destNode = -1;
                            }
                        }
                    }
                }

                if (destNode >= 0 && sinkNodes.contains(Integer.valueOf(destNode))) {
                    double systemResponseTime = ssjSim.time() - forkInfo.parentSystemArrivalTime;
                    systemResponseTimeTally[forkInfo.parentClassId].add(systemResponseTime);
                    systemCompletedCustomers[forkInfo.parentClassId]++;
                } else if (destNode >= 0) {
                    routeParentJobFromJoin(destNode, destClassId, forkInfo);
                }
            }
        }

        // ---------- routeParentJobFromJoin (Kotlin lines 11783-11876) ----------

        private void routeParentJobFromJoin(int destNode, int destClassId, ForkJobInfo forkInfo) {
            double currentTime = ssjSim.time();

            int currentNode = destNode;
            int currentClass = destClassId;
            int maxIterations = 100;

            while (maxIterations > 0) {
                maxIterations--;

                if (sinkNodes.contains(Integer.valueOf(currentNode))) {
                    double systemResponseTime = currentTime - forkInfo.parentSystemArrivalTime;
                    systemResponseTimeTally[forkInfo.parentClassId].add(systemResponseTime);
                    systemCompletedCustomers[forkInfo.parentClassId]++;
                    return;
                }

                if (forkNodes.contains(Integer.valueOf(currentNode))) {
                    long newParentJobId = nextJobId++;
                    handleForkArrival(currentNode, newParentJobId, currentClass,
                            forkInfo.parentSystemArrivalTime);
                    return;
                }

                if (joinNodes.contains(Integer.valueOf(currentNode))) {
                    RoutingResult routingResult = selectDestinationWithClassSwitch(currentNode, currentClass);
                    if (routingResult.destNode < 0) return;
                    currentNode = routingResult.destNode;
                    currentClass = routingResult.destClassId;
                    continue;
                }

                if (placeNodes.contains(Integer.valueOf(currentNode))) {
                    handlePlaceArrival(currentNode, currentClass, forkInfo.parentSystemArrivalTime);
                    return;
                }

                if (transitionNodes.contains(Integer.valueOf(currentNode))) {
                    checkAndFireTransitions();
                    return;
                }

                int queueIdx = serviceNodes.indexOf(Integer.valueOf(currentNode));
                if (queueIdx >= 0) {
                    Customer customer = new Customer(
                            currentClass,
                            classPrio[currentClass],
                            forkInfo.parentSystemArrivalTime,
                            currentTime,
                            siroRng.nextDouble());
                    arriveAtQueue(queueIdx, customer);
                    return;
                }

                if (loggerNodes.contains(Integer.valueOf(currentNode))) {
                    logJobPassage(currentNode, currentClass, forkInfo.parentJobId);
                    RoutingResult result = selectDestinationWithClassSwitch(currentNode, currentClass);
                    if (result.destNode < 0) return;
                    currentNode = result.destNode;
                    currentClass = result.destClassId;
                    continue;
                }

                if (routerNodes.contains(Integer.valueOf(currentNode))) {
                    RoutingResult result = selectDestinationWithClassSwitch(currentNode, currentClass);
                    if (result.destNode < 0) return;
                    currentNode = result.destNode;
                    currentClass = result.destClassId;
                    continue;
                }

                if (classSwitchNodes.contains(Integer.valueOf(currentNode))) {
                    RoutingResult result = selectDestinationWithClassSwitch(currentNode, currentClass);
                    if (result.destNode < 0) return;
                    currentNode = result.destNode;
                    currentClass = result.destClassId;
                    continue;
                }

                return;
            }
        }

        // ---------- handleDestinationIfFork (Kotlin lines 11882-11888) ----------

        private boolean handleDestinationIfFork(int destNode, int classId, double systemArrivalTime) {
            if (forkNodes.contains(Integer.valueOf(destNode))) {
                long parentJobId = nextJobId++;
                handleForkArrival(destNode, parentJobId, classId, systemArrivalTime);
                return true;
            }
            return false;
        }

        // ---------- isJoinNode (Kotlin lines 11896-11898) ----------

        private boolean isJoinNode(int destNode) {
            return joinNodes.contains(Integer.valueOf(destNode));
        }

        // ---------- handleUnknownJoinArrival (Kotlin lines 11904-11935) ----------

        private void handleUnknownJoinArrival(int joinNodeIdx, int classId, double systemArrivalTime) {
            int joinListIdx = joinNodes.indexOf(Integer.valueOf(joinNodeIdx));
            if (joinListIdx < 0) return;

            int forkListIdx = joinToForkMap[joinListIdx];
            if (forkListIdx < 0) return;

            int forkNodeIdx = forkNodes.get(forkListIdx);

            Long matchingParentId = null;
            ForkJobInfo matchingForkInfo = null;
            for (Map.Entry<Long, ForkJobInfo> entry : forkJobInfoMap.entrySet()) {
                ForkJobInfo info = entry.getValue();
                if (info.forkNodeIdx == forkNodeIdx
                        && Math.abs(info.parentSystemArrivalTime - systemArrivalTime) < 1e-9) {
                    matchingParentId = entry.getKey();
                    matchingForkInfo = info;
                    break;
                }
            }

            if (matchingForkInfo != null) {
                ForkedJob syntheticForkedJob = new ForkedJob(
                        nextForkedJobId++,
                        matchingParentId.longValue(),
                        classId,
                        classPrio[classId],
                        systemArrivalTime,
                        ssjSim.time(),
                        0.0);
                handleJoinArrival(joinNodeIdx, syntheticForkedJob);
            }
        }

        // ---------- findFreeServer impl (Kotlin lines 11944-11977) ----------

        private int findFreeServerImpl(int queueIdx) {
            if (!hasSetupDelayoff[queueIdx]) {
                for (int i = 0; i < serverBusy[queueIdx].length; i++) {
                    if (!serverBusy[queueIdx][i]) return i;
                }
                return -1;
            }

            // First pass: Find ACTIVE idle servers
            for (int i = 0; i < serverBusy[queueIdx].length; i++) {
                if (!serverBusy[queueIdx][i] && serverState[queueIdx][i] == ServerState.ACTIVE) {
                    return i;
                }
            }

            // Second pass: Find DELAYOFF servers
            for (int i = 0; i < serverBusy[queueIdx].length; i++) {
                if (serverState[queueIdx][i] == ServerState.DELAYOFF) {
                    return i;
                }
            }

            // Third pass: Find OFF servers
            for (int i = 0; i < serverBusy[queueIdx].length; i++) {
                if (serverState[queueIdx][i] == ServerState.OFF) {
                    return i;
                }
            }

            return -1;
        }

        // ---------- findFreeServerForClass impl (start, Kotlin lines 11994-12000) ----------
        // NOTE: full implementation completes in PART 9 (Kotlin lines 12000+); this stub
        // handles the homogeneous case and falls through to subsequent translation chunks
        // for heterogeneous queues.

        private ServerSelection findFreeServerForClassImpl(int queueIdx, int classId) {
            if (numServerTypes[queueIdx] == 0) {
                int serverId = findFreeServer(queueIdx);
                return new ServerSelection(serverId, -1);
            }
            // Heterogeneous case - delegate to PART 9 helper (defined in next chunk)
            return findFreeServerForClassHetero(queueIdx, classId);
        }

        /** Heterogeneous server selection (Kotlin lines 12001-12019). */
        private ServerSelection findFreeServerForClassHetero(int queueIdx, int classId) {
            // Get compatible server types with available capacity
            List<Integer> compatibleTypes = getCompatibleServerTypes(queueIdx, classId);
            if (compatibleTypes.isEmpty()) {
                return new ServerSelection(-1, -1);
            }

            // Apply scheduling policy to select server type
            HeteroSchedPolicy policy = (heteroSchedPolicy != null && queueIdx < heteroSchedPolicy.length
                    && heteroSchedPolicy[queueIdx] != null)
                    ? heteroSchedPolicy[queueIdx]
                    : HeteroSchedPolicy.ORDER;
            int selectedType = selectServerType(queueIdx, classId, compatibleTypes, policy);

            if (selectedType < 0) {
                return new ServerSelection(-1, -1);
            }

            // Find free server within selected type
            int serverId = findFreeServerOfType(queueIdx, selectedType);
            return new ServerSelection(serverId, selectedType);
        }

        // ==================== PART 9 TRANSLATION (Kotlin lines 12001-13500) ====================

        /**
         * Gets list of server types compatible with given class that have available capacity.
         * Kotlin line 12024.
         */
        private List<Integer> getCompatibleServerTypes(int queueIdx, int classId) {
            int nTypes = numServerTypes[queueIdx];
            List<Integer> result = new ArrayList<Integer>();
            for (int typeId = 0; typeId < nTypes; typeId++) {
                if (serverCompat[queueIdx][typeId][classId]) {
                    int totalServers = serversPerType[queueIdx][typeId];
                    int busyServers = busyCountPerType[queueIdx][typeId];
                    if (busyServers < totalServers) {
                        result.add(typeId);
                    }
                }
            }
            return result;
        }

        /**
         * Finds a free server within a specific server type.
         * Handles setup/delayoff state if enabled.
         * Kotlin line 12045.
         */
        private int findFreeServerOfType(int queueIdx, int typeId) {
            // Find the range of server IDs for this type
            int startId = 0;
            for (int t = 0; t < typeId; t++) {
                startId += serversPerType[queueIdx][t];
            }
            int endId = startId + serversPerType[queueIdx][typeId];

            if (!hasSetupDelayoff[queueIdx]) {
                // Simple case: find first free server of this type
                for (int serverId = startId; serverId < endId; serverId++) {
                    if (!serverBusy[queueIdx][serverId]) {
                        return serverId;
                    }
                }
            } else {
                // With setup/delayoff: prioritize ACTIVE > DELAYOFF > OFF
                // First pass: ACTIVE idle servers
                for (int serverId = startId; serverId < endId; serverId++) {
                    if (!serverBusy[queueIdx][serverId]
                            && serverState[queueIdx][serverId] == ServerState.ACTIVE) {
                        return serverId;
                    }
                }
                // Second pass: DELAYOFF servers
                for (int serverId = startId; serverId < endId; serverId++) {
                    if (serverState[queueIdx][serverId] == ServerState.DELAYOFF) {
                        return serverId;
                    }
                }
                // Third pass: OFF servers
                for (int serverId = startId; serverId < endId; serverId++) {
                    if (serverState[queueIdx][serverId] == ServerState.OFF) {
                        return serverId;
                    }
                }
            }
            return -1;
        }

        /**
         * Selects server type according to the heterogeneous scheduling policy.
         * Kotlin line 12087.
         */
        private int selectServerType(int queueIdx, int classId,
                                     List<Integer> compatibleTypes,
                                     HeteroSchedPolicy policy) {
            if (compatibleTypes.isEmpty()) return -1;
            if (compatibleTypes.size() == 1) return compatibleTypes.get(0);

            switch (policy) {
                case ORDER:
                    return selectServerTypeORDER(compatibleTypes);
                case ALIS:
                    return selectServerTypeALIS(queueIdx, compatibleTypes);
                case ALFS:
                    return selectServerTypeALFS(queueIdx, compatibleTypes);
                case FAIRNESS:
                    return selectServerTypeFAIRNESS(queueIdx, compatibleTypes);
                case FSF:
                    return selectServerTypeFSF(queueIdx, classId, compatibleTypes);
                case RAIS:
                    return selectServerTypeRAIS(compatibleTypes);
                default:
                    return selectServerTypeORDER(compatibleTypes);
            }
        }

        /** ORDER: First compatible server type in definition order. Kotlin line 12107. */
        private int selectServerTypeORDER(List<Integer> compatibleTypes) {
            if (compatibleTypes.isEmpty()) return -1;
            int min = Integer.MAX_VALUE;
            for (int t : compatibleTypes) {
                if (t < min) min = t;
            }
            return (min == Integer.MAX_VALUE) ? -1 : min;
        }

        /** ALIS: Round-robin among compatible types. Kotlin line 12112. */
        private int selectServerTypeALIS(int queueIdx, List<Integer> compatibleTypes) {
            List<Integer> order = serverTypeOrder[queueIdx];
            for (int i = 0; i < order.size(); i++) {
                int typeId = order.get(i);
                if (compatibleTypes.contains(typeId)) {
                    // Move used type to end for round-robin
                    order.remove(Integer.valueOf(typeId));
                    order.add(typeId);
                    return typeId;
                }
            }
            return compatibleTypes.isEmpty() ? -1 : compatibleTypes.get(0);
        }

        /** ALFS: Prefer types with fewer compatible classes (least flexible first). Kotlin line 12126. */
        private int selectServerTypeALFS(int queueIdx, List<Integer> compatibleTypes) {
            int[] sortedTypes = alfsOrder[queueIdx];
            for (int i = 0; i < sortedTypes.length; i++) {
                int typeId = sortedTypes[i];
                if (compatibleTypes.contains(typeId)) {
                    return typeId;
                }
            }
            return compatibleTypes.isEmpty() ? -1 : compatibleTypes.get(0);
        }

        /** FAIRNESS: Simple round-robin among compatible types. Kotlin line 12138. */
        private int selectServerTypeFAIRNESS(int queueIdx, List<Integer> compatibleTypes) {
            return selectServerTypeALIS(queueIdx, compatibleTypes);
        }

        /** FSF: Fastest Server First — pick type with highest service rate for this class. Kotlin line 12144. */
        private int selectServerTypeFSF(int queueIdx, int classId, List<Integer> compatibleTypes) {
            int bestType = -1;
            double bestRate = -1.0;
            for (int i = 0; i < compatibleTypes.size(); i++) {
                int typeId = compatibleTypes.get(i);
                double rate = heteroMus[queueIdx][typeId][classId];
                if (rate > bestRate && rate < Double.MAX_VALUE) {
                    bestRate = rate;
                    bestType = typeId;
                }
            }
            if (bestType >= 0) return bestType;
            return compatibleTypes.isEmpty() ? -1 : compatibleTypes.get(0);
        }

        /** RAIS: Random Available Idle Server. Kotlin line 12160. */
        private int selectServerTypeRAIS(List<Integer> compatibleTypes) {
            if (compatibleTypes.isEmpty()) return -1;
            int idx = (int) (routingRng.nextDouble() * compatibleTypes.size());
            if (idx < 0) idx = 0;
            if (idx > compatibleTypes.size() - 1) idx = compatibleTypes.size() - 1;
            return compatibleTypes.get(idx);
        }

        /**
         * Generates a setup (cold start) time for a server at the given queue and class.
         * Kotlin line 12473.
         */
        private double generateSetupTime(int queueIdx, int classId) {
            RandomVariateGen gen = setupGens[queueIdx][classId];
            return (gen != null) ? gen.nextDouble() : 0.0;
        }

        /**
         * Generates a delayoff (teardown) time for a server at the given queue and class.
         * Kotlin line 12485.
         */
        private double generateDelayoffTime(int queueIdx, int classId) {
            RandomVariateGen gen = delayoffGens[queueIdx][classId];
            return (gen != null) ? gen.nextDouble() : 0.0;
        }

        /** Holder for a transition-mode/firing-time race result. */
        private static final class ModeRaceResult {
            final TransitionModeInfo winningMode;
            final double winningTime;
            ModeRaceResult(TransitionModeInfo winningMode, double winningTime) {
                this.winningMode = winningMode;
                this.winningTime = winningTime;
            }
        }

        /**
         * Race semantics: pick winner among enabled modes.
         * Kotlin line 12740.
         */
        private ModeRaceResult selectModeByRace(int transListIdx, List<TransitionModeInfo> modes) {
            if (modes.isEmpty()) return new ModeRaceResult(null, 0.0);
            if (modes.size() == 1) {
                TransitionModeInfo mode = modes.get(0);
                double time = (mode.timingStrategy == TimingStrategy.IMMEDIATE)
                        ? 0.0
                        : sampleTransitionFiringTime(transListIdx, mode.modeIdx);
                return new ModeRaceResult(mode, time);
            }

            double minTime = Double.MAX_VALUE;
            TransitionModeInfo winningMode = null;

            for (int i = 0; i < modes.size(); i++) {
                TransitionModeInfo mode = modes.get(i);
                double firingTime;
                if (mode.timingStrategy == TimingStrategy.IMMEDIATE) {
                    // Immediate transitions fire at time 0, use weight for tie-breaking
                    firingTime = -mode.weight;  // More negative = higher priority (larger weight wins)
                } else {
                    // Timed transition: sample the firing time
                    firingTime = sampleTransitionFiringTime(transListIdx, mode.modeIdx);
                }

                if (firingTime < minTime) {
                    minTime = firingTime;
                    winningMode = mode;
                }
            }

            // For immediate transitions with negative times, use 0 for actual firing
            double actualTime = (minTime < 0) ? 0.0 : minTime;
            return new ModeRaceResult(winningMode, actualTime);
        }

        /**
         * Selects a transition mode by priority and probabilistic weight tie-break.
         * Kotlin line 12772.
         */
        private TransitionModeInfo selectTransitionMode(List<TransitionModeInfo> enabledModes) {
            if (enabledModes.isEmpty()) return null;
            if (enabledModes.size() == 1) return enabledModes.get(0);

            // Find highest priority among enabled modes
            int maxPriority = Integer.MIN_VALUE;
            for (int i = 0; i < enabledModes.size(); i++) {
                int p = enabledModes.get(i).priority;
                if (p > maxPriority) maxPriority = p;
            }

            List<TransitionModeInfo> highPriorityModes = new ArrayList<TransitionModeInfo>();
            for (int i = 0; i < enabledModes.size(); i++) {
                if (enabledModes.get(i).priority == maxPriority) {
                    highPriorityModes.add(enabledModes.get(i));
                }
            }

            if (highPriorityModes.size() == 1) return highPriorityModes.get(0);

            // Select by weight (probabilistic) among equally-prioritized modes
            double totalWeight = 0.0;
            for (int i = 0; i < highPriorityModes.size(); i++) {
                totalWeight += highPriorityModes.get(i).weight;
            }
            if (totalWeight <= 0) return highPriorityModes.get(0);

            double rand = siroRng.nextDouble() * totalWeight;
            double cumWeight = 0.0;
            for (int i = 0; i < highPriorityModes.size(); i++) {
                cumWeight += highPriorityModes.get(i).weight;
                if (rand <= cumWeight) return highPriorityModes.get(i);
            }
            return highPriorityModes.get(highPriorityModes.size() - 1);
        }

        /**
         * Check if a transition mode is enabled (all enabling conditions met,
         * no inhibiting conditions violated, and server availability).
         * Kotlin line 12799.
         */
        private boolean isTransitionModeEnabled(int transListIdx, TransitionModeInfo mode) {
            // Check server availability
            if (transitionInService[transListIdx][mode.modeIdx] >= mode.numServers) {
                return false;
            }

            // Check enabling conditions: each input place must have required available tokens
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int required = mode.enablingConditions[placeListIdx][classIdx];
                    int available = placeAvail(placeListIdx, classIdx);
                    if (required > 0 && available < required) {
                        return false;
                    }
                }
            }

            // Check inhibiting conditions: each input place must have fewer tokens than threshold
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int maxTokens = mode.inhibitingConditions[placeListIdx][classIdx];
                    int currentTokens = placeAvail(placeListIdx, classIdx);
                    if (maxTokens < Integer.MAX_VALUE && currentTokens >= maxTokens) {
                        return false;
                    }
                }
            }

            return true;
        }

        /**
         * Fire a transition: consume input tokens, track in transit, schedule completion.
         * Kotlin line 12833.
         * preSampledFiringTime may be null (use NaN sentinel via Double object).
         */
        private void fireTransition(int transListIdx, TransitionModeInfo mode,
                                    Double preSampledFiringTime) {
            double currentTime = ssjSim.time();

            // Update time-weighted token counts and consume input tokens from places
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                // Update statistics before consuming tokens
                double elapsed = currentTime - lastPlaceUpdateTime[placeListIdx];
                if (elapsed > 0) {
                    for (int k = 0; k < numClasses; k++) {
                        totalPlaceTokenTime[placeListIdx][k] +=
                                placeTokens[placeListIdx][k] * elapsed;
                    }
                    lastPlaceUpdateTime[placeListIdx] = currentTime;
                }

                // Consume tokens and track completions
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int required = mode.enablingConditions[placeListIdx][classIdx];
                    if (required > 0) {
                        consumePlaceTokens(placeListIdx, classIdx, required);
                        placeCompletions[placeListIdx][classIdx] += required;
                    }
                }
            }

            // Update tokens in transit statistics before adding new tokens
            double transitElapsed = currentTime - lastTransitUpdateTime[transListIdx];
            if (transitElapsed > 0) {
                for (int k = 0; k < numClasses; k++) {
                    totalTransitTokenTime[transListIdx][k] +=
                            tokensInTransit[transListIdx][k] * transitElapsed;
                }
                lastTransitUpdateTime[transListIdx] = currentTime;
            }

            // Add consumed tokens to "in transit" and track per-place transit
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                for (int classIdx = 0; classIdx < numClasses; classIdx++) {
                    int required = mode.enablingConditions[placeListIdx][classIdx];
                    if (required > 0) {
                        tokensInTransit[transListIdx][classIdx] += required;
                        // Update place transit statistics before adding tokens
                        double placeTransitElapsed =
                                currentTime - lastPlaceTransitUpdateTime[placeListIdx];
                        if (placeTransitElapsed > 0) {
                            for (int k = 0; k < numClasses; k++) {
                                placeTransitTokenTime[placeListIdx][k] +=
                                        placeTokensInTransit[placeListIdx][k] * placeTransitElapsed;
                            }
                            lastPlaceTransitUpdateTime[placeListIdx] = currentTime;
                        }
                        // Track that these tokens from this place are now in transit
                        placeTokensInTransit[placeListIdx][classIdx] += required;
                    }
                }
            }

            // Increment in-service count
            transitionInService[transListIdx][mode.modeIdx]++;

            // Schedule firing completion
            if (mode.timingStrategy == TimingStrategy.IMMEDIATE) {
                // Immediate transition: fire at current time (schedule with 0 delay)
                new TransitionFiring(transListIdx, mode.modeIdx).schedule(0.0);
            } else {
                // Timed transition: use pre-sampled time if available, otherwise sample new
                double firingTime = (preSampledFiringTime != null)
                        ? preSampledFiringTime.doubleValue()
                        : sampleTransitionFiringTime(transListIdx, mode.modeIdx);
                new TransitionFiring(transListIdx, mode.modeIdx).schedule(firingTime);
            }
        }

        /**
         * Sample firing time for a timed transition.
         * Kotlin line 12903.
         */
        /**
         * Cancel every in-flight dependent firing clock and free its server, so the
         * timed top-up in checkAndFireTransitions redraws it at the current marking.
         * Exact for exponential firing by the memoryless property (the residual of a
         * canceled exponential redrawn at the new rate is the correct competing clock).
         */
        private void resampleDependentFirings() {
            if (inflightDependent.isEmpty()) return;
            List<TransitionFiring> pending = new ArrayList<TransitionFiring>(inflightDependent);
            inflightDependent.clear();
            for (int i = 0; i < pending.size(); i++) {
                TransitionFiring tf = pending.get(i);
                tf.cancel();
                if (transitionInService[tf.transListIdx][tf.modeIdx] > 0) {
                    transitionInService[tf.transListIdx][tf.modeIdx]--;
                }
            }
        }

        private double sampleTransitionFiringTime(int transListIdx, int modeIdx) {
            RandomVariateGen gen = transitionFiringGens[transListIdx][modeIdx];
            double delay;
            if (gen != null) {
                delay = gen.nextDouble();
            } else {
                // Fallback: exponential with rate 1
                delay = -FastMath.log(siroRng.nextDouble());
            }
            // Marking-dependent firing-rate multiplier g(marking): effective rate is
            // rate_base*g, so the memoryless delay is divided by g. Sampled at the
            // current marking and resampled on every marking change (see
            // resampleDependentFirings), which is exact for exponential firing.
            SerializableFunction<Matrix, Double> g = transitionModes[transListIdx].get(modeIdx).firingDep;
            if (g != null) {
                double mult = g.apply(currentMarkingMatrix());
                if (mult > 0) {
                    delay = delay / mult;
                } else {
                    delay = Double.POSITIVE_INFINITY; // zero rate: never fires at this marking
                }
            }
            return delay;
        }

        /**
         * Node-indexed marking matrix (numNodes x numClasses) built from the
         * current place tokens, the argument passed to a firing-rate dependence
         * handle. Only place rows are nonzero, matching the JSON lattice the
         * handle was tabulated over.
         */
        private Matrix currentMarkingMatrix() {
            Matrix m = new Matrix(numNodes, numClasses);
            m.zero();
            for (int placeListIdx = 0; placeListIdx < placeNodes.size(); placeListIdx++) {
                int nodeIdx = placeNodes.get(placeListIdx);
                for (int c = 0; c < numClasses; c++) {
                    m.set(nodeIdx, c, placeTokens[placeListIdx][c]);
                }
            }
            return m;
        }

        /**
         * Transition firing completion event. When the firing delay elapses,
         * tokens are removed from transit and produced at output nodes.
         * Kotlin line 12917.
         */
        private final class TransitionFiring extends SimEvent {
            @Override
            protected double slotPhase() {
                return SLOT_PHASE_INTERNAL;
            }

            private final int transListIdx;
            private final int modeIdx;

            TransitionFiring(int transListIdx, int modeIdx) {
                this.transListIdx = transListIdx;
                this.modeIdx = modeIdx;
            }

            @Override
            public void actions() {
                TransitionModeInfo mode = transitionModes[transListIdx].get(modeIdx);

                // This clock fired, so it is no longer a resampleable in-flight clock.
                if (mode.firingDep != null) inflightDependent.remove(this);
                // This clock instance is no longer in flight.
                transitionInService[transListIdx][modeIdx]--;

                // Atomic-firing race semantics: the required tokens were not reserved
                // when this clock was scheduled, so re-check the mode is still enabled
                // in the current marking. If a competing firing consumed the tokens in
                // the meantime the clock is simply discarded (enabling-memory policy;
                // exact for exponential firing times), and checkAndFireTransitions will
                // reschedule the mode once it becomes enabled again.
                if (transitionModeEnablingDegree(mode) >= 1) {
                    fireAtomic(transListIdx, mode);
                }

                // Re-evaluate the marking: fire immediate transitions and top up
                // timed clocks for whatever is now enabled.
                checkAndFireTransitions();
            }
        }

        /**
         * Route a token from a transition to a destination node.
         * Kotlin line 12983.
         */
        private void routeTokenFromTransition(int transNodeIdx, int destNodeIdx, int classId) {
            double currentTime = ssjSim.time();

            if (placeNodes.contains(destNodeIdx)) {
                // Destination is a Place - add token with proper statistics tracking
                int placeListIdx = placeNodes.indexOf(destNodeIdx);
                if (placeListIdx >= 0) {
                    // Update time-weighted token count before adding new token
                    double elapsed = currentTime - lastPlaceUpdateTime[placeListIdx];
                    if (elapsed > 0) {
                        for (int k = 0; k < numClasses; k++) {
                            totalPlaceTokenTime[placeListIdx][k] +=
                                    placeTokens[placeListIdx][k] * elapsed;
                        }
                        lastPlaceUpdateTime[placeListIdx] = currentTime;
                    }
                    // Add token to place (marking). For a queueing place the token enters the
                    // embedded queue and becomes available to output transitions only after
                    // service completion; for an ordinary place it is immediately available.
                    placeTokens[placeListIdx][classId]++;
                    if (isQueueingPlace[placeListIdx]) {
                        placeWaiting[placeListIdx].addLast(Integer.valueOf(classId));
                        tryStartPlaceService(placeListIdx);
                    }
                }
                // Note: checkAndFireTransitions() is called at the end of TransitionFiring
            } else if (sinkNodes.contains(destNodeIdx)) {
                // Token leaves system
                systemCompletedCustomers[classId]++;
            } else if (serviceNodes.contains(destNodeIdx)) {
                // Route to service node (queue)
                int queueIdx = serviceNodes.indexOf(destNodeIdx);
                Customer customer = new Customer(
                        classId, classPrio[classId], currentTime, currentTime,
                        siroRng.nextDouble(),
                        -1.0, -1L,
                        currentTime + classDeadline[classId], -1, null);
                arriveAtQueue(queueIdx, customer);
            } else if (forkNodes.contains(destNodeIdx)) {
                // Route to fork
                long parentJobId = nextJobId++;
                handleForkArrival(destNodeIdx, parentJobId, classId, currentTime);
            } else if (transitionNodes.contains(destNodeIdx)) {
                // Token goes directly to another transition
                checkAndFireTransitions();
            }
        }

        /** Check if a node is a Place node. Kotlin line 13030. */
        private boolean isPlaceNode(int nodeIdx) {
            return placeNodes.contains(nodeIdx);
        }

        /** Check if a node is a Transition node. Kotlin line 13035. */
        private boolean isTransitionNode(int nodeIdx) {
            return transitionNodes.contains(nodeIdx);
        }

        /**
         * Gets the number of servers compatible with a given class at a queue.
         * For heterogeneous queues, returns sum of servers across compatible types.
         * Kotlin line 13070.
         */
        /**
         * Peak service-rate scaling at a station, in units of its nominal rate, or
         * 0 when the station has no load or class dependence (the caller then falls
         * back to the measured busy-time tally).
         *
         * <p>For load dependence this is max_n alpha(n) over the tabulated range.
         * For class dependence it is max over the reachable per-class population
         * lattice and over the classes, since utilization is a per-station quantity
         * and the whole station shares one normalizer. Both mirror
         * max(lldscaling(ist,:)) in solver_ncld.</p>
         */
        private double getPeakScaling(int queueIdx) {
            double peak = 0.0;
            if (hasCd && cdFunctions != null && queueIdx < cdFunctions.length
                    && cdFunctions[queueIdx] != null) {
                int K = numClasses;
                int[] maxc = new int[K];
                for (int r = 0; r < K; r++) {
                    double nj = sn.njobs.get(r);
                    // Open classes are unbounded; the tabulated handle saturates, so
                    // probe the same cutoff the JSON writer materializes them on.
                    maxc[r] = (Double.isFinite(nj) && nj > 0) ? (int) Math.round(nj) : 10;
                }
                int[] n = new int[K];
                while (true) {
                    int tot = 0;
                    for (int r = 0; r < K; r++) tot += n[r];
                    if (tot > 0) {
                        Matrix nv = new Matrix(1, K);
                        for (int r = 0; r < K; r++) nv.set(0, r, n[r]);
                        Matrix bv = cdFunctions[queueIdx].apply(nv);
                        for (int r = 0; r < bv.getNumElements(); r++) {
                            double v = bv.get(r);
                            if (Double.isFinite(v) && v > peak) peak = v;
                        }
                    }
                    int s = K - 1;
                    while (s >= 0 && n[s] == maxc[s]) { n[s] = 0; s--; }
                    if (s < 0) break;
                    n[s]++;
                }
                return peak;
            }
            if (isLoadDependent != null && queueIdx < isLoadDependent.length
                    && isLoadDependent[queueIdx] && lldScaling != null
                    && lldScaling[queueIdx] != null) {
                for (int i = 0; i < lldScaling[queueIdx].length; i++) {
                    double v = lldScaling[queueIdx][i];
                    if (Double.isFinite(v) && v > peak) peak = v;
                }
                return peak;
            }
            return 0.0;
        }

        private int getCompatibleServerCount(int queueIdx, int classId) {
            if (numServerTypes[queueIdx] > 0) {
                int count = 0;
                for (int typeId = 0; typeId < numServerTypes[queueIdx]; typeId++) {
                    if (serverCompat[queueIdx][typeId][classId]) {
                        count += serversPerType[queueIdx][typeId];
                    }
                }
                return (count > 0) ? count : numServers[queueIdx];
            }
            return numServers[queueIdx];
        }

        /**
         * Get effective utilization (service time only, excluding blocking time).
         * Kotlin line 13088.
         */
        public double getEffectiveUtilization(int queueIdx, int classId) {
            if (isDelayNode.get(queueIdx)) {
                double simTime = getActualSimTime();
                if (simTime > 0) {
                    return totalQueueTime[queueIdx][classId] / simTime;
                }
                return 0.0;
            }
            int compatibleServers = getCompatibleServerCount(queueIdx, classId);
            // Use MSER-5 truncated busy time if available
            if (mserEnabled && busyTimeObservations != null && observationTimes != null
                    && queueIdx < busyTimeObservations.length
                    && classId < busyTimeObservations[queueIdx].length) {
                List<Double> busyObs = busyTimeObservations[queueIdx][classId];
                int truncationIdx = mserTruncationBatch * effectiveMserBatchSize;
                if (truncationIdx < busyObs.size() && truncationIdx < observationTimes.size()) {
                    double startBusy = busyObs.get(truncationIdx);
                    double endBusy = totalBusyTime[queueIdx][classId];
                    double startTime = observationTimes.get(truncationIdx);
                    double endTime = ssjSim.time();
                    double elapsed = endTime - startTime;
                    if (elapsed > 0 && compatibleServers > 0) {
                        return (endBusy - startBusy) / (elapsed * compatibleServers);
                    }
                }
            }
            // Fallback
            double simTime = getActualSimTime();
            if (simTime > 0 && compatibleServers > 0) {
                return totalBusyTime[queueIdx][classId] / (simTime * compatibleServers);
            }
            return 0.0;
        }

        /** Returns arrival rate at a service node for a given class. Kotlin line 13225. */
        public double getArrivalRate(int queueIdx, int classId) {
            double simTime = getActualSimTime();
            if (simTime > 0) {
                return ((double) arrivedCustomers[queueIdx][classId]) / simTime;
            }
            return 0.0;
        }

        // ==================== OBM Confidence Interval Functions ====================

        /**
         * Compute overlap adjustment factor for OBM variance estimation.
         * Kotlin line 13273.
         */
        private double computeOverlapAdjustmentFactor(double overlap) {
            if (overlap <= 0.0) {
                return 1.0;  // Non-overlapping
            } else if (overlap >= 0.5) {
                return 4.0 / 3.0;  // 50% overlap (standard OBM)
            } else {
                return 1.0 + (overlap / 0.5) * (4.0 / 3.0 - 1.0);
            }
        }

        /** Triple of (grandMean, stdError, df) returned by batch-means routines. */
        private static final class StatTriple {
            final double grandMean;
            final double stdError;
            final int df;
            StatTriple(double grandMean, double stdError, int df) {
                this.grandMean = grandMean;
                this.stdError = stdError;
                this.df = df;
            }
        }

        /**
         * Compute OBM (Overlapping Batch Means) statistics with configurable overlap.
         * Kotlin line 13287.
         */
        private StatTriple computeOBMStatisticsInternal(List<Double> observations, int batchSize) {
            int n = observations.size();
            if (n < batchSize * 2) return null;

            // Use configurable overlap fraction
            int stepSize = (int) (batchSize * (1.0 - effectiveObmOverlap));
            if (stepSize < 1) stepSize = 1;
            int numBatches = (n - batchSize) / stepSize + 1;
            if (numBatches < 2) return null;

            // Compute overlapping batch means
            double[] batchMeans = new double[numBatches];
            for (int i = 0; i < numBatches; i++) {
                int startIdx = i * stepSize;
                double sum = 0.0;
                for (int j = 0; j < batchSize; j++) {
                    if (startIdx + j < n) {
                        sum += observations.get(startIdx + j);
                    }
                }
                batchMeans[i] = sum / batchSize;
            }

            // Compute grand mean
            double grandMean = 0.0;
            for (int i = 0; i < numBatches; i++) {
                grandMean += batchMeans[i];
            }
            grandMean /= numBatches;

            // Compute sum of squared differences
            double sumSquaredDiff = 0.0;
            for (int i = 0; i < numBatches; i++) {
                double diff = batchMeans[i] - grandMean;
                sumSquaredDiff += diff * diff;
            }

            // OBM variance with configurable overlap adjustment factor
            double overlapAdjustment = computeOverlapAdjustmentFactor(effectiveObmOverlap);
            double batchMeanVariance = overlapAdjustment * sumSquaredDiff / (numBatches - 1);
            double stdError = Math.sqrt(batchMeanVariance / numBatches);

            // Effective degrees of freedom for OBM
            int effectiveDf = (int) ((numBatches - 1) / overlapAdjustment);
            if (effectiveDf < 1) effectiveDf = 1;

            return new StatTriple(grandMean, stdError, effectiveDf);
        }

        /**
         * Compute BM (Non-overlapping Batch Means) statistics.
         * Kotlin line 13340.
         */
        private StatTriple computeBMStatisticsInternal(List<Double> observations, int batchSize) {
            int n = observations.size();
            int numBatches = n / batchSize;
            if (numBatches < 2) return null;

            // Compute non-overlapping batch means
            double[] batchMeans = new double[numBatches];
            for (int i = 0; i < numBatches; i++) {
                int startIdx = i * batchSize;
                double sum = 0.0;
                for (int j = 0; j < batchSize; j++) {
                    sum += observations.get(startIdx + j);
                }
                batchMeans[i] = sum / batchSize;
            }

            // Compute grand mean
            double grandMean = 0.0;
            for (int i = 0; i < numBatches; i++) {
                grandMean += batchMeans[i];
            }
            grandMean /= numBatches;

            // Sample variance of batch means (no overlap adjustment needed)
            double sumSquaredDiff = 0.0;
            for (int i = 0; i < numBatches; i++) {
                double diff = batchMeans[i] - grandMean;
                sumSquaredDiff += diff * diff;
            }
            double batchMeanVariance = sumSquaredDiff / (numBatches - 1);
            double stdError = Math.sqrt(batchMeanVariance / numBatches);

            int df = numBatches - 1;
            return new StatTriple(grandMean, stdError, df);
        }

        /**
         * Compute Heidelberger-Welch spectral analysis statistics.
         * Kotlin line 13391.
         */
        private StatTriple computeSpectralStatisticsInternal(List<Double> observations, int batchSize) {
            int n = observations.size();
            int numBatches = n / batchSize;
            if (numBatches < 4) return computeBMStatisticsInternal(observations, batchSize);

            // Step 1: Compute non-overlapping batch means
            double[] batchMeans = new double[numBatches];
            for (int i = 0; i < numBatches; i++) {
                int startIdx = i * batchSize;
                double sum = 0.0;
                for (int j = 0; j < batchSize; j++) {
                    sum += observations.get(startIdx + j);
                }
                batchMeans[i] = sum / batchSize;
            }

            // Compute grand mean
            double grandMean = 0.0;
            for (int i = 0; i < numBatches; i++) {
                grandMean += batchMeans[i];
            }
            grandMean /= numBatches;

            // Step 2: Compute periodogram via direct DFT (O(m^2), m is small)
            int numFreqs = numBatches / 2;
            if (numFreqs < 2) return computeBMStatisticsInternal(observations, batchSize);

            double[] periodogram = new double[numFreqs];
            for (int k = 1; k <= numFreqs; k++) {
                double realPart = 0.0;
                double imagPart = 0.0;
                double freq = 2.0 * Math.PI * k / numBatches;
                for (int t = 0; t < numBatches; t++) {
                    double centered = batchMeans[t] - grandMean;
                    realPart += centered * Math.cos(freq * t);
                    imagPart += centered * Math.sin(freq * t);
                }
                periodogram[k - 1] = (realPart * realPart + imagPart * imagPart) / numBatches;
            }

            // Step 3: Log-periodogram regression on lowest fraction of frequencies
            int numRegPts = (int) (numFreqs * effectiveSpectralLowFreqFrac);
            if (numRegPts < 3) numRegPts = 3;
            if (numRegPts < 3) return computeBMStatisticsInternal(observations, batchSize);

            // Check for zero periodogram ordinates (would make log undefined)
            for (int k = 0; k < numRegPts; k++) {
                if (periodogram[k] <= 0.0) {
                    return computeBMStatisticsInternal(observations, batchSize);
                }
            }

            // Fit quadratic: log(I(f_k)) = beta_0 + beta_1*f_k + beta_2*f_k^2
            // Using Apache Commons Math OLS regression
            org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression ols =
                    new org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression();
            double[] yData = new double[numRegPts];
            double[][] xData = new double[numRegPts][2];
            for (int k = 0; k < numRegPts; k++) {
                double fk = 2.0 * Math.PI * (k + 1) / numBatches;
                yData[k] = Math.log(periodogram[k]);
                xData[k][0] = fk;
                xData[k][1] = fk * fk;
            }

            double[] beta;
            try {
                ols.newSampleData(yData, xData);
                beta = ols.estimateRegressionParameters();  // [intercept, beta1, beta2]
            } catch (Exception e) {
                // Regression failed (e.g., singular matrix) — fall back to BM
                return computeBMStatisticsInternal(observations, batchSize);
            }

            // Step 4: S(0) = exp(beta_0)  (intercept is the log-spectral-density at f=0)
            double s0 = Math.exp(beta[0]);
            if (Double.isNaN(s0) || Double.isInfinite(s0) || s0 <= 0.0) {
                return computeBMStatisticsInternal(observations, batchSize);
            }

            // Step 5: Variance of sample mean = 2*pi*S(0)/m
            double varMean = 2.0 * Math.PI * s0 / numBatches;
            double stdError = Math.sqrt(varMean);

            // Step 6: Degrees of freedom = numRegPts - 3 (3 regression parameters)
            int df = numRegPts - 3;
            if (df < 1) df = 1;

            return new StatTriple(grandMean, stdError, df);
        }

        /**
         * Compute batch means statistics using the configured CI method.
         * Kotlin line 13482.
         */
        private StatTriple computeCIStatistics(List<Double> observations, int batchSize) {
            if ("obm".equals(effectiveCiMethod)) {
                return computeOBMStatisticsInternal(observations, batchSize);
            } else if ("bm".equals(effectiveCiMethod)) {
                return computeBMStatisticsInternal(observations, batchSize);
            } else if ("spectral".equals(effectiveCiMethod)) {
                return computeSpectralStatisticsInternal(observations, batchSize);
            }
            return null;  // "none" - no CI computation
        }

    }

    // ==================== Top-Level OBM Test Helper Functions ====================

    /**
     * Compute OBM (Overlapping Batch Means) statistics with 50% overlap.
     * Exposed as static helper for testing.
     * @param observations List of observations to analyze
     * @param batchSize Size of each batch
     * @return jline.util.Triple(grandMean, stdError, effectiveDf) or null if insufficient data
     */
    public static org.apache.commons.lang3.tuple.Triple<Double, Double, Integer> computeOBMStatistics(
            List<Double> observations, int batchSize) {
        int n = observations.size();
        if (n < batchSize * 2) return null;

        int halfBatch = batchSize / 2;
        int numBatches = (n - batchSize) / halfBatch + 1;
        if (numBatches < 2) return null;

        // Compute overlapping batch means
        double[] batchMeans = new double[numBatches];
        for (int i = 0; i < numBatches; i++) {
            int startIdx = i * halfBatch;
            double sum = 0.0;
            for (int j = 0; j < batchSize; j++) {
                sum += observations.get(startIdx + j);
            }
            batchMeans[i] = sum / batchSize;
        }

        // Compute grand mean
        double grandMean = 0.0;
        for (int i = 0; i < numBatches; i++) {
            grandMean += batchMeans[i];
        }
        grandMean /= numBatches;

        // Compute sum of squared differences
        double sumSquaredDiff = 0.0;
        for (int i = 0; i < numBatches; i++) {
            double diff = batchMeans[i] - grandMean;
            sumSquaredDiff += diff * diff;
        }

        // OBM variance with 4/3 overlap adjustment factor
        double overlapAdjustment = 4.0 / 3.0;
        double batchMeanVariance = overlapAdjustment * sumSquaredDiff / (numBatches - 1);
        double stdError = Math.sqrt(batchMeanVariance / numBatches);

        // Effective degrees of freedom for OBM
        int effectiveDf = Math.max(1, (int) ((numBatches - 1) / overlapAdjustment));

        return org.apache.commons.lang3.tuple.Triple.of(grandMean, stdError, effectiveDf);
    }

    /**
     * Get t-distribution critical value for given confidence level and degrees of freedom.
     * Exposed as static helper for testing.
     * @param confintLevel Confidence level (e.g., 0.95 for 95% CI)
     * @param df Degrees of freedom
     * @return Critical value from t-distribution
     */
    public static double getTCriticalValue(double confintLevel, int df) {
        // Lookup based on confidence level directly for better precision
        double[] tTable95 = new double[]{
            12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
            2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
            2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042
        };
        double[] tTable90 = new double[]{
            6.314, 2.920, 2.353, 2.132, 2.015, 1.943, 1.895, 1.860, 1.833, 1.812,
            1.796, 1.782, 1.771, 1.761, 1.753, 1.746, 1.740, 1.734, 1.729, 1.725,
            1.721, 1.717, 1.714, 1.711, 1.708, 1.706, 1.703, 1.701, 1.699, 1.697
        };
        double[] tTable99 = new double[]{
            63.657, 9.925, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.250, 3.169,
            3.106, 3.055, 3.012, 2.977, 2.947, 2.921, 2.898, 2.878, 2.861, 2.845,
            2.831, 2.819, 2.807, 2.797, 2.787, 2.779, 2.771, 2.763, 2.756, 2.750
        };

        // Table index is df - 1 (df=1 is at index 0)
        int tableIdx = Math.min(df, 30) - 1;
        if (tableIdx < 0) return 1.96;

        // Select table based on confidence level (checking higher confidence first)
        if (confintLevel >= 0.99) {
            return tableIdx < tTable99.length ? tTable99[tableIdx] : 2.576;
        } else if (confintLevel >= 0.95) {
            return tableIdx < tTable95.length ? tTable95[tableIdx] : 1.96;
        } else if (confintLevel >= 0.90) {
            return tableIdx < tTable90.length ? tTable90[tableIdx] : 1.645;
        } else {
            return 1.96;
        }
    }

    /**
     * Compute standard (non-overlapping) batch means statistics.
     * Exposed as static helper for testing comparison with OBM.
     * @param observations List of observations to analyze
     * @param batchSize Size of each batch
     * @return jline.util.Triple(grandMean, stdError, effectiveDf) or null if insufficient data
     */
    public static org.apache.commons.lang3.tuple.Triple<Double, Double, Integer> computeStandardBatchMeansStatistics(
            List<Double> observations, int batchSize) {
        int n = observations.size();
        int numBatches = n / batchSize;
        if (numBatches < 2) return null;

        // Compute non-overlapping batch means
        double[] batchMeans = new double[numBatches];
        for (int i = 0; i < numBatches; i++) {
            int startIdx = i * batchSize;
            double sum = 0.0;
            for (int j = 0; j < batchSize; j++) {
                sum += observations.get(startIdx + j);
            }
            batchMeans[i] = sum / batchSize;
        }

        // Compute grand mean
        double grandMean = 0.0;
        for (int i = 0; i < numBatches; i++) {
            grandMean += batchMeans[i];
        }
        grandMean /= numBatches;

        // Compute sum of squared differences
        double sumSquaredDiff = 0.0;
        for (int i = 0; i < numBatches; i++) {
            double diff = batchMeans[i] - grandMean;
            sumSquaredDiff += diff * diff;
        }

        // Standard batch means variance (no overlap adjustment)
        double batchMeanVariance = sumSquaredDiff / (numBatches - 1);
        double stdError = Math.sqrt(batchMeanVariance / numBatches);

        return org.apache.commons.lang3.tuple.Triple.of(grandMean, stdError, numBatches - 1);
    }

    /**
     * Compute Heidelberger-Welch spectral analysis statistics.
     * Exposed as static helper for testing.
     *
     * Uses log-periodogram regression at low frequencies to estimate the spectral
     * density at frequency zero, which accounts for autocorrelation between batches.
     *
     * @param observations List of observations to analyze
     * @param batchSize Size of each batch
     * @param lowFreqFrac Fraction of lowest frequencies to use for regression (default 0.25)
     * @return jline.util.Triple(grandMean, stdError, df) or null if insufficient data or regression fails
     */
    public static org.apache.commons.lang3.tuple.Triple<Double, Double, Integer> computeSpectralStatistics(
            List<Double> observations, int batchSize, double lowFreqFrac) {
        int n = observations.size();
        int numBatches = n / batchSize;
        if (numBatches < 4) return null;

        // Step 1: Compute non-overlapping batch means
        double[] batchMeans = new double[numBatches];
        for (int i = 0; i < numBatches; i++) {
            int startIdx = i * batchSize;
            double sum = 0.0;
            for (int j = 0; j < batchSize; j++) {
                sum += observations.get(startIdx + j);
            }
            batchMeans[i] = sum / batchSize;
        }

        // Compute grand mean
        double grandMean = 0.0;
        for (int i = 0; i < batchMeans.length; i++) {
            grandMean += batchMeans[i];
        }
        grandMean /= numBatches;

        // Step 2: Compute periodogram via direct DFT
        int numFreqs = numBatches / 2;
        if (numFreqs < 2) return null;

        double[] periodogram = new double[numFreqs];
        for (int k = 1; k <= numFreqs; k++) {
            double realPart = 0.0;
            double imagPart = 0.0;
            double freq = 2.0 * Math.PI * k / numBatches;
            for (int t = 0; t < numBatches; t++) {
                double centered = batchMeans[t] - grandMean;
                realPart += centered * Math.cos(freq * t);
                imagPart += centered * Math.sin(freq * t);
            }
            periodogram[k - 1] = (realPart * realPart + imagPart * imagPart) / numBatches;
        }

        // Step 3: Log-periodogram regression on lowest frequencies
        int numRegPts = Math.max(3, (int) (numFreqs * lowFreqFrac));
        if (numRegPts < 3) return null;

        // Check for zero periodogram ordinates
        for (int k = 0; k < numRegPts; k++) {
            if (periodogram[k] <= 0.0) return null;
        }

        // Fit quadratic: log(I(f_k)) = beta_0 + beta_1*f_k + beta_2*f_k^2
        org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression ols =
                new org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression();
        double[] yData = new double[numRegPts];
        double[][] xData = new double[numRegPts][2];
        for (int k = 0; k < numRegPts; k++) {
            double fk = 2.0 * Math.PI * (k + 1) / numBatches;
            yData[k] = Math.log(periodogram[k]);
            xData[k][0] = fk;
            xData[k][1] = fk * fk;
        }

        double[] beta;
        try {
            ols.newSampleData(yData, xData);
            beta = ols.estimateRegressionParameters();
        } catch (Exception e) {
            return null;
        }

        // S(0) = exp(beta_0)
        double s0 = Math.exp(beta[0]);
        if (Double.isNaN(s0) || Double.isInfinite(s0) || s0 <= 0.0) return null;

        // Variance of sample mean = 2*pi*S(0)/m
        double varMean = 2.0 * Math.PI * s0 / numBatches;
        double stdError = Math.sqrt(varMean);

        // Degrees of freedom = numRegPts - 3
        int df = Math.max(1, numRegPts - 3);

        return org.apache.commons.lang3.tuple.Triple.of(grandMean, stdError, df);
    }

    /**
     * Overload for default lowFreqFrac = 0.25.
     */
    public static org.apache.commons.lang3.tuple.Triple<Double, Double, Integer> computeSpectralStatistics(
            List<Double> observations, int batchSize) {
        return computeSpectralStatistics(observations, batchSize, 0.25);
    }
}
