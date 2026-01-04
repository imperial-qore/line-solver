/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * @file Solver_ssj.kt
 * @brief SSJ-based discrete event simulation engine for queueing networks and Petri nets.
 *
 * @details This file contains the core DES simulation engine built on the SSJ
 * (Stochastic Simulation in Java) library. It provides both steady-state and
 * transient analysis capabilities for complex queueing network models.
 *
 * @section ssj_architecture Architecture Overview
 *
 * The simulation engine is organized around the SSJSimulator class which manages:
 * - Event scheduling and processing using SSJ's event-driven framework
 * - Customer/token lifecycle from arrival to departure
 * - Queue management with multiple scheduling disciplines
 * - Statistics collection for performance metrics
 *
 * @section ssj_events Event Types
 *
 * The simulator uses the following event types:
 * - **SourceArrival**: New job arrival from external source
 * - **ServiceCompletion**: Job finishes service at a queue
 * - **DelayCompletion**: Job finishes at infinite-server station
 * - **ForkArrival/JoinArrival**: Fork-join synchronization events
 * - **TransitionFiring**: Petri net transition fires
 * - **PSServiceCheck**: Processor sharing time slice events
 *
 * @section ssj_statistics Statistics Collection
 *
 * The simulator collects:
 * - Time-weighted queue lengths (for average queue length)
 * - Cumulative busy time (for utilization)
 * - Job counts and response times (for throughput)
 * - Per-class metrics for multiclass systems
 *
 * @section ssj_features Supported Features
 *
 * See SolverDES documentation for complete feature list.
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */

package jline.solvers.des.handlers

import jline.GlobalConstants
import jline.VerboseLevel
import jline.api.mam.bmap_sample
import jline.api.mam.mmap_sample
import jline.api.mam.map_sample
import jline.io.lineTempName
import jline.io.line_warning
import jline.lang.NetworkStruct
import jline.lang.constant.DropStrategy
import jline.lang.constant.HeteroSchedPolicy
import jline.lang.constant.JoinStrategy
import jline.lang.constant.NodeType
import jline.lang.constant.PollingType
import jline.lang.constant.RemovalPolicy
import jline.api.sn.snGetArvRFromTput
import jline.lang.constant.ProcessType
import jline.lang.constant.ReplacementStrategy
import jline.lang.constant.RoutingStrategy
import jline.lang.constant.SchedStrategy
import jline.lang.constant.SignalType
import jline.lang.constant.TimingStrategy
import jline.lang.nodes.Cache
import jline.lang.nodes.Queue
import jline.lang.nodes.Source
import jline.lang.nodeparam.CacheNodeParam
import jline.lang.nodeparam.QueueNodeParam
import jline.lang.nodeparam.TransitionNodeParam
import jline.lang.processes.DiscreteDistribution
import jline.lang.processes.Replayer
import jline.solvers.SolverOptions
import jline.solvers.des.DESOptions
import jline.solvers.des.DESResult
import jline.streaming.Collector
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import umontreal.ssj.randvar.*
import umontreal.ssj.rng.MRG32k3a
import umontreal.ssj.simevents.Event
import umontreal.ssj.simevents.Sim
import umontreal.ssj.stat.Tally
import java.io.File
import java.io.PrintWriter
import java.util.LinkedList
import java.util.PriorityQueue
import java.util.Random

/**
 * Custom random variate generator for deterministic (constant) values.
 * Used when simulating deterministic distributions in SSJ.
 */
class ConstantGen(stream: umontreal.ssj.rng.RandomStream, private val value: Double) :
    umontreal.ssj.randvar.RandomVariateGen(stream, null) {
    override fun nextDouble(): Double = value
}

/**
 * @brief Steady-state queueing network simulation using SSJ library.
 *
 * @details This function performs steady-state discrete event simulation of
 * queueing networks and stochastic Petri nets using the SSJ library.
 *
 * @section solver_ssj_model Supported Model Elements
 * - **Source nodes**: Poisson or general arrival processes
 * - **Queue nodes**: Multiple scheduling disciplines (FCFS, PS, LCFS, etc.)
 * - **Delay nodes**: Infinite server stations (M/M/∞)
 * - **Sink nodes**: Job departure points
 * - **Fork/Join nodes**: Parallel processing with synchronization
 * - **Router nodes**: Multiple routing strategies
 * - **ClassSwitch nodes**: Dynamic class switching
 * - **Place/Transition nodes**: Stochastic Petri nets
 *
 * @section solver_ssj_warmup Warmup and Estimation
 * The simulation uses event-count based stopping, running until the specified
 * number of service completion events (options.samples) is reached. By default,
 * 20% of events are used as warmup period before collecting statistics, with
 * MSER-5 enabled to automatically determine the optimal warmup truncation point.
 *
 * @section solver_ssj_metrics Collected Metrics
 * - Queue length: Time-weighted average number of jobs
 * - Utilization: Fraction of time servers are busy
 * - Response time: Average time from arrival to departure
 * - Throughput: Average job completion rate
 *
 * @param sn Network structure representing the queueing network model
 * @param options Solver configuration options (samples = max events, seed, verbose, etc.)
 * @return DESResult with steady-state performance metrics
 *
 * @see solver_ssj_transient For transient analysis
 * @see SSJSimulator Core simulation engine
 */
fun solver_ssj(sn: NetworkStruct, options: SolverOptions, stream: Collector? = null): DESResult {
    val simulator = SSJSimulator(sn, options, null, stream)
    val maxEvents = options.samples.toDouble()
    simulator.simulate(maxEvents)
    return simulator.getDESResult()
}

/**
 * @brief Transient analysis of queueing network using SSJ library.
 *
 * @details This function performs transient discrete event simulation,
 * collecting time-series metrics from time 0 without warmup period.
 * Useful for analyzing system behavior during startup or after perturbations.
 *
 * @section solver_ssj_transient_output Transient Output
 * The result contains time-indexed matrices:
 * - **QNt[i][k]**: Queue length at station i, class k over time
 * - **UNt[i][k]**: Utilization at station i, class k over time
 * - **TNt[i][k]**: Throughput at station i, class k over time
 * - **t**: Vector of time points for the measurements
 *
 * @section solver_ssj_transient_sampling Time Sampling
 * Metrics are sampled at regular intervals determined by the simulation.
 * The number of samples is proportional to the time horizon.
 *
 * @section solver_ssj_transient_init Initial Conditions
 * If options.init_sol is provided (non-null and non-empty), the simulation
 * will start with jobs distributed according to this matrix instead of
 * placing all closed class jobs at reference stations. The init_sol matrix
 * should be [1 x (M*K)] in row-major order where M = stations, K = classes.
 *
 * @param sn Network structure representing the queueing network model
 * @param options Solver configuration options (timespan defines the analysis interval)
 * @return DESResult with transient time-series metrics in QNt, UNt, TNt arrays
 *
 * @see solver_ssj For steady-state analysis
 */
fun solver_ssj_transient(sn: NetworkStruct, options: SolverOptions, stream: Collector? = null): DESResult {
    val initSol = options.init_sol
    val simulator = SSJSimulator(sn, options, initSol, stream)
    val timespan = options.timespan
    if (timespan != null && timespan.size >= 2) {
        val endTime = timespan[1]
        simulator.simulateTransient(endTime)
    }
    return simulator.getTransientDESResult()
}

/**
 * @brief Core discrete event simulation engine using SSJ library.
 *
 * @details This class implements the main simulation logic for queueing networks
 * and stochastic Petri nets. It manages the event calendar, customer lifecycle,
 * queue operations, and statistics collection.
 *
 * @section simulator_components Key Components
 *
 * @subsection simulator_nodes Node Management
 * The simulator classifies and tracks different node types:
 * - Source/Sink: External arrival/departure points
 * - Queue/Delay: Service stations with finite/infinite servers
 * - Fork/Join: Parallel processing synchronization
 * - Router/ClassSwitch: Routing and class transformation
 * - Place/Transition: Petri net elements
 *
 * @subsection simulator_queues Queue Management
 * Each service station maintains:
 * - Wait queue: Jobs waiting for service (ordered by scheduling discipline)
 * - Server pool: Available servers for processing
 * - Per-class statistics: Queue lengths, busy times, completions
 *
 * @subsection simulator_events Event Processing
 * Events are processed in chronological order:
 * 1. Arrival events add jobs to queues
 * 2. Service completion events route jobs to next destination
 * 3. PS events manage time-sliced processor sharing
 * 4. Transition events fire Petri net transitions
 *
 * @subsection simulator_stats Statistics Collection
 * Statistics are collected using time-weighted averaging:
 * - Queue length: Integral of queue size over time
 * - Utilization: Integral of busy servers over time
 * - Throughput: Count of completed jobs over time
 *
 * @section simulator_debug Debug Output
 * When VerboseLevel is DEBUG, detailed event traces are saved to a CSV file named
 * "ssj_trace_{timestamp}.csv" in the current working directory.
 *
 * @param sn Network structure containing model topology and parameters
 * @param options Solver options (samples, seed, verbose level, etc.)
 * @param initSol Optional initial queue length matrix [1 x (M*K)] in row-major order.
 *        If provided, jobs will be distributed according to this matrix at simulation start.
 */
internal class SSJSimulator(
    private val sn: NetworkStruct,
    private val options: SolverOptions,
    private val initSol: Matrix? = null,
    private val stream: Collector? = null
) {
    // Network parameters derived from sn
    private val numClasses = sn.nclasses
    private val numStations = sn.nstations
    private val numNodes = sn.nnodes
    private val seed = options.seed.toLong()

    /**
     * Server state for setup and delayoff support.
     * Models serverless/FaaS cold start and teardown behavior.
     */
    private enum class ServerState {
        /** Server is completely off, requires setup before servicing jobs */
        OFF,
        /** Server is warming up (cold start), cannot serve jobs yet */
        SETUP,
        /** Server is active and can serve jobs (may be idle or busy) */
        ACTIVE,
        /** Server is in teardown phase before turning off */
        DELAYOFF
    }

    // Node classification
    private val sourceNodes = mutableListOf<Int>()
    private val sourceStations = mutableListOf<Int>()
    private val serviceNodes = mutableListOf<Int>()
    private val serviceStations = mutableListOf<Int>()
    private val isDelayNode = mutableListOf<Boolean>()
    private val sinkNodes = mutableListOf<Int>()
    private val loggerNodes = mutableListOf<Int>()
    private val routerNodes = mutableListOf<Int>()
    private val classSwitchNodes = mutableListOf<Int>()
    private val forkNodes = mutableListOf<Int>()
    private val joinNodes = mutableListOf<Int>()
    private val joinStations = mutableListOf<Int>()
    private val placeNodes = mutableListOf<Int>()
    private val transitionNodes = mutableListOf<Int>()
    private val cacheNodes = mutableListOf<Int>()

    // Cache state information
    private data class CacheStateInfo(
        val nodeIdx: Int,
        val numItems: Int,
        val levelCapacities: IntArray,
        val replacementStrategy: ReplacementStrategy,
        val levels: Array<LinkedList<Int>>,  // Each level is a list of item indices
        val accessSampler: RandomVariateGen?,  // Sampler for item access patterns
        val hitClass: IntArray,   // [classIdx] -> class to switch to on hit
        val missClass: IntArray,  // [classIdx] -> class to switch to on miss
        var totalHits: LongArray,   // [classIdx] -> total hits
        var totalMisses: LongArray  // [classIdx] -> total misses
    )
    private lateinit var cacheStates: Array<CacheStateInfo?>

    // Place token storage: [placeNodeIdx] -> token counts per class
    private lateinit var placeTokens: Array<IntArray>

    // Transition mode information
    private data class TransitionModeInfo(
        val modeIdx: Int,
        val modeName: String,
        val timingStrategy: TimingStrategy,
        val priority: Int,
        val weight: Double,
        val numServers: Int,
        val enablingConditions: Array<IntArray>,  // [inputPlaceIdx][classIdx] -> required tokens
        val inhibitingConditions: Array<IntArray>, // [inputPlaceIdx][classIdx] -> max tokens (or Int.MAX_VALUE)
        val firingOutcomes: Array<IntArray>        // [outputNodeIdx][classIdx] -> tokens produced
    )

    // Transition parameters: [transitionNodeIdx] -> list of modes
    private lateinit var transitionModes: Array<List<TransitionModeInfo>>

    // Transition firing distributions: [transitionNodeIdx][modeIdx] -> distribution generator
    private lateinit var transitionFiringGens: Array<Array<RandomVariateGen?>>

    // Tokens currently being processed by transitions (for multi-server)
    // [transitionListIdx][modeIdx] -> count of in-service server instances
    private lateinit var transitionInService: Array<IntArray>

    // Tokens in transit at transitions: [transListIdx][classIdx] -> token count
    private lateinit var tokensInTransit: Array<IntArray>
    // Time-weighted tokens in transit: [transListIdx][classIdx] -> accumulated time * tokens
    private lateinit var totalTransitTokenTime: Array<DoubleArray>
    private lateinit var lastTransitUpdateTime: DoubleArray

    // Place statistics tracking: [placeListIdx][classIdx]
    private lateinit var totalPlaceTokenTime: Array<DoubleArray>   // Time-weighted token count
    private lateinit var placeCompletions: Array<IntArray>          // Tokens departing place
    private lateinit var lastPlaceUpdateTime: DoubleArray           // Last update time per place
    // Time-weighted tokens in transit that originated from each place: [placeListIdx][classIdx]
    private lateinit var placeTransitTokenTime: Array<DoubleArray>
    // Current tokens in transit that originated from each place: [placeListIdx][classIdx]
    private lateinit var placeTokensInTransit: Array<IntArray>
    private lateinit var lastPlaceTransitUpdateTime: DoubleArray

    private data class ForkJobInfo(
        val parentJobId: Long,
        val parentClassId: Int,
        val parentSystemArrivalTime: Double,
        val forkNodeIdx: Int,
        val totalTasks: Int,
        var completedTasks: Int = 0,
        var firstJoinArrivalTime: Double = -1.0,
        val forkedJobJoinArrivalTimes: MutableList<Double> = mutableListOf(),
        val forkedJobJoinClasses: MutableList<Int> = mutableListOf()  // Track class of each forked job at Join
    )

    private data class ForkedJob(
        val forkJobId: Long,
        val parentJobId: Long,
        val classId: Int,
        val priority: Int,
        val systemArrivalTime: Double,
        val queueArrivalTime: Double,
        val randomRank: Double
    )

    // Maps parent job ID to fork info (for join synchronization)
    private val forkJobInfoMap = mutableMapOf<Long, ForkJobInfo>()

    // Maps forked job ID to parent job ID (for tracking)
    private val forkedJobParentMap = mutableMapOf<Long, Long>()

    // Fork parameters: [forkNodeIdx] -> fanOut (tasks per link)
    private lateinit var forkFanOut: IntArray

    // Join parameters: [joinNodeIdx] -> corresponding fork node index
    private lateinit var joinToForkMap: IntArray
    // Fork parameters: [forkListIdx] -> corresponding join node index (in joinNodes list)
    private lateinit var forkToJoinMap: IntArray

    // Join strategies: [joinNodeIdx][classId] -> JoinStrategy
    private lateinit var joinStrategies: Array<Array<JoinStrategy>>

    // Join required counts: [joinNodeIdx][classId] -> required tasks (-1 = all)
    private lateinit var joinRequired: Array<IntArray>

    // Join statistics tracking: [joinListIdx][classIdx]
    private lateinit var totalJoinQueueTime: Array<DoubleArray>    // Time-weighted queue length
    private lateinit var joinCompletions: Array<IntArray>           // Completed synchronizations
    private lateinit var lastJoinUpdateTime: Array<DoubleArray>     // Last update time per join/class
    private lateinit var currentJoinQueueLength: Array<IntArray>    // Current parent jobs in fork-join section (JMT semantics)
    private lateinit var joinResponseTimeTally: Array<Array<Tally>> // Response time at Join
    private lateinit var arrivedAtJoin: Array<IntArray>              // Arrived forked jobs (for arrival rate)

    // Counter for generating unique forked job IDs
    private var nextForkedJobId: Long = 0L

    // Class switch matrix storage: [nodeIdx] -> row-major probability matrix
    private lateinit var classSwitchMatrices: Array<Array<DoubleArray>?>  // [nodeIdx][fromClass][toClass]

    // Routing strategy tracking
    private lateinit var nodeRoutingStrategies: Array<Array<RoutingStrategy?>>  // [nodeIdx][classIdx]
    private lateinit var roundRobinCounters: IntArray  // [nodeIdx] - for RROBIN
    private lateinit var wrrobinWeights: Array<Array<DoubleArray?>>  // [nodeIdx][classIdx] -> weights by destNodeIdx
    private lateinit var rroutlinks: Array<Array<IntArray?>>  // [nodeIdx][classIdx] -> ordered list of dest node indices
    private var kchoicesK: Int = 2  // Default K for power-of-K-choices

    // Derived arrays
    private lateinit var lambdas: Array<DoubleArray>
    private lateinit var mus: Array<DoubleArray>
    private lateinit var numServers: IntArray
    private lateinit var bufferCapacities: IntArray
    private lateinit var classCapacities: Array<IntArray>

    // Finite capacity regions
    private val numRegions = sn.nregions
    private val fcRegionIndices = MutableList(numStations) { -1 }
    private val fcRegionGlobalMax = MutableList(numRegions) { Int.MAX_VALUE }
    private val fcRegionClassMax = Matrix(numRegions, numClasses)
    // Per-class drop rules: fcRegionDropRule[f][r] = true if class r in region f should drop, false if should block
    private val fcRegionDropRule = Array(numRegions) { BooleanArray(numClasses) { true } }

    // Blocking policy support (BAS and BBS)
    // Station-level drop rules: stationDropRule[svcIdx][classId] = DropStrategy ID
    private lateinit var stationDropRule: Array<IntArray>

    // =================================================================================
    // BBS (Blocking Before Service) - Server blocking when destination is full
    // The SERVER at the source is blocked until the destination has capacity.
    // This prevents the source from processing any new jobs while blocked.
    // =================================================================================
    private data class BBSBlockedServer(
        val customer: Customer,           // The job that completed service
        val destQueueIdx: Int,            // Destination queue index
        val destClassId: Int,             // Destination class ID
        val sourceQueueIdx: Int,          // Source queue index
        val serverId: Int,                // Server ID at source that is blocked
        val sourceClassId: Int,           // Class ID at source (for stats)
        val blockStartTime: Double        // When blocking started (for stats)
    )

    // BBS blocked servers: source queue index -> list of blocked servers (FIFO order)
    private val bbsBlockedServers: MutableMap<Int, MutableList<BBSBlockedServer>> = mutableMapOf()

    // BBS: Reverse index - destination queue index -> set of source queues with blocked servers
    private val bbsDestinationToSources: MutableMap<Int, MutableSet<Int>> = mutableMapOf()

    // =================================================================================
    // BAS (Blocking After Service) - Server blocking when destination is full
    // When a job completes service and finds the destination full, it waits at the
    // SOURCE occupying the server space. The SERVER is blocked until the destination
    // has capacity, at which point the job moves and the server resumes processing.
    // This matches JMT behavior and the classical BAS definition.
    // =================================================================================
    private data class BASWaitingJob(
        val customer: Customer,           // The job waiting to enter destination
        val destQueueIdx: Int,            // Destination queue index
        val destClassId: Int,             // Class ID at destination
        val sourceQueueIdx: Int,          // Source queue index
        val serverId: Int,                // Server ID at source that is blocked
        val sourceClassId: Int,           // Class ID at source (for stats)
        val arrivalTime: Double           // When job arrived at waiting buffer (for FIFO ordering)
    )

    // BAS blocked jobs: source queue index -> list of blocked jobs (FIFO order)
    // Jobs physically reside at source, blocking the server, counted in destination's qlen
    private val basOutgoingBuffer: MutableMap<Int, MutableList<BASWaitingJob>> = mutableMapOf()

    // BAS: Reverse index - destination queue index -> set of source queues with blocked servers
    private val basDestinationToSources: MutableMap<Int, MutableSet<Int>> = mutableMapOf()

    // Counts
    private var numServiceNodes = 0
    private var numSources = 0

    // Tracking for max queue length (used internally)
    private var maxQueueLengthReached = 0

    // Class priorities (from sn.classprio)
    private val classPrio: IntArray = IntArray(numClasses) { k ->
        sn.classprio.get(k).toInt()
    }

    // Class deadlines (from sn.classdeadline) - relative deadline from arrival
    private val classDeadline: DoubleArray = DoubleArray(numClasses) { k ->
        if (sn.classdeadline != null) {
            sn.classdeadline.get(k)
        } else {
            Double.POSITIVE_INFINITY
        }
    }

    // Scheduling strategies per service node
    private lateinit var schedStrategies: Array<SchedStrategy>

    // Polling scheduling state (for stations with SchedStrategy.POLLING)
    private lateinit var isPollingStation: BooleanArray           // [svcIdx] -> true if polling
    private lateinit var pollingType: Array<PollingType?>         // [svcIdx] -> polling type (GATED, EXHAUSTIVE, KLIMITED)
    private lateinit var pollingK: IntArray                        // [svcIdx] -> K value for KLIMITED
    private lateinit var pollingCurrentClass: IntArray             // [svcIdx] -> currently served class index
    private lateinit var pollingJobsServedInRound: IntArray        // [svcIdx] -> jobs served in current round
    private lateinit var pollingGateSize: IntArray                 // [svcIdx] -> gate size for GATED (jobs at start of visit)
    private lateinit var pollingInSwitchover: BooleanArray         // [svcIdx] -> true if currently in switchover
    private lateinit var pollingQueues: Array<Array<LinkedList<Customer>>>  // [svcIdx][classIdx] -> per-class queue
    private lateinit var pollingSwitchoverGens: Array<Array<RandomVariateGen?>>  // [svcIdx][classIdx] -> switchover time generators

    // Class type classification (open vs closed)
    private lateinit var isOpenClass: BooleanArray       // [k] -> true if class k is open
    private lateinit var isClosedClass: BooleanArray     // [k] -> true if class k is closed
    private lateinit var closedClassPopulation: IntArray // [k] -> population if closed, 0 otherwise
    private lateinit var referenceStation: IntArray      // [k] -> reference station index for closed class k

    // Signal class detection (for G-networks with negative customers)
    private lateinit var isSignalClass: BooleanArray     // [k] -> true if class k is a Signal
    private lateinit var isNegativeSignal: BooleanArray  // [k] -> true if class k is a NEGATIVE signal
    private var hasNegativeSignals: Boolean = false      // Quick check flag for signal handling
    // Batch removal configuration (for extended G-networks)
    private lateinit var signalRemovalDist: Array<DiscreteDistribution?>  // [k] -> removal distribution for signal class k
    private lateinit var signalRemovalPolicy: Array<RemovalPolicy?>       // [k] -> removal policy for signal class k
    private lateinit var isCatastropheSignal: BooleanArray                // [k] -> true if class k is a Catastrophe
    private var hasCatastropheSignals: Boolean = false                       // Quick check flag for catastrophe handling

    // REPLY signal detection (for synchronous call blocking semantics)
    private lateinit var isReplySignal: BooleanArray     // [k] -> true if class k is a REPLY signal
    private var hasReplySignals: Boolean = false         // Quick check flag for reply signal handling

    // Synchronous call configuration - maps class to expected reply signal class
    private lateinit var synchCallReplyClass: IntArray   // [k] -> expected reply class index, -1 if none

    // Server blocking state for synchronous calls
    private lateinit var serverBlocked: Array<BooleanArray>  // [queueIdx][serverId] -> true if blocked waiting for reply

    // Pending reply tracking for synchronous call semantics
    private data class PendingReply(
        val jobId: Long,
        val originalClassId: Int,
        val queueIdx: Int,
        val serverId: Int,
        val blockStartTime: Double
    )
    private val pendingReplyMap = mutableMapOf<Long, PendingReply>()

    // Job ID counter for reply tracking
    private var nextJobId: Long = 0L

    // FCFS comparator (ignores priority, only uses arrival time)
    private val fcfsComparator = Comparator<Customer> { c1, c2 ->
        c1.queueArrivalTime.compareTo(c2.queueArrivalTime)  // FCFS: earlier arrival first
    }

    // LCFS comparator (ignores priority, only uses arrival time)
    private val lcfsComparator = Comparator<Customer> { c1, c2 ->
        c2.queueArrivalTime.compareTo(c1.queueArrivalTime)  // LCFS: later arrival first
    }

    // SIRO comparator (uses random rank)
    private val siroComparator = Comparator<Customer> { c1, c2 ->
        c1.randomRank.compareTo(c2.randomRank)
    }

    // Priority comparator (uses priority first, then FCFS for tie-breaking)
    // LINE/JMT convention: lower priority value = higher priority (priority 0 is highest)
    private val priorityComparator = Comparator<Customer> { c1, c2 ->
        when {
            c1.priority != c2.priority -> c1.priority.compareTo(c2.priority)  // Lower value = higher priority
            else -> c1.queueArrivalTime.compareTo(c2.queueArrivalTime)        // FCFS: earlier arrival first
        }
    }

    // Priority comparator with LCFS tie-breaking (for LCFSPRIO variants)
    // LINE/JMT convention: lower priority value = higher priority (priority 0 is highest)
    private val priorityLcfsComparator = Comparator<Customer> { c1, c2 ->
        when {
            c1.priority != c2.priority -> c1.priority.compareTo(c2.priority)  // Lower value = higher priority
            else -> c2.queueArrivalTime.compareTo(c1.queueArrivalTime)        // LCFS: later arrival first
        }
    }

    // SJF comparator: Shortest Job First
    private val sjfComparator = Comparator<Customer> { c1, c2 ->
        if (c1.serviceTime != c2.serviceTime) {
            c1.serviceTime.compareTo(c2.serviceTime)
        } else {
            c1.queueArrivalTime.compareTo(c2.queueArrivalTime) // FCFS tie-break
        }
    }

    // LJF comparator: Longest Job First
    private val ljfComparator = Comparator<Customer> { c1, c2 ->
        if (c1.serviceTime != c2.serviceTime) {
            c2.serviceTime.compareTo(c1.serviceTime)
        } else {
            c1.queueArrivalTime.compareTo(c2.queueArrivalTime) // FCFS tie-break
        }
    }

    // EDD/EDF comparator: Earliest Deadline First (earlier deadline = higher priority)
    private val eddComparator = Comparator<Customer> { c1, c2 ->
        when {
            c1.absoluteDeadline != c2.absoluteDeadline ->
                c1.absoluteDeadline.compareTo(c2.absoluteDeadline)  // Earlier deadline first
            else -> c1.queueArrivalTime.compareTo(c2.queueArrivalTime)  // FCFS tie-break
        }
    }

    // EDF uses same comparator as EDD (preemption handled by server type)
    private val edfComparator = eddComparator

    // SRPT comparator: Shortest Remaining Processing Time first
    // Uses actual remaining work, not class-based expected times
    private fun createSRPTComparator(queueIdx: Int): Comparator<Customer> {
        return Comparator { c1, c2 ->
            // Look up remaining work from preemption history
            val key1 = Triple(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime)
            val key2 = Triple(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime)

            val record1 = preemptedJobHistory[key1]
            val record2 = preemptedJobHistory[key2]

            // Use remaining work if available, otherwise use serviceTime (for new arrivals)
            val remaining1 = record1?.remainingWork ?: c1.serviceTime
            val remaining2 = record2?.remainingWork ?: c2.serviceTime

            if (remaining1 != remaining2) {
                remaining1.compareTo(remaining2) // Shortest remaining first
            } else {
                c1.queueArrivalTime.compareTo(c2.queueArrivalTime) // FCFS tie-break (oldest first)
            }
        }
    }

    // SRPTPRIO comparator: Priority first, then SRPT (actual remaining work)
    private fun createSRPTPRIOComparator(queueIdx: Int): Comparator<Customer> {
        return Comparator { c1, c2 ->
            // Look up remaining work from preemption history
            val key1 = Triple(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime)
            val key2 = Triple(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime)

            val record1 = preemptedJobHistory[key1]
            val record2 = preemptedJobHistory[key2]

            // Use remaining work if available, otherwise use serviceTime (for new arrivals)
            val remaining1 = record1?.remainingWork ?: c1.serviceTime
            val remaining2 = record2?.remainingWork ?: c2.serviceTime

            when {
                c1.priority != c2.priority ->
                    c2.priority.compareTo(c1.priority)  // Higher priority first
                remaining1 != remaining2 -> remaining1.compareTo(remaining2)  // Shortest remaining first
                else -> c1.queueArrivalTime.compareTo(c2.queueArrivalTime)  // FCFS tie-break (oldest first)
            }
        }
    }

    // Returns appropriate comparator based on scheduling strategy
    private fun getComparatorForStrategy(strategy: SchedStrategy, queueIdx: Int): Comparator<Customer> {
        return when (strategy) {
            SchedStrategy.SIRO -> siroComparator
            SchedStrategy.LCFS, SchedStrategy.LCFSPR, SchedStrategy.LCFSPI -> lcfsComparator
            SchedStrategy.HOL, SchedStrategy.FCFSPRIO,
            SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO -> priorityComparator
            SchedStrategy.LCFSPRIO, SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO -> priorityLcfsComparator
            SchedStrategy.SEPT -> createSEPTComparator(queueIdx)
            SchedStrategy.LEPT -> createLEPTComparator(queueIdx)
            SchedStrategy.SJF -> sjfComparator
            SchedStrategy.LJF -> ljfComparator
            SchedStrategy.EDD -> eddComparator
            SchedStrategy.EDF -> edfComparator
            SchedStrategy.SRPT -> createSRPTComparator(queueIdx)
            SchedStrategy.SRPTPRIO -> createSRPTPRIOComparator(queueIdx)
            SchedStrategy.PSJF -> createPSJFComparator(queueIdx)
            SchedStrategy.FB -> createFBComparator(queueIdx)
            SchedStrategy.LRPT -> createLRPTComparator(queueIdx)
            SchedStrategy.SETF -> createSETFComparator(queueIdx)
            else -> fcfsComparator  // FCFS, PS, etc. - ignore priorities
        }
    }

    // SETF comparator: Shortest Elapsed Time First (non-preemptive FB)
    // Jobs with least attained service (smallest age) get priority
    // Unlike FB, this is non-preemptive - jobs run to completion once started
    private fun createSETFComparator(queueIdx: Int): Comparator<Customer> {
        return Comparator { c1, c2 ->
            // For SETF, newly arriving jobs have zero attained service
            // so we use FCFS for tie-breaking. If a job was previously served
            // (multiserver case), we use the preemption history.
            val key1 = Triple(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime)
            val key2 = Triple(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime)

            val record1 = preemptedJobHistory[key1]
            val record2 = preemptedJobHistory[key2]

            // Attained service = elapsed time (or 0 if no history)
            val attained1 = record1?.elapsedTime ?: 0.0
            val attained2 = record2?.elapsedTime ?: 0.0

            // Shorter attained service first (least progress)
            if (attained1 != attained2) {
                attained1.compareTo(attained2)
            } else {
                c1.queueArrivalTime.compareTo(c2.queueArrivalTime) // FCFS tie-break
            }
        }
    }

    // PSJF comparator: Preemptive Shortest Job First
    // Priority based on original job size (not remaining size)
    // For preempted jobs, look up original total from preemption history
    private fun createPSJFComparator(queueIdx: Int): Comparator<Customer> {
        return Comparator { c1, c2 ->
            // Look up original total from preemption history
            val key1 = Triple(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime)
            val key2 = Triple(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime)

            val record1 = preemptedJobHistory[key1]
            val record2 = preemptedJobHistory[key2]

            // Use original total if available, otherwise use serviceTime (for new arrivals)
            val originalSize1 = record1?.originalTotal ?: c1.serviceTime
            val originalSize2 = record2?.originalTotal ?: c2.serviceTime

            if (originalSize1 != originalSize2) {
                originalSize1.compareTo(originalSize2) // Smaller original size first
            } else {
                c1.queueArrivalTime.compareTo(c2.queueArrivalTime) // FCFS tie-break
            }
        }
    }

    // LRPT comparator: Longest Remaining Processing Time first
    // Uses actual remaining work: for preempted jobs from history, for new arrivals use serviceTime
    // Tie-break: SPT (Shortest Processing Time based on original service time sample)
    private fun createLRPTComparator(queueIdx: Int): Comparator<Customer> {
        return Comparator { c1, c2 ->
            // Look up remaining work from preemption history
            val key1 = Triple(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime)
            val key2 = Triple(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime)

            val record1 = preemptedJobHistory[key1]
            val record2 = preemptedJobHistory[key2]

            // Use remaining work if available, otherwise use serviceTime (for new arrivals)
            val remaining1 = record1?.remainingWork ?: c1.serviceTime
            val remaining2 = record2?.remainingWork ?: c2.serviceTime

            if (remaining1 != remaining2) {
                remaining2.compareTo(remaining1) // Longest remaining first
            } else {
                // SPT tie-break: shortest original service time first
                val original1 = record1?.originalTotal ?: c1.serviceTime
                val original2 = record2?.originalTotal ?: c2.serviceTime
                original1.compareTo(original2)
            }
        }
    }

    // FB/LAS comparator: Feedback / Least Attained Service
    // Jobs with least attained service (smallest age) get priority
    private fun createFBComparator(queueIdx: Int): Comparator<Customer> {
        return Comparator { c1, c2 ->
            // Look up attained service from preemption history
            val key1 = Triple(queueIdx, c1.systemArrivalTime, c1.queueArrivalTime)
            val key2 = Triple(queueIdx, c2.systemArrivalTime, c2.queueArrivalTime)

            val record1 = preemptedJobHistory[key1]
            val record2 = preemptedJobHistory[key2]

            // Attained service = elapsed time (or 0 if no history)
            val attained1 = record1?.elapsedTime ?: 0.0
            val attained2 = record2?.elapsedTime ?: 0.0

            // Shorter attained service first (least progress)
            if (attained1 != attained2) {
                attained1.compareTo(attained2)
            } else {
                c1.queueArrivalTime.compareTo(c2.queueArrivalTime) // FCFS tie-break
            }
        }
    }

    // SEPT comparator: Shortest Expected Processing Time first
    // EPT = 1.0 / mu. Higher mu = Lower EPT.
    private fun createSEPTComparator(queueIdx: Int): Comparator<Customer> {
        return Comparator { c1, c2 ->
            val mu1 = mus[queueIdx][c1.classId]
            val mu2 = mus[queueIdx][c2.classId]
            // Higher mu comes first
            if (mu1 != mu2) {
                mu2.compareTo(mu1)
            } else {
                c1.queueArrivalTime.compareTo(c2.queueArrivalTime) // FCFS tie-break
            }
        }
    }

    // LEPT comparator: Longest Expected Processing Time first
    // EPT = 1.0 / mu. Lower mu = Higher EPT.
    private fun createLEPTComparator(queueIdx: Int): Comparator<Customer> {
        return Comparator { c1, c2 ->
            val mu1 = mus[queueIdx][c1.classId]
            val mu2 = mus[queueIdx][c2.classId]
            // Lower mu comes first
            if (mu1 != mu2) {
                mu1.compareTo(mu2)
            } else {
                c1.queueArrivalTime.compareTo(c2.queueArrivalTime) // FCFS tie-break
            }
        }
    }

    // Returns true if this is a PS-family scheduling strategy
    private fun isPSScheduling(strategy: SchedStrategy): Boolean {
        return when (strategy) {
            SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
            SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO -> true
            else -> false
        }
    }

    // Returns true if this PS variant uses priorities
    private fun isPSWithPriority(strategy: SchedStrategy): Boolean {
        return when (strategy) {
            SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO -> true
            else -> false
        }
    }

    // Returns true if this PS variant uses weights (DPS/GPS)
    private fun isPSWithWeights(strategy: SchedStrategy): Boolean {
        return when (strategy) {
            SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO -> true
            else -> false
        }
    }

    // Returns true if this is a preemptive LCFS scheduling strategy
    private fun isPreemptiveLCFSScheduling(strategy: SchedStrategy): Boolean {
        return when (strategy) {
            SchedStrategy.LCFSPR, SchedStrategy.LCFSPI,
            SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO -> true
            else -> false
        }
    }

    // Returns true if this scheduling strategy uses LCFS ordering (including variants)
    private fun isLCFSFamily(strategy: SchedStrategy): Boolean {
        return when (strategy) {
            SchedStrategy.LCFS, SchedStrategy.LCFSPR, SchedStrategy.LCFSPI,
            SchedStrategy.LCFSPRIO, SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO -> true
            else -> false
        }
    }

    // Returns true if this is an SRPT (Shortest Remaining Processing Time) scheduling strategy
    private fun isSRPTScheduling(strategy: SchedStrategy): Boolean {
        return when (strategy) {
            SchedStrategy.SRPT, SchedStrategy.SRPTPRIO -> true
            else -> false
        }
    }

    // Returns true if this is a PSJF (Preemptive Shortest Job First) scheduling strategy
    private fun isPSJFScheduling(strategy: SchedStrategy): Boolean {
        return strategy == SchedStrategy.PSJF
    }

    // Returns true if this is an FB/LAS (Feedback/Least Attained Service) scheduling strategy
    private fun isFBScheduling(strategy: SchedStrategy): Boolean {
        return strategy == SchedStrategy.FB
    }

    // Returns true if this is an LRPT (Longest Remaining Processing Time) scheduling strategy
    private fun isLRPTScheduling(strategy: SchedStrategy): Boolean {
        return strategy == SchedStrategy.LRPT
    }

    // Returns true if this is a size-based preemptive scheduling strategy (SRPT, PSJF, FB, LRPT)
    private fun isSizeBasedPreemptiveScheduling(strategy: SchedStrategy): Boolean {
        return isSRPTScheduling(strategy) || isPSJFScheduling(strategy) ||
               isFBScheduling(strategy) || isLRPTScheduling(strategy)
    }

    // Random number generators
    private lateinit var arrivalGens: Array<Array<RandomVariateGen?>>  // [sourceIdx][classIdx]
    private lateinit var serviceGens: Array<Array<RandomVariateGen?>>  // [svcIdx][classIdx]
    private lateinit var setupGens: Array<Array<RandomVariateGen?>>    // [svcIdx][classIdx] - for setup time
    private lateinit var delayoffGens: Array<Array<RandomVariateGen?>> // [svcIdx][classIdx] - for delayoff time
    private lateinit var routingRng: MRG32k3a
    private lateinit var siroRng: Random

    // PH distribution support for service
    private lateinit var serviceProcessType: Array<Array<ProcessType>>  // [svcIdx][classIdx]
    private lateinit var serviceProc: Array<Array<MatrixCell?>>         // [svcIdx][classIdx] - D0, D1 for PH
    private lateinit var serviceRng: Array<Array<Random?>>              // [svcIdx][classIdx] - for PH sampling

    // PH distribution support for arrivals
    private lateinit var arrivalProcessType: Array<Array<ProcessType>>  // [sourceIdx][classIdx]
    private lateinit var arrivalProc: Array<Array<MatrixCell?>>         // [sourceIdx][classIdx] - D0, D1 for PH
    private lateinit var arrivalRng: Array<Array<Random?>>              // [sourceIdx][classIdx] - for PH sampling

    // BMAP support: cache next batch size for batch arrivals
    private lateinit var arrivalBatchSize: Array<IntArray>              // [sourceIdx][classIdx] - batch size for next arrival

    // Service node state: waitQueue[svcIdx] contains customers waiting (empty for Delay nodes)
    private lateinit var waitQueues: Array<PriorityQueue<Customer>>
    // Server state: serverBusy[svcIdx][serverIdx] (not used for Delay nodes)
    private lateinit var serverBusy: Array<BooleanArray>
    // customersInService[svcIdx]
    private lateinit var customersInService: IntArray

    // Setup and delayoff state tracking
    // Server state per service node: serverState[svcIdx][serverIdx]
    private lateinit var serverState: Array<Array<ServerState>>
    // Which service nodes have setup/delayoff enabled: hasSetupDelayoff[svcIdx]
    private lateinit var hasSetupDelayoff: BooleanArray
    // Track which class triggered setup for each server: serverLastClass[svcIdx][serverIdx]
    private lateinit var serverLastClass: Array<IntArray>
    // Pending delayoff events for cancellation: pendingDelayoffEvents[svcIdx][serverIdx]
    private lateinit var pendingDelayoffEvents: Array<Array<Event?>>

    // Track in-service jobs for signal-based removal (G-networks)
    // Map from (queueIdx, serverId) to (Customer, Departure Event) for cancellation
    private data class InServiceJob(val customer: Customer, val departureEvent: Event)
    private val inServiceJobs: MutableMap<Pair<Int, Int>, InServiceJob> = mutableMapOf()

    // Track Delay node departures for signal-based removal
    // Map from unique job ID to (queueIdx, classId, DelayDeparture Event)
    private data class DelayJob(val queueIdx: Int, val customer: Customer, val departureEvent: Event)
    private val delayJobs: MutableMap<Long, DelayJob> = mutableMapOf()
    private var nextDelayJobId: Long = 0L

    // ==================== Heterogeneous Server Support ====================
    // These variables support heterogeneous multiserver queues where different server types
    // can have different service rates and class compatibility constraints.

    /** Number of server types per service node: [svcIdx] -> count of server types (0 if homogeneous) */
    private lateinit var numServerTypes: IntArray

    /** Servers per type: serversPerType[svcIdx][typeId] -> number of servers of that type */
    private lateinit var serversPerType: Array<IntArray>

    /** Server-class compatibility: serverCompat[svcIdx][typeId][classId] -> true if compatible */
    private lateinit var serverCompat: Array<Array<BooleanArray>>

    /** Busy count per type: busyCountPerType[svcIdx][typeId] -> number of busy servers of that type */
    private lateinit var busyCountPerType: Array<IntArray>

    /** Maps server ID to type ID: serverToType[svcIdx][serverId] -> typeId */
    private lateinit var serverToType: Array<IntArray>

    /** Heterogeneous scheduling policy per service node: [svcIdx] -> policy (null if homogeneous) */
    private lateinit var heteroSchedPolicy: Array<HeteroSchedPolicy?>

    /** Heterogeneous service generators: heteroServiceGens[svcIdx][typeId][classId] -> generator */
    private lateinit var heteroServiceGens: Array<Array<Array<RandomVariateGen?>>>

    /** Heterogeneous service rates: heteroMus[svcIdx][typeId][classId] -> rate */
    private lateinit var heteroMus: Array<Array<DoubleArray>>

    /** Heterogeneous process types: heteroServiceProcType[svcIdx][typeId][classId] -> ProcessType */
    private lateinit var heteroServiceProcType: Array<Array<Array<ProcessType>>>

    /** Heterogeneous PH processes: heteroServiceProc[svcIdx][typeId][classId] -> D0/D1 matrices */
    private lateinit var heteroServiceProc: Array<Array<Array<MatrixCell?>>>

    /** Heterogeneous service RNG: heteroServiceRng[svcIdx][typeId][classId] -> Random for PH sampling */
    private lateinit var heteroServiceRng: Array<Array<Array<Random?>>>

    /** For ALIS/FAIRNESS: server type order for round-robin: serverTypeOrder[svcIdx] -> mutable list of type indices */
    private lateinit var serverTypeOrder: Array<MutableList<Int>>

    /** For ALFS: pre-sorted server types by flexibility (least compatible classes first): alfsOrder[svcIdx] -> sorted type indices */
    private lateinit var alfsOrder: Array<IntArray>

    // Statistics per service node and class
    private lateinit var responseTimeTally: Array<Array<Tally>>      // [svcIdx][classIdx]
    private lateinit var responseTimeSamples: Array<Array<MutableList<Double>>>  // [svcIdx][classIdx] - individual samples for CDF
    private lateinit var completedCustomers: Array<IntArray>         // [svcIdx][classIdx]
    private lateinit var totalQueueTime: Array<DoubleArray>          // [svcIdx][classIdx] - time-weighted queue length
    private lateinit var lastQueueUpdateTime: Array<DoubleArray>     // [svcIdx][classIdx]
    private lateinit var currentQueueLength: Array<IntArray>         // [svcIdx][classIdx]

    // Utilization tracking per service node and class
    private lateinit var totalBusyTime: Array<DoubleArray>           // [svcIdx][classIdx] - cumulative busy time (service only)
    private lateinit var lastBusyUpdateTime: Array<DoubleArray>      // [svcIdx][classIdx]
    private lateinit var currentBusyServers: Array<IntArray>         // [svcIdx][classIdx] - servers busy with this class

    // Blocking time tracking for REPLY signals (synchronous call semantics)
    // When a server blocks waiting for REPLY, this time should be counted towards utilization and queue length
    private lateinit var totalBlockingTime: Array<DoubleArray>       // [svcIdx][classIdx] - cumulative blocking time
    private lateinit var currentBlockedServers: Array<IntArray>      // [svcIdx][classIdx] - servers blocked waiting for SYNCH reply (at source)

    // Blocking policy tracking - blocked jobs count at DESTINATION for queue length calculation
    private lateinit var basBlockedAtDest: Array<IntArray>           // [destQueueIdx][classIdx] - jobs BAS-blocked waiting to ENTER this destination

    // DEBUG: BAS event counters
    private var basBlockCount = 0
    private var basUnblockCount = 0
    private var q1QueueTimeAtBlock = 0.0
    private var q1QueueTimeAfterBlock = 0.0
    private lateinit var bbsBlockedAtDest: Array<IntArray>           // [destQueueIdx][classIdx] - jobs BBS-blocked waiting to ENTER this destination
    private lateinit var fcrBlockedAtDest: Array<IntArray>           // [destQueueIdx][classIdx] - jobs FCR-blocked waiting to ENTER this destination

    // Setup time tracking per service node and class (for setup/delayoff feature)
    private lateinit var totalSetupTime: Array<DoubleArray>          // [svcIdx][classIdx] - cumulative setup time
    private lateinit var lastSetupUpdateTime: Array<DoubleArray>     // [svcIdx][classIdx] - last update timestamp
    private lateinit var currentServersInSetup: Array<IntArray>      // [svcIdx][classIdx] - servers currently in setup phase

    // Delayoff time tracking per service node and class (for setup/delayoff feature)
    private lateinit var totalDelayoffTime: Array<DoubleArray>       // [svcIdx][classIdx] - cumulative delayoff time
    private lateinit var lastDelayoffUpdateTime: Array<DoubleArray>  // [svcIdx][classIdx] - last update timestamp
    private lateinit var currentServersInDelayoff: Array<IntArray>   // [svcIdx][classIdx] - servers currently in delayoff phase

    // System-level statistics
    private lateinit var systemResponseTimeTally: Array<Tally>       // [classIdx]
    private lateinit var systemTardinessTally: Array<Tally>          // [classIdx] - system tardiness
    private lateinit var tardinessTally: Array<Array<Tally>>         // [svcIdx][classIdx] - per-queue tardiness
    private lateinit var systemCompletedCustomers: IntArray          // [classIdx]

    // Dropped customers due to finite buffer
    private lateinit var droppedCustomers: Array<IntArray>           // [svcIdx][classIdx]

    // Arrived customers at each station (including dropped) for arrival rate calculation
    private lateinit var arrivedCustomers: Array<IntArray>           // [svcIdx][classIdx]

    // ==================== Impatience Statistics ====================

    // Reneging (timer-based abandonment) statistics
    private lateinit var renegedCustomers: Array<IntArray>           // [svcIdx][classIdx] - customers who abandoned due to patience expiry
    private lateinit var totalRenegingWaitTime: Array<DoubleArray>   // [svcIdx][classIdx] - cumulative wait time before reneging

    // Balking (queue-length based refusal) statistics
    private lateinit var balkedCustomers: Array<IntArray>            // [svcIdx][classIdx] - customers who refused to join

    // Retrial/Orbit statistics
    private lateinit var orbitJobs: Array<MutableList<OrbitJob>>     // [svcIdx] -> jobs currently in orbit
    private lateinit var retriedCustomers: Array<IntArray>           // [svcIdx][classIdx] - successful retries (re-entered queue)
    private lateinit var maxRetriesExceeded: Array<IntArray>         // [svcIdx][classIdx] - dropped after exceeding max retries
    private lateinit var currentOrbitSize: Array<IntArray>           // [svcIdx][classIdx] - current jobs in orbit
    private lateinit var totalOrbitTime: Array<DoubleArray>          // [svcIdx][classIdx] - time-weighted orbit population
    private lateinit var lastOrbitUpdateTime: Array<DoubleArray>     // [svcIdx][classIdx] - last orbit stats update time

    // Impatience configuration per station-class (cached from NetworkStruct)
    private lateinit var hasPatienceConfig: Array<BooleanArray>      // [svcIdx][classIdx] - true if reneging configured
    private lateinit var patienceGens: Array<Array<RandomVariateGen?>> // [svcIdx][classIdx] - patience time generators
    private lateinit var hasBalkingConfig: Array<BooleanArray>       // [svcIdx][classIdx] - true if balking configured
    private lateinit var hasRetrialConfig: Array<BooleanArray>       // [svcIdx][classIdx] - true if retrial configured
    private lateinit var retrialGens: Array<Array<RandomVariateGen?>>  // [svcIdx][classIdx] - retrial delay generators
    private lateinit var retrialMaxAttemptsConfig: Array<IntArray>   // [svcIdx][classIdx] - max retrial attempts (-1 = unlimited)

    // Track impatient customers in queue for reneging event cancellation
    // Key: (queueIdx, systemArrivalTime.toBits(), classId)
    private val waitingImpatientCustomers: MutableMap<Triple<Int, Long, Int>, ImpatientCustomer> = mutableMapOf()

    // PS (Processor Sharing) scheduling state
    private lateinit var psJobsInService: Array<MutableList<PSCustomer>>  // [svcIdx] -> jobs in PS service
    private lateinit var psLastUpdateTime: DoubleArray                     // [svcIdx] -> last state update time for remaining work
    private lateinit var psLastBusyUpdateTime: DoubleArray                 // [svcIdx] -> last busy stats update time

    // Preemptive LCFS scheduling state
    private lateinit var preemptiveJobsInService: Array<MutableList<PreemptiveCustomer>>  // [svcIdx] -> jobs in preemptive service
    private lateinit var isPreemptiveScheduling: BooleanArray                              // [svcIdx] -> true if using preemptive scheduling

    /**
     * Record for storing complete preemption state including:
     * - remainingWork: how much service time is left (residual time)
     * - originalTotal: the original totalServiceRequirement (for accurate statistics)
     * - elapsedTime: how much time was already spent in service before preemption
     * - distType: process type for distribution tracking (optional, for PH distribution support)
     * - phaseParam: shape parameter for Erlang/Coxian (k for Erlang, optional for other PH)
     */
    private data class PreemptionRecord(
        val remainingWork: Double,
        val originalTotal: Double,
        val elapsedTime: Double,
        val distType: ProcessType? = null,
        val phaseParam: Int? = null  // k for Erlang (number of phases)
    )
    private lateinit var preemptedJobHistory: MutableMap<Triple<Int, Double, Double>, PreemptionRecord>  // Map<(queueIdx, systemArrivalTime, queueArrivalTime), PreemptionRecord>

    // DPS/GPS weights per station per class (cached from sn.schedparam)
    private lateinit var schedWeights: Array<DoubleArray>                  // [svcIdx][classIdx]

    // Load-dependent service support
    private lateinit var lldScaling: Array<DoubleArray?>                  // [svcIdx] -> scaling factors array (null if not load-dependent)
    private lateinit var isLoadDependent: BooleanArray                     // [svcIdx] -> true if station has load-dependent service

    // Limited joint-dependent (LJD) service support - scaling indexed by population vector (n1, n2, ..., nK)
    private var ljdScaling: Array<DoubleArray?>? = null                   // [svcIdx] -> linearized scaling table (null if no LJD)
    private var ljdCutoffs: Array<IntArray?>? = null                      // [svcIdx] -> per-class cutoffs [N1, N2, ..., NK]
    private var hasLjd: Boolean = false

    // Limited joint class-dependent (LJCD) service support - per-class scaling tables
    private var ljcdScaling: Array<Array<DoubleArray?>?>? = null          // [svcIdx][classIdx] -> linearized scaling table
    private var ljcdCutoffs: Array<IntArray?>? = null                     // [svcIdx] -> per-class cutoffs [N1, N2, ..., NK]
    private var hasLjcd: Boolean = false

    // Finite capacity region state tracking
    private lateinit var currentJobsInRegion: Array<IntArray>        // [regionIdx][classIdx] - current jobs in region
    private lateinit var droppedByRegion: Array<IntArray>            // [regionIdx][classIdx] - dropped due to region capacity

    // FCR time-weighted metrics tracking (for QNfcr, TNfcr, etc.)
    private lateinit var totalRegionJobTime: Array<DoubleArray>       // [regionIdx][classIdx] - time-weighted sum for avg queue length
    private lateinit var lastRegionUpdateTime: DoubleArray            // [regionIdx] - last time region stats were updated
    private lateinit var regionCompletions: Array<IntArray>           // [regionIdx][classIdx] - departures from region (for throughput)
    private lateinit var regionResponseTimeTally: Array<Array<Tally>> // [regionIdx][classIdx] - response time tallies

    // FCR class weights - for weighted occupation computation (JMT parity)
    private lateinit var fcRegionClassWeights: Array<DoubleArray>  // [regionIdx][classIdx] - weight per job class

    // FCR arrival rate tracking - using inter-arrival times (correct implementation, unlike JMT's buggy version)
    private lateinit var lastRegionArrivalTime: Array<DoubleArray>     // [regionIdx][classIdx] - last arrival time
    private lateinit var regionArrivalCount: Array<IntArray>           // [regionIdx][classIdx] - count of arrivals
    private lateinit var regionInterArrivalTimeSum: Array<DoubleArray> // [regionIdx][classIdx] - sum of inter-arrival times

    // FCR blocking (waiting queue) support - customers blocked at region entry
    // Data class to hold a blocked customer along with its destination queue
    private data class BlockedCustomer(
        val customer: Customer,
        val destQueueIdx: Int
    )
    // Queue of customers blocked due to FCR capacity (FIFO order)
    private lateinit var fcRegionBlockedQueue: Array<java.util.LinkedList<BlockedCustomer>>
    // Count of blocked customers per region per class (for tracking blocked jobs - NOT included in FCR QLen per JMT semantics)
    private lateinit var blockedInRegion: Array<IntArray>  // [regionIdx][classIdx]

    // Trace writer for DEBUG verbose level
    private var traceWriter: PrintWriter? = null
    private var traceEnabled: Boolean = false

    // Logger node support - CSV output in JMT-compatible format
    // Column names (fixed order, matching JMT): LOGGERNAME, TIMESTAMP, JOB_ID, CLASS_ID, INTERARRIVAL_SAMECLASS, INTERARRIVAL_ANYCLASS, SIMUL_START_TIME
    private val loggerWriters = mutableMapOf<Int, java.io.BufferedWriter>()
    private val loggerConfigs = mutableMapOf<Int, LoggerConfig>()
    private val loggerLastJobTimePerClass = mutableMapOf<Int, DoubleArray>()  // [loggerNodeIdx][classId] -> last job time
    private val loggerLastJobTimeAny = mutableMapOf<Int, Double>()            // [loggerNodeIdx] -> last job time (any class)
    private var simulationStartTime: String = ""

    // Logger configuration data class
    private data class LoggerConfig(
        val nodeIdx: Int,
        val fileName: String,
        val filePath: String,
        val logLoggerName: Boolean,
        val logTimestamp: Boolean,
        val logJobID: Boolean,
        val logJobClass: Boolean,
        val logTimeSameClass: Boolean,
        val logTimeAnyClass: Boolean,
        val logStartTime: Boolean,
        val loggerName: String,
        val delimiter: String = ",",
        val decimalSeparator: String = "."
    )

    // Trace sampler for Replayer distribution - cycles through trace data sequentially
    private class TraceSampler(private val data: DoubleArray) {
        private var index = 0

        fun nextSample(): Double {
            if (data.isEmpty()) return 1e-9  // Fallback for empty trace
            val value = data[index]
            index = (index + 1) % data.size  // Cycle back to start at EOF
            return if (value <= 0) 1e-9 else value  // Ensure positive value
        }

        fun reset() {
            index = 0
        }
    }

    // Trace samplers for arrival and service distributions
    private lateinit var arrivalTraceSamplers: Array<Array<TraceSampler?>>  // [sourceIdx][classIdx]
    private lateinit var serviceTraceSamplers: Array<Array<TraceSampler?>>  // [svcIdx][classIdx]

    // Routing result: destination node and (possibly switched) class
    private data class RoutingResult(val destNode: Int, val destClassId: Int)

    // Customer data class
    private data class Customer(
        val classId: Int,
        val priority: Int,              // Class priority (lower value = higher priority, matching LINE/JMT convention)
        val systemArrivalTime: Double,  // When customer entered the system
        val queueArrivalTime: Double,   // When customer arrived at current queue
        val randomRank: Double,         // Random rank for SIRO scheduling
        var serviceTime: Double = -1.0, // Service time (used for SJF/LJF), initialized to -1.0
        val jobId: Long = -1L,          // Unique job ID for reply signal tracking in synchronous calls
        var absoluteDeadline: Double = Double.POSITIVE_INFINITY,  // Absolute deadline (arrival time + class deadline)
        var assignedServerType: Int = -1  // Server type assigned to this job (-1 = homogeneous or not yet assigned)
    )

    // PS Customer data class for Processor Sharing scheduling
    private data class PSCustomer(
        val classId: Int,
        val priority: Int,                      // Class priority (lower value = higher priority, matching LINE/JMT convention)
        val systemArrivalTime: Double,          // When customer entered the system
        val queueArrivalTime: Double,           // When customer arrived at current queue
        val totalServiceRequirement: Double,    // Original service time (for statistics)
        var remainingServiceWork: Double,       // Remaining service work (mutable)
        var scheduledDepartureEvent: Event?,    // Reference for cancellation (mutable)
        val forkedJob: ForkedJob? = null,       // Fork-join tracking (null if not a forked job)
        var assignedServerType: Int = -1        // Server type assigned to this job (-1 = homogeneous)
    )

    // Preemptive LCFS Customer data class for preemptive scheduling
    private data class PreemptiveCustomer(
        val classId: Int,
        val priority: Int,                      // Class priority (lower value = higher priority, matching LINE/JMT convention)
        val systemArrivalTime: Double,          // When customer entered the system
        val queueArrivalTime: Double,           // When customer arrived at current queue
        val randomRank: Double,                 // Random rank for SIRO scheduling
        val totalServiceRequirement: Double,    // Original service time (for statistics)
        var remainingServiceWork: Double,       // Remaining work (for LCFSPR)
        var elapsedServiceTime: Double,         // Time already spent in service
        var serviceStartTime: Double,           // When current service period started
        var scheduledDepartureEvent: Event?,    // For cancellation (mutable)
        val serverId: Int,                      // Which server is handling this job
        var assignedServerType: Int = -1        // Server type assigned to this job (-1 = homogeneous)
    )

    // ==================== Impatience Data Classes ====================

    /**
     * Impatient customer tracking for reneging (timer-based abandonment).
     * When a customer joins the queue, a reneging event is scheduled.
     * If the customer starts service before the event fires, the event is cancelled.
     */
    private data class ImpatientCustomer(
        val customer: Customer,
        val patienceTime: Double,           // Sampled patience time
        val patienceDeadline: Double,       // Absolute time when patience expires (arrival + patience)
        var renegingEvent: Event? = null    // Scheduled reneging event (for cancellation)
    )

    /**
     * Orbit job tracking for retrial (customer retry after rejection).
     * When a customer is rejected due to capacity, they enter an "orbit"
     * and periodically attempt to re-enter the queue.
     */
    private data class OrbitJob(
        val customer: Customer,
        val orbitEntryTime: Double,         // When customer entered the orbit
        val retrialAttempts: Int,           // Number of retrial attempts so far
        val maxAttempts: Int,               // Maximum attempts (-1 = unlimited)
        val sourceQueueIdx: Int,            // Queue that rejected the customer
        var retrialEvent: Event? = null     // Scheduled retrial event
    )

    /** Warm-up period as fraction of total simulation time (fallback if MSER-5 not used) */
    private var effectiveWarmupFraction = DESOptions.DEFAULT_WARMUP_FRAC

    /** Flag indicating whether warmup period is complete */
    private var warmupDone = false

    /** Time when warmup ended (for calculating actual simulation time) */
    private var warmupEndTime = 0.0

    // MSER-5 transient detection
    /** Batch size for MSER (configurable via options.mserbatch) */
    private var effectiveMserBatchSize = DESOptions.DEFAULT_MSER_BATCH

    // Configurable CI parameters (read from DESOptions)
    /** CI method: "obm", "bm", or "none" */
    private var effectiveCiMethod = "obm"

    /** OBM overlap fraction (0.0-1.0) */
    private var effectiveObmOverlap = DESOptions.DEFAULT_OBM_OVERLAP

    /** Minimum batch size for CI */
    private var effectiveCiMinBatch = DESOptions.DEFAULT_CI_MIN_BATCH

    /** Minimum observations for CI */
    private var effectiveCiMinObs = DESOptions.DEFAULT_CI_MIN_OBS

    /** Sampling interval for collecting queue length observations */
    private var mserSamplingInterval = 0.0

    /** Queue length observations for MSER-5: [svcIdx][classIdx] -> list of observations */
    private lateinit var queueLengthObservations: Array<Array<MutableList<Double>>>

    /** Throughput observations for MSER-5: [svcIdx][classIdx] -> list of completed counts at sample times */
    private lateinit var throughputObservations: Array<Array<MutableList<Int>>>

    /** Place completion observations for MSER-5: [placeListIdx][classIdx] -> list of completed counts at sample times */
    private lateinit var placeCompletionObservations: Array<Array<MutableList<Int>>>

    /** Time stamps for observations */
    private lateinit var observationTimes: MutableList<Double>

    /** MSER-5 truncation point (batch index) determined at end of simulation */
    private var mserTruncationBatch = 0

    /** Flag indicating if MSER-5 is enabled */
    private var mserEnabled = true

    // Event-count based stopping
    /** Total number of service completions (event count) */
    private var totalEventCount: Long = 0

    /** Maximum events to simulate (from options.samples) */
    private var maxEvents: Long = 0

    /** Event count threshold for warmup completion */
    private var warmupEventThreshold: Long = 0

    /** Event count at last MSER sample */
    private var lastMserEventCount: Long = 0

    /** Event interval between MSER samples */
    private var mserEventInterval: Long = 0

    /** Event count at last streaming push */
    private var lastStreamEventCount: Long = 0

    /** Last simulation time when streaming was pushed */
    private var lastStreamTime: Double = 0.0

    // ==================== Variance Reduction Support ====================

    /** Flag indicating if antithetic variates are enabled */
    private var useAntitheticVariates = false

    /** Flag indicating if control variates are enabled */
    private var useControlVariates = false

    /** Antithetic arrival generators: [sourceIdx][classIdx] */
    private lateinit var antitheticArrivalGens: Array<Array<RandomVariateGen?>>

    /** Antithetic service generators: [svcIdx][classIdx] */
    private lateinit var antitheticServiceGens: Array<Array<RandomVariateGen?>>

    /** Antithetic RNG for PH/MAP arrivals: [sourceIdx][classIdx] */
    private lateinit var antitheticArrivalRng: Array<Array<Random?>>

    /** Antithetic RNG for PH/MAP services: [svcIdx][classIdx] */
    private lateinit var antitheticServiceRng: Array<Array<Random?>>

    // Control variates tracking
    /** Sum of sampled arrival times: [sourceIdx][classIdx] */
    private lateinit var arrivalSampleSum: Array<DoubleArray>

    /** Count of arrival samples: [sourceIdx][classIdx] */
    private lateinit var arrivalSampleCount: Array<LongArray>

    /** Expected mean arrival time (1/lambda): [sourceIdx][classIdx] */
    private lateinit var arrivalExpectedMean: Array<DoubleArray>

    /** Sum of sampled service times: [svcIdx][classIdx] */
    private lateinit var serviceSampleSum: Array<DoubleArray>

    /** Count of service samples: [svcIdx][classIdx] */
    private lateinit var serviceSampleCount: Array<LongArray>

    /** Expected mean service time (1/mu): [svcIdx][classIdx] */
    private lateinit var serviceExpectedMean: Array<DoubleArray>

    // Convergence detection state
    /** Flag indicating if convergence checking is enabled */
    private var convergenceEnabled = true

    /** Convergence tolerance - stop when (CI half-width / mean) < tolerance */
    private var convergenceTolerance = 0.05

    /** Minimum number of batches before checking convergence */
    private var convergenceMinBatches = 20

    /** Number of events between convergence checks */
    private var convergenceCheckInterval: Long = 0

    /** Event count at last convergence check */
    private var lastConvergenceCheckEventCount: Long = 0

    /** Batch statistics for queue lengths: [svcIdx][classIdx] -> list of batch means */
    private lateinit var queueBatchMeans: Array<Array<MutableList<Double>>>

    /** Batch statistics for utilizations: [svcIdx][classIdx] -> list of batch means */
    private lateinit var utilBatchMeans: Array<Array<MutableList<Double>>>

    /** Batch statistics for response times: [svcIdx][classIdx] -> list of batch means */
    private lateinit var respTimeBatchMeans: Array<Array<MutableList<Double>>>

    /** Batch statistics for throughputs: [svcIdx][classIdx] -> list of batch means */
    private lateinit var throughputBatchMeans: Array<Array<MutableList<Double>>>

    /** Current batch observation count within a batch */
    private var currentBatchObservations = 0

    /** Cumulative queue time at start of current batch */
    private lateinit var batchStartQueueTime: Array<DoubleArray>

    /** Cumulative busy time at start of current batch */
    private lateinit var batchStartBusyTime: Array<DoubleArray>

    /** Completed customers at start of current batch */
    private lateinit var batchStartCompletions: Array<IntArray>

    /** Response time sum within current batch */
    private lateinit var currentBatchRespTimeSum: Array<DoubleArray>

    /** Response time count within current batch */
    private lateinit var currentBatchRespTimeCount: Array<IntArray>

    /** Time at start of current batch */
    private var batchStartTime = 0.0

    /** Flag indicating if simulation has converged */
    private var hasConverged = false

    /** Stopping reason: "convergence", "max_events", or "max_time" */
    private var stoppingReason = "max_events"

    /** Final CI half-widths for queue lengths */
    private lateinit var finalQNCI: Array<DoubleArray>

    /** Final CI half-widths for utilizations */
    private lateinit var finalUNCI: Array<DoubleArray>

    /** Final CI half-widths for response times */
    private lateinit var finalRNCI: Array<DoubleArray>

    /** Final CI half-widths for throughputs */
    private lateinit var finalTNCI: Array<DoubleArray>

    /** Final relative precision for queue lengths */
    private lateinit var finalQNRelPrec: Array<DoubleArray>

    /** Final relative precision for utilizations */
    private lateinit var finalUNRelPrec: Array<DoubleArray>

    /** Final relative precision for response times */
    private lateinit var finalRNRelPrec: Array<DoubleArray>

    /** Final relative precision for throughputs */
    private lateinit var finalTNRelPrec: Array<DoubleArray>

    init {
        // Identify nodes
        for (i in 0 until numNodes) {
            when (sn.nodetype[i]) {
                NodeType.Source -> {
                    sourceNodes.add(i)
                    sourceStations.add(sn.nodeToStation[i].toInt())
                }
                NodeType.Queue -> {
                    serviceNodes.add(i)
                    serviceStations.add(sn.nodeToStation[i].toInt())
                    isDelayNode.add(false)
                }
                NodeType.Delay -> {
                    serviceNodes.add(i)
                    serviceStations.add(sn.nodeToStation[i].toInt())
                    isDelayNode.add(true)
                }
                NodeType.Sink -> {
                    sinkNodes.add(i)
                }
                NodeType.Logger -> {
                    loggerNodes.add(i)
                }
                NodeType.Router -> {
                    routerNodes.add(i)
                }
                NodeType.ClassSwitch -> {
                    classSwitchNodes.add(i)
                }
                NodeType.Fork -> {
                    forkNodes.add(i)
                }
                NodeType.Join -> {
                    joinNodes.add(i)
                    joinStations.add(sn.nodeToStation[i].toInt())
                }
                NodeType.Place -> {
                    placeNodes.add(i)
                }
                NodeType.Transition -> {
                    transitionNodes.add(i)
                }
                NodeType.Cache -> {
                    cacheNodes.add(i)
                }
                else -> { /* ignore other node types */ }
            }
        }

        // Classify classes as open or closed based on sn.njobs
        // Open class: njobs = Inf, Closed class: njobs = finite population
        isOpenClass = BooleanArray(numClasses) { k ->
            java.lang.Double.isInfinite(sn.njobs.get(k))
        }
        isClosedClass = BooleanArray(numClasses) { k ->
            !java.lang.Double.isInfinite(sn.njobs.get(k))
        }
        closedClassPopulation = IntArray(numClasses) { k ->
            if (isClosedClass[k]) sn.njobs.get(k).toInt() else 0
        }
        referenceStation = IntArray(numClasses) { k ->
            sn.refstat.get(k).toInt()
        }

        // Signal class detection from sn.issignal and sn.signaltype
        isSignalClass = BooleanArray(numClasses) { k ->
            sn.issignal != null && sn.issignal.get(k, 0) > 0
        }
        isNegativeSignal = BooleanArray(numClasses) { k ->
            isSignalClass[k] && sn.signaltype != null && sn.signaltype[k] == SignalType.NEGATIVE
        }
        hasNegativeSignals = isNegativeSignal.any { it }

        // Batch removal configuration from sn.signalRemovalDist, sn.signalRemovalPolicy, sn.isCatastrophe
        signalRemovalDist = Array(numClasses) { k ->
            sn.signalRemovalDist?.getOrNull(k)
        }
        signalRemovalPolicy = Array(numClasses) { k ->
            sn.signalRemovalPolicy?.getOrNull(k)
        }
        isCatastropheSignal = BooleanArray(numClasses) { k ->
            sn.isCatastrophe != null && sn.isCatastrophe.get(k, 0) > 0
        }
        hasCatastropheSignals = isCatastropheSignal.any { it }

        // REPLY signal detection from sn.signaltype
        isReplySignal = BooleanArray(numClasses) { k ->
            isSignalClass[k] && sn.signaltype != null && sn.signaltype[k] == SignalType.REPLY
        }
        hasReplySignals = isReplySignal.any { it }

        // Synchronous call reply class mapping from sn.syncreply
        synchCallReplyClass = if (sn.syncreply != null) {
            IntArray(numClasses) { k -> sn.syncreply.get(k, 0).toInt() }
        } else {
            IntArray(numClasses) { -1 }
        }

        val hasOpenClasses = isOpenClass.any { it }
        val hasClosedClasses = isClosedClass.any { it }
        val hasPetriNet = placeNodes.isNotEmpty() || transitionNodes.isNotEmpty()
        val hasCache = cacheNodes.isNotEmpty()

        // Validate network structure based on class types
        if (serviceNodes.isEmpty() && !hasPetriNet && !hasCache) {
            throw RuntimeException("solver_ssj: Network must have Queue, Delay, Place, Transition, or Cache nodes")
        }
        if (hasOpenClasses && (sourceNodes.isEmpty() || sinkNodes.isEmpty())) {
            throw RuntimeException("solver_ssj: Network with open classes must have Source and Sink nodes")
        }

        numServiceNodes = serviceNodes.size
        numSources = sourceNodes.size

        // Extract arrival rates
        lambdas = Array(numSources) { srcIdx ->
            DoubleArray(numClasses) { k ->
                val rate = sn.rates[sourceStations[srcIdx], k]
                if (rate <= 0 || rate.isNaN() || rate.isInfinite()) 0.0 else rate
            }
        }

        // Extract service rates (with fallback to heterogeneous rates)
        mus = Array(numServiceNodes) { svcIdx ->
            val stationIdx = serviceStations[svcIdx]
            val station = sn.stations[stationIdx]
            DoubleArray(numClasses) { k ->
                var rate = sn.rates[stationIdx, k]
                // If standard rate is not available, try heterogeneous rates
                if ((rate <= 0 || rate.isNaN() || rate.isInfinite()) && sn.heterorates != null) {
                    val heteroRates = sn.heterorates?.get(station)
                    if (heteroRates != null) {
                        // Use the first available heterogeneous rate for this class
                        for ((typeId, classRates) in heteroRates) {
                            val heteroRate = classRates[k]
                            if (heteroRate != null && heteroRate > 0 && heteroRate < Double.MAX_VALUE) {
                                rate = heteroRate
                                break
                            }
                        }
                    }
                }
                if (rate <= 0 || rate.isNaN() || rate.isInfinite()) Double.MAX_VALUE else rate
            }
        }

        // Extract number of servers (with support for heterogeneous server types)
        numServers = IntArray(numServiceNodes) { svcIdx ->
            if (isDelayNode[svcIdx]) {
                Int.MAX_VALUE
            } else {
                val stationIdx = serviceStations[svcIdx]
                var nservers = sn.nservers[stationIdx].toInt()
                // If nservers is 0 or 1 (default), check if there are heterogeneous server types
                if (nservers <= 1 && sn.serverspertype != null) {
                    val station = sn.stations[stationIdx]
                    val serversMatrix = sn.serverspertype?.get(station)
                    if (serversMatrix != null && serversMatrix.length() > 0) {
                        // Sum across all server types
                        var total = 0
                        for (t in 0 until serversMatrix.length()) {
                            total += serversMatrix.get(t).toInt()
                        }
                        if (total > 0) {
                            nservers = total
                        }
                    }
                }
                nservers
            }
        }

        // Extract buffer capacities
        bufferCapacities = IntArray(numServiceNodes) { svcIdx ->
            val capValue = sn.cap[serviceStations[svcIdx]].toInt()
            if (capValue <= 0 || capValue == Int.MAX_VALUE) {
                Int.MAX_VALUE
            } else {
                capValue
            }
        }

        // Extract class capacities
        classCapacities = Array(numServiceNodes) { svcIdx ->
            IntArray(numClasses) { k ->
                val capValue = sn.classcap.get(serviceStations[svcIdx], k)
                when {
                    capValue.isNaN() || capValue <= 0 || capValue >= Int.MAX_VALUE -> Int.MAX_VALUE
                    else -> capValue.toInt()
                }
            }
        }

        // Extract scheduling strategies for each service node
        schedStrategies = Array(numServiceNodes) { svcIdx ->
            val station = sn.stations[serviceStations[svcIdx]]
            sn.sched[station] ?: SchedStrategy.FCFS
        }

        // Extract station-level drop rules (for BAS blocking)
        // Java: sn.droprule is Map<Station, Map<JobClass, DropStrategy>>
        stationDropRule = Array(numServiceNodes) { svcIdx ->
            IntArray(numClasses) { k ->
                val stationIdx = serviceStations[svcIdx]  // Global station index in sn.stations
                if (sn.droprule != null) {
                    val station = sn.stations[stationIdx]
                    val jobClass = sn.jobclasses[k]
                    val dropStrategy = sn.droprule.get(station)?.get(jobClass)
                    val dropStrategyID = dropStrategy?.getID() ?: DropStrategy.Drop.getID()
                    dropStrategyID
                } else {
                    DropStrategy.Drop.getID()  // Default: drop
                }
            }
        }

        // Initialize FC Regions
        fcRegionClassMax.fill(Int.MAX_VALUE.toDouble())
        for (f in 0 until numRegions) {
            val regionMatrix = sn.region?.get(f) ?: continue

            // Read per-class drop rules
            if (sn.regionrule != null && f < sn.regionrule.getNumRows()) {
                for (r in 0 until numClasses) {
                    if (r < sn.regionrule.getNumCols()) {
                        val dropRuleId = sn.regionrule.get(f, r).toInt()
                        fcRegionDropRule[f][r] = (dropRuleId == DropStrategy.Drop.getID())
                    }
                }
            }

            for (i in 0 until numStations) {
                val globalMax = regionMatrix.get(i, numClasses).toInt()
                if (globalMax >= 0) {
                    if (fcRegionIndices[i] < 0) {
                        fcRegionIndices[i] = f
                    }
                    if (globalMax < fcRegionGlobalMax[f]) {
                        fcRegionGlobalMax[f] = globalMax
                    }
                    for (r in 0 until numClasses) {
                        val classMax = regionMatrix.get(i, r).toInt()
                        if (classMax >= 0 && classMax < fcRegionClassMax.get(f, r)) {
                            fcRegionClassMax.set(f, r, classMax.toDouble())
                        }
                    }
                }
            }
        }

        // Initialize Fork/Join parameters
        initializeForkJoinParams()
        initializePlaceTransitionParams()
        initializeCacheParams()
    }

    /**
     * Initialize Place and Transition node parameters from network structure.
     */
    private fun initializePlaceTransitionParams() {
        val numPlaces = placeNodes.size
        val numTransitions = transitionNodes.size

        if (numPlaces == 0 && numTransitions == 0) {
            // No Petri net nodes, initialize empty arrays
            placeTokens = Array(0) { IntArray(0) }
            transitionModes = Array(0) { emptyList() }
            transitionFiringGens = Array(0) { arrayOf() }
            transitionInService = Array(0) { IntArray(0) }
            totalPlaceTokenTime = Array(0) { DoubleArray(0) }
            placeCompletions = Array(0) { IntArray(0) }
            lastPlaceUpdateTime = DoubleArray(0)
            placeTransitTokenTime = Array(0) { DoubleArray(0) }
            placeTokensInTransit = Array(0) { IntArray(0) }
            lastPlaceTransitUpdateTime = DoubleArray(0)
            return
        }

        // Initialize place tokens from initial state
        // Debug: print state map info
        if (options.verbose == VerboseLevel.DEBUG) {
            println("DES: Initializing place tokens, numPlaces=$numPlaces, numClasses=$numClasses")
            println("DES: sn.state is null: ${sn.state == null}")
            if (sn.state != null) {
                println("DES: sn.state size: ${sn.state.size}")
                for ((key, value) in sn.state) {
                    println("DES:   Key: ${key.name} (${System.identityHashCode(key)}), Value: $value, length=${value.length()}")
                }
            }
            println("DES: placeNodes: $placeNodes")
            for (i in placeNodes.indices) {
                val nodeIdx = placeNodes[i]
                val node = sn.nodes[nodeIdx]
                println("DES:   placeNodes[$i] = nodeIdx=$nodeIdx, node=${node.name} (${System.identityHashCode(node)})")
            }
        }

        placeTokens = Array(numPlaces) { placeListIdx ->
            val placeNodeIdx = placeNodes[placeListIdx]
            val node = sn.nodes[placeNodeIdx]
            IntArray(numClasses) { classIdx ->
                // Get initial state from sn.state if available
                val stateMap = sn.state
                if (stateMap != null && node is jline.lang.nodes.StatefulNode) {
                    val state = stateMap[node]
                    if (options.verbose == VerboseLevel.DEBUG) {
                        println("DES: Looking up state for node ${node.name} (${System.identityHashCode(node)}): found=${state != null}")
                    }
                    if (state != null && state.length() > classIdx) {
                        val tokens = state.get(classIdx).toInt()
                        if (options.verbose == VerboseLevel.DEBUG) {
                            println("DES:   class $classIdx: $tokens tokens")
                        }
                        tokens
                    } else {
                        if (options.verbose == VerboseLevel.DEBUG) {
                            println("DES:   class $classIdx: state null or too short, returning 0")
                        }
                        0
                    }
                } else {
                    if (options.verbose == VerboseLevel.DEBUG) {
                        println("DES:   stateMap null or node not StatefulNode, returning 0")
                    }
                    0
                }
            }
        }

        // Initialize Place statistics tracking
        totalPlaceTokenTime = Array(numPlaces) { DoubleArray(numClasses) { 0.0 } }
        placeCompletions = Array(numPlaces) { IntArray(numClasses) { 0 } }
        lastPlaceUpdateTime = DoubleArray(numPlaces) { 0.0 }
        // Initialize place transit tracking (tokens in transit that originated from each place)
        placeTransitTokenTime = Array(numPlaces) { DoubleArray(numClasses) { 0.0 } }
        placeTokensInTransit = Array(numPlaces) { IntArray(numClasses) { 0 } }
        lastPlaceTransitUpdateTime = DoubleArray(numPlaces) { 0.0 }

        // Initialize transition parameters
        transitionModes = Array(numTransitions) { transListIdx ->
            val transNodeIdx = transitionNodes[transListIdx]
            val node = sn.nodes[transNodeIdx]
            val param = sn.nodeparam?.get(node)

            if (param is TransitionNodeParam) {
                val nmodes = param.nmodes
                (0 until nmodes).map { modeIdx ->
                    val modeName = param.modenames?.getOrNull(modeIdx) ?: "Mode$modeIdx"
                    val timing = param.timing?.getOrNull(modeIdx) ?: TimingStrategy.TIMED
                    val priority = param.firingprio?.get(modeIdx)?.toInt() ?: 0
                    val weight = param.fireweight?.get(modeIdx) ?: 1.0
                    val servers = param.nmodeservers?.get(modeIdx)?.toInt() ?: Int.MAX_VALUE

                    // Extract enabling conditions for this mode
                    val enabling = Array(numPlaces) { pIdx ->
                        val placeNodeIdx = placeNodes[pIdx]
                        IntArray(numClasses) { cIdx ->
                            val enabMatrix = param.enabling?.getOrNull(modeIdx)
                            if (enabMatrix != null && placeNodeIdx < enabMatrix.getNumRows() && cIdx < enabMatrix.getNumCols()) {
                                enabMatrix.get(placeNodeIdx, cIdx).toInt()
                            } else {
                                0
                            }
                        }
                    }

                    // Extract inhibiting conditions for this mode
                    // Inhibiting means "block if tokens >= threshold", so:
                    // - 0 or NaN = no inhibition (use MAX_VALUE)
                    // - Infinity = no inhibition (use MAX_VALUE)
                    // - Positive value = inhibit when tokens >= value
                    val inhibiting = Array(numPlaces) { pIdx ->
                        val placeNodeIdx = placeNodes[pIdx]
                        IntArray(numClasses) { cIdx ->
                            val inhibMatrix = param.inhibiting?.getOrNull(modeIdx)
                            if (inhibMatrix != null && placeNodeIdx < inhibMatrix.getNumRows() && cIdx < inhibMatrix.getNumCols()) {
                                val value = inhibMatrix.get(placeNodeIdx, cIdx)
                                // Treat 0, NaN, or Infinity as no inhibition
                                if (java.lang.Double.isInfinite(value) || java.lang.Double.isNaN(value) || value <= 0) {
                                    Int.MAX_VALUE
                                } else {
                                    value.toInt()
                                }
                            } else {
                                Int.MAX_VALUE  // No inhibition
                            }
                        }
                    }

                    // Extract firing outcomes for this mode
                    val firing = Array(numNodes) { nIdx ->
                        IntArray(numClasses) { cIdx ->
                            val fireMatrix = param.firing?.getOrNull(modeIdx)
                            if (fireMatrix != null && nIdx < fireMatrix.getNumRows() && cIdx < fireMatrix.getNumCols()) {
                                fireMatrix.get(nIdx, cIdx).toInt()
                            } else {
                                0
                            }
                        }
                    }

                    TransitionModeInfo(
                        modeIdx = modeIdx,
                        modeName = modeName,
                        timingStrategy = timing,
                        priority = priority,
                        weight = weight,
                        numServers = servers,
                        enablingConditions = enabling,
                        inhibitingConditions = inhibiting,
                        firingOutcomes = firing
                    )
                }
            } else {
                emptyList()
            }
        }

        // Initialize firing distribution generators
        transitionFiringGens = Array(numTransitions) { transListIdx ->
            val modes = transitionModes[transListIdx]
            Array(modes.size) { modeIdx ->
                val transNodeIdx = transitionNodes[transListIdx]
                val node = sn.nodes[transNodeIdx]
                val param = sn.nodeparam?.get(node) as? TransitionNodeParam

                if (param != null && modes[modeIdx].timingStrategy == TimingStrategy.TIMED) {
                    // Get firing process type and parameters
                    val procType = param.firingprocid?.values?.toList()?.getOrNull(modeIdx)
                    val procCell = param.firingproc?.values?.toList()?.getOrNull(modeIdx)

                    if (procType != null && procCell != null) {
                        createDistributionGenerator(procType, procCell, transListIdx, modeIdx)
                    } else {
                        // Default: exponential with rate 1
                        ExponentialGen(MRG32k3a(), 1.0)
                    }
                } else {
                    null  // Immediate transition
                }
            }
        }

        // Initialize in-service counters
        transitionInService = Array(numTransitions) { transListIdx ->
            IntArray(transitionModes[transListIdx].size) { 0 }
        }

        // Initialize tokens in transit tracking
        tokensInTransit = Array(numTransitions) { IntArray(numClasses) { 0 } }
        totalTransitTokenTime = Array(numTransitions) { DoubleArray(numClasses) { 0.0 } }
        lastTransitUpdateTime = DoubleArray(numTransitions) { 0.0 }
    }

    /**
     * Create a distribution generator for transition firing time.
     */
    private fun createDistributionGenerator(
        procType: ProcessType,
        procCell: MatrixCell,
        transListIdx: Int,
        modeIdx: Int
    ): RandomVariateGen? {
        val stream = MRG32k3a()

        return when (procType) {
            ProcessType.EXP -> {
                // D1 contains the rate matrix
                val d1 = procCell.get(1)
                val rate = if (d1 != null && d1.length() > 0) d1.get(0, 0) else 1.0
                ExponentialGen(stream, rate)
            }
            ProcessType.ERLANG -> {
                // For phase-type representation: D0 is k×k transition matrix
                // k = number of phases (rows in D0), phaseRate = -D0[0,0]
                val d0 = procCell.get(0)
                val k = if (d0 != null && d0.getNumRows() > 0) d0.getNumRows() else 1
                val phaseRate = if (d0 != null && d0.length() > 0) FastMath.abs(d0.get(0, 0)) else 1.0
                ErlangGen(stream, k, phaseRate)
            }
            ProcessType.HYPEREXP -> {
                // HyperExp representation:
                // D0 = [[-lambda1, 0], [0, -lambda2]] (negative rates on diagonal)
                // D1 = [[lambda1*p, lambda1*(1-p)], [lambda2*p, lambda2*(1-p)]]
                val d0 = procCell.get(0)
                val d1 = procCell.get(1)
                // Extract rates from D0 diagonal (negate since stored as negative)
                val lambda1 = if (d0 != null && d0.length() > 0) -d0.get(0, 0) else 1.0
                val lambda2 = if (d0 != null && d0.getNumRows() > 1) -d0.get(1, 1) else 1.0
                // Extract p from D1: D1[0,0] = lambda1*p, so p = D1[0,0] / lambda1
                val p = if (d1 != null && d1.length() > 0 && lambda1 > 0) d1.get(0, 0) / lambda1 else 0.5
                HyperExponentialDistGen(stream, doubleArrayOf(p, 1.0 - p), doubleArrayOf(lambda1, lambda2))
            }
            ProcessType.PH, ProcessType.APH, ProcessType.COXIAN, ProcessType.COX2, ProcessType.MAP, ProcessType.MMPP2, ProcessType.ME, ProcessType.RAP -> {
                // Phase-type, ME, RAP, and MAP/MMPP: use map_sample from jline.api.mam
                null  // Will use map_sample directly
            }
            ProcessType.IMMEDIATE -> {
                null  // No delay
            }
            else -> {
                ExponentialGen(stream, 1.0)  // Default fallback
            }
        }
    }

    /**
     * Custom generator for hyper-exponential distribution (mixture of exponentials).
     */
    private class HyperExponentialDistGen(
        stream: umontreal.ssj.rng.RandomStream,
        private val probs: DoubleArray,
        private val rates: DoubleArray
    ) : RandomVariateGen(stream, null) {
        override fun nextDouble(): Double {
            val u = stream.nextDouble()
            var cumProb = 0.0
            for (i in probs.indices) {
                cumProb += probs[i]
                if (u <= cumProb) {
                    return -FastMath.log(stream.nextDouble()) / rates[i]
                }
            }
            return -FastMath.log(stream.nextDouble()) / rates.last()
        }
    }

    /**
     * Custom generator wrapper for PH/MAP distributions using map_sample.
     */
    private inner class MapSampleGen(
        stream: umontreal.ssj.rng.RandomStream,
        private val proc: jline.util.matrix.MatrixCell,
        private val rng: Random
    ) : RandomVariateGen(stream, null) {
        override fun nextDouble(): Double {
            if (proc.size() >= 2) {
                val D0 = proc.get(0)
                val D1 = proc.get(1)
                val samples = map_sample(D0, D1, 1, rng)
                return samples[0]
            }
            return 0.0
        }
    }

    /**
     * Initialize Fork and Join node parameters from network structure.
     */
    private fun initializeForkJoinParams() {
        val numForks = forkNodes.size
        val numJoins = joinNodes.size

        // Initialize Fork parameters
        forkFanOut = IntArray(numForks) { 1 }

        // Extract fanOut from ForkNodeParam for each fork node
        for (forkListIdx in 0 until numForks) {
            val forkNodeIdx = forkNodes[forkListIdx]
            val node = sn.nodes[forkNodeIdx]
            val param = sn.nodeparam?.get(node)
            if (param is jline.lang.nodeparam.ForkNodeParam) {
                val fanOut = param.fanOut
                forkFanOut[forkListIdx] = if (fanOut.isNaN() || fanOut <= 0) 1 else fanOut.toInt()
            }
        }

        // Initialize Join parameters
        joinToForkMap = IntArray(numJoins) { -1 }
        forkToJoinMap = IntArray(numForks) { -1 }
        joinStrategies = Array(numJoins) { Array(numClasses) { JoinStrategy.STD } }
        joinRequired = Array(numJoins) { IntArray(numClasses) { -1 } }  // -1 means all tasks

        // Build Fork-Join pairing from sn.fj matrix
        val fjMatrix = sn.fj
        if (fjMatrix != null) {
            for (joinListIdx in 0 until numJoins) {
                val joinNodeIdx = joinNodes[joinListIdx]
                // Find the corresponding fork node from sn.fj matrix
                for (forkListIdx in 0 until numForks) {
                    val forkNodeIdx = forkNodes[forkListIdx]
                    if (fjMatrix.get(forkNodeIdx, joinNodeIdx) > 0) {
                        joinToForkMap[joinListIdx] = forkListIdx
                        forkToJoinMap[forkListIdx] = joinListIdx  // Reverse mapping for JMT queue semantics
                        break
                    }
                }
            }
        }

        // Extract join strategy and required counts from JoinNodeParam
        for (joinListIdx in 0 until numJoins) {
            val joinNodeIdx = joinNodes[joinListIdx]
            val node = sn.nodes[joinNodeIdx]
            val param = sn.nodeparam?.get(node)
            if (param is jline.lang.nodeparam.JoinNodeParam) {
                // Extract strategy per class
                param.joinStrategy?.forEach { (jobClass, strategy) ->
                    val classIdx = sn.jobclasses.indexOf(jobClass)
                    if (classIdx >= 0 && classIdx < numClasses) {
                        joinStrategies[joinListIdx][classIdx] = strategy
                    }
                }
                // Extract required counts per class (from fanIn, which is joinRequired)
                param.fanIn?.forEach { (jobClass, required) ->
                    val classIdx = sn.jobclasses.indexOf(jobClass)
                    if (classIdx >= 0 && classIdx < numClasses) {
                        joinRequired[joinListIdx][classIdx] = if (required < 0) -1 else required.toInt()
                    }
                }
            }
        }

        // Initialize Join statistics arrays
        totalJoinQueueTime = Array(numJoins) { DoubleArray(numClasses) { 0.0 } }
        joinCompletions = Array(numJoins) { IntArray(numClasses) { 0 } }
        lastJoinUpdateTime = Array(numJoins) { DoubleArray(numClasses) { 0.0 } }
        currentJoinQueueLength = Array(numJoins) { IntArray(numClasses) { 0 } }
        joinResponseTimeTally = Array(numJoins) { joinListIdx ->
            Array(numClasses) { k -> Tally("Join response time J$joinListIdx C$k") }
        }
        arrivedAtJoin = Array(numJoins) { IntArray(numClasses) { 0 } }
    }

    /**
     * Initialize Cache node parameters from network structure.
     * Sets up cache levels, replacement strategies. Access samplers are deferred until runtime.
     */
    private fun initializeCacheParams() {
        val numCaches = cacheNodes.size
        if (numCaches == 0) {
            cacheStates = Array(numNodes) { null }
            return
        }

        cacheStates = Array(numNodes) { null }

        for (cacheListIdx in 0 until numCaches) {
            val cacheNodeIdx = cacheNodes[cacheListIdx]
            val node = sn.nodes[cacheNodeIdx]

            if (node is Cache) {
                val numItems = node.numberOfItems
                val numLevels = node.getnLevels()
                val itemLevelCap = node.itemLevelCap  // Matrix with level capacities

                // Build level capacities from the itemLevelCap matrix
                val levelCapacities = IntArray(numLevels) { level ->
                    if (itemLevelCap != null && level < itemLevelCap.length()) {
                        itemLevelCap.get(level).toInt()
                    } else {
                        1  // Default capacity
                    }
                }

                // Get replacement strategy
                val replacementStrategy = node.replacementStrategy ?: ReplacementStrategy.LRU

                // Initialize cache levels - all items start outside the cache
                val levels = Array(numLevels) { LinkedList<Int>() }

                // Get hit/miss class mappings from CacheClassSwitcher
                val hitClassArr = IntArray(numClasses) { k -> k }
                val missClassArr = IntArray(numClasses) { k -> k }

                val cacheServer = node.cacheServer
                if (cacheServer != null) {
                    val hitClassMatrix = cacheServer.hitClass
                    val missClassMatrix = cacheServer.missClass
                    for (k in 0 until numClasses) {
                        if (hitClassMatrix != null && k < hitClassMatrix.numCols) {
                            val hitIdx = hitClassMatrix.get(0, k).toInt()
                            if (hitIdx >= 0) hitClassArr[k] = hitIdx
                        }
                        if (missClassMatrix != null && k < missClassMatrix.numCols) {
                            val missIdx = missClassMatrix.get(0, k).toInt()
                            if (missIdx >= 0) missClassArr[k] = missIdx
                        }
                    }
                }

                // Access sampler is null initially, will be created lazily on first access
                cacheStates[cacheNodeIdx] = CacheStateInfo(
                    nodeIdx = cacheNodeIdx,
                    numItems = numItems,
                    levelCapacities = levelCapacities,
                    replacementStrategy = replacementStrategy,
                    levels = levels,
                    accessSampler = null,  // Deferred initialization
                    hitClass = hitClassArr,
                    missClass = missClassArr,
                    totalHits = LongArray(numClasses) { 0L },
                    totalMisses = LongArray(numClasses) { 0L }
                )
            }
        }
    }

    /**
     * Initialize cache access samplers. Must be called after routingRng is initialized.
     */
    private fun initializeCacheSamplers() {
        for (cacheNodeIdx in cacheNodes) {
            val cacheState = cacheStates[cacheNodeIdx] ?: continue
            if (cacheState.accessSampler == null) {
                val node = sn.nodes[cacheNodeIdx]
                if (node is Cache) {
                    val sampler = createAccessSampler(node, cacheState.numItems)
                    // Create new state with sampler (since data class fields are val)
                    cacheStates[cacheNodeIdx] = cacheState.copy(accessSampler = sampler)
                }
            }
        }
    }

    /**
     * Create an access sampler for cache items based on the popularity distribution.
     * Uses the Zipf or other popularity distribution set on the cache via setRead().
     */
    private fun createAccessSampler(cache: Cache, numItems: Int): RandomVariateGen? {
        // Try to get access probabilities from the accessProb matrix if it exists
        val accessProbs = cache.accessProb
        if (accessProbs != null && accessProbs.isNotEmpty()) {
            // accessProb is a 2D array [level][class], use first available
            for (levelProbs in accessProbs) {
                if (levelProbs != null) {
                    for (classProb in levelProbs) {
                        if (classProb != null && classProb.length() >= numItems) {
                            val probs = DoubleArray(numItems) { i -> classProb.get(i) }
                            return DiscreteAccessSampler(routingRng, probs)
                        }
                    }
                }
            }
        }

        // Try to use popularity distribution (e.g., Zipf)
        val popularityDist = cache.popularityGet(0, 0)
        if (popularityDist != null) {
            // Build list of item indices (1-indexed for Zipf)
            val itemIndices = (1..numItems).map { it.toDouble() }
            // Get PMF for all items - works with Zipf and other discrete distributions
            when (popularityDist) {
                is jline.lang.processes.Zipf -> {
                    val pmfMatrix = popularityDist.evalPMF(itemIndices)
                    if (pmfMatrix.length() >= numItems) {
                        val probs = DoubleArray(numItems) { i -> pmfMatrix.get(i) }
                        return DiscreteAccessSampler(routingRng, probs)
                    }
                }
                is jline.lang.processes.DiscreteDistribution -> {
                    val pmfMatrix = popularityDist.evalPMF(itemIndices)
                    if (pmfMatrix.length() >= numItems) {
                        val probs = DoubleArray(numItems) { i -> pmfMatrix.get(i) }
                        return DiscreteAccessSampler(routingRng, probs)
                    }
                }
            }
        }

        // Default to uniform access if no distribution specified
        return DiscreteAccessSampler(routingRng, DoubleArray(numItems) { 1.0 / numItems })
    }

    /**
     * Discrete sampler for cache item accesses.
     */
    private class DiscreteAccessSampler(
        private val rng: MRG32k3a,
        private val probs: DoubleArray
    ) : RandomVariateGen(rng, null) {
        private val cdf: DoubleArray

        init {
            // Build CDF
            cdf = DoubleArray(probs.size)
            var sum = 0.0
            for (i in probs.indices) {
                sum += probs[i]
                cdf[i] = sum
            }
            // Normalize
            if (sum > 0) {
                for (i in cdf.indices) {
                    cdf[i] /= sum
                }
            }
        }

        override fun nextDouble(): Double {
            val u = rng.nextDouble()
            for (i in cdf.indices) {
                if (u <= cdf[i]) return i.toDouble()
            }
            return (cdf.size - 1).toDouble()
        }

        fun nextItem(): Int = nextDouble().toInt()
    }

    /**
     * Process a cache access for a job, returning the new class after hit/miss determination.
     * @param cacheNodeIdx The cache node index
     * @param currentClass The current job class
     * @return The new class after cache access (hit class or miss class)
     */
    private fun processCacheAccess(cacheNodeIdx: Int, currentClass: Int): Int {
        val cacheState = cacheStates[cacheNodeIdx] ?: return currentClass

        // Sample which item is being accessed
        val itemIdx = (cacheState.accessSampler as? DiscreteAccessSampler)?.nextItem()
            ?: (routingRng.nextDouble() * cacheState.numItems).toInt().coerceIn(0, cacheState.numItems - 1)

        // Check if item is in cache (search all levels)
        var hitLevel = -1
        for (level in cacheState.levels.indices) {
            if (cacheState.levels[level].contains(itemIdx)) {
                hitLevel = level
                break
            }
        }

        return if (hitLevel >= 0) {
            // Cache hit
            cacheState.totalHits[currentClass]++
            handleCacheHit(cacheState, itemIdx, hitLevel)
            cacheState.hitClass[currentClass]
        } else {
            // Cache miss
            cacheState.totalMisses[currentClass]++
            handleCacheMiss(cacheState, itemIdx)
            cacheState.missClass[currentClass]
        }
    }

    /**
     * Handle a cache hit - move item according to replacement policy.
     */
    private fun handleCacheHit(cacheState: CacheStateInfo, itemIdx: Int, hitLevel: Int) {
        when (cacheState.replacementStrategy) {
            ReplacementStrategy.LRU -> {
                // Move to front of level 0 (MRU position)
                cacheState.levels[hitLevel].remove(itemIdx)
                cacheState.levels[0].addFirst(itemIdx)
                rebalanceLevelsAfterHit(cacheState, 0)
            }
            ReplacementStrategy.FIFO, ReplacementStrategy.SFIFO -> {
                // FIFO: don't move on hit, item stays where it is
            }
            ReplacementStrategy.RR -> {
                // Random replacement: don't move on hit
            }
            else -> {
                // Default: LRU behavior
                cacheState.levels[hitLevel].remove(itemIdx)
                cacheState.levels[0].addFirst(itemIdx)
                rebalanceLevelsAfterHit(cacheState, 0)
            }
        }
    }

    /**
     * Handle a cache miss - insert item according to replacement policy.
     */
    private fun handleCacheMiss(cacheState: CacheStateInfo, itemIdx: Int) {
        when (cacheState.replacementStrategy) {
            ReplacementStrategy.LRU -> {
                // Insert at front of level 0
                cacheState.levels[0].addFirst(itemIdx)
                rebalanceLevelsAfterMiss(cacheState)
            }
            ReplacementStrategy.FIFO -> {
                // Insert at back of level 0 (FIFO order)
                cacheState.levels[0].addLast(itemIdx)
                rebalanceLevelsAfterMissFIFO(cacheState)
            }
            ReplacementStrategy.SFIFO -> {
                // Strict FIFO: insert at front of level 0
                cacheState.levels[0].addFirst(itemIdx)
                rebalanceLevelsAfterMissFIFO(cacheState)
            }
            ReplacementStrategy.RR -> {
                // Random replacement
                cacheState.levels[0].addFirst(itemIdx)
                rebalanceLevelsAfterMissRR(cacheState)
            }
            else -> {
                // Default: LRU behavior
                cacheState.levels[0].addFirst(itemIdx)
                rebalanceLevelsAfterMiss(cacheState)
            }
        }
    }

    /**
     * Rebalance cache levels after a hit (LRU).
     * Cascades overflow from level to level.
     */
    private fun rebalanceLevelsAfterHit(cacheState: CacheStateInfo, startLevel: Int) {
        for (level in startLevel until cacheState.levels.size) {
            while (cacheState.levels[level].size > cacheState.levelCapacities[level]) {
                val evicted = cacheState.levels[level].removeLast()
                if (level + 1 < cacheState.levels.size) {
                    cacheState.levels[level + 1].addFirst(evicted)
                }
                // If last level overflows, item is evicted from cache entirely
            }
        }
    }

    /**
     * Rebalance cache levels after a miss (LRU).
     */
    private fun rebalanceLevelsAfterMiss(cacheState: CacheStateInfo) {
        rebalanceLevelsAfterHit(cacheState, 0)
    }

    /**
     * Rebalance cache levels after a miss (FIFO).
     */
    private fun rebalanceLevelsAfterMissFIFO(cacheState: CacheStateInfo) {
        for (level in 0 until cacheState.levels.size) {
            while (cacheState.levels[level].size > cacheState.levelCapacities[level]) {
                val evicted = cacheState.levels[level].removeFirst()  // FIFO: remove oldest
                if (level + 1 < cacheState.levels.size) {
                    cacheState.levels[level + 1].addLast(evicted)
                }
            }
        }
    }

    /**
     * Rebalance cache levels after a miss (Random Replacement).
     */
    private fun rebalanceLevelsAfterMissRR(cacheState: CacheStateInfo) {
        for (level in 0 until cacheState.levels.size) {
            while (cacheState.levels[level].size > cacheState.levelCapacities[level]) {
                // Random eviction
                val levelSize = cacheState.levels[level].size
                val evictIdx = (routingRng.nextDouble() * levelSize).toInt().coerceIn(0, levelSize - 1)
                val evicted = cacheState.levels[level].removeAt(evictIdx)
                if (level + 1 < cacheState.levels.size) {
                    cacheState.levels[level + 1].addFirst(evicted)
                }
            }
        }
    }

    /**
     * Run the simulation until the specified number of events (service completions).
     * Uses MSER-5 for automatic transient detection, or falls back to fixed warmup if disabled.
     * @param maxEventCount Maximum number of service completion events to simulate
     */
    fun simulate(maxEventCount: Double) {
        Sim.init()

        // Initialize transient detection and CI configuration from DESOptions
        initializeTransientAndCIConfig()

        // Initialize variance reduction mode from options
        val variatesMode = options.config?.variates ?: "none"
        useAntitheticVariates = variatesMode == "antithetic" || variatesMode == "both"
        useControlVariates = variatesMode == "control" || variatesMode == "both"

        if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
            if (useAntitheticVariates) println("DES: Antithetic variates enabled")
            if (useControlVariates) println("DES: Control variates enabled")
        }

        // Initialize event-count based stopping
        maxEvents = maxEventCount.toLong()
        totalEventCount = 0
        warmupEventThreshold = (maxEvents * effectiveWarmupFraction).toLong()

        // Initialize trace logging if DEBUG verbose level
        initializeTracing()

        // Initialize Logger nodes for CSV output
        initializeLoggers()

        // Initialize routing strategies
        initializeRouting()

        // Initialize random generators
        initializeGenerators()

        // Initialize cache samplers (requires routingRng to be initialized)
        initializeCacheSamplers()

        // Initialize queues and statistics
        initializeState()

        // Validate closed class routing (no routes to Sink)
        validateClosedClassRouting()

        // Initialize closed class populations at reference stations
        initClosedClassPopulations()

        if (mserEnabled) {
            // MSER-5 mode: collect observations throughout simulation, determine truncation at the end
            initializeMSEREventBased(maxEvents)
            warmupDone = true  // Statistics are collected from start
            warmupEndTime = 0.0
            // Note: MSER sampling is triggered by event count in checkEventCountStop()
        } else {
            // Legacy fixed warmup mode - warmup triggered by event count
            warmupDone = false
            // Note: Warmup completion is triggered by event count in checkEventCountStop()
        }

        // Initialize convergence checking (runs alongside MSER)
        initializeConvergence(maxEvents)

        // Note: Progress reporting is triggered by event count in checkEventCountStop()

        // Schedule initial arrivals from each source for OPEN classes only
        // Closed classes have fixed population already initialized at reference stations
        for (srcIdx in 0 until numSources) {
            for (k in 0 until numClasses) {
                if (isOpenClass[k]) {
                    val interarrivalTime = generateInterarrivalTime(srcIdx, k)
                    if (interarrivalTime > 0) {
                        ExternalArrival(srcIdx, k).schedule(interarrivalTime)
                    }
                }
            }
        }

        // Run simulation (stops when maxEvents reached via checkEventCountStop())
        Sim.start()
    }

    fun getDESResult(): DESResult {
        val QN = Matrix(numStations, numClasses)
        QN.fill(0.0)
        val UN = Matrix(numStations, numClasses)
        UN.fill(0.0)
        val RN = Matrix(numStations, numClasses)
        RN.fill(0.0)
        val TN = Matrix(numStations, numClasses)
        TN.fill(0.0)
        val AN = Matrix(numStations, numClasses)
        AN.fill(0.0)
        val TardN = Matrix(numStations, numClasses)
        TardN.fill(0.0)
        val SysTardN = Matrix(1, numClasses)
        SysTardN.fill(0.0)
        val CN = Matrix(1, numClasses)
        CN.fill(0.0)
        val XN = Matrix(1, numClasses)
        XN.fill(0.0)

        // Build set of classes that can receive jobs via class-switching
        // A class can receive jobs if it's in the same chain as another class with jobs
        val classSwitchClasses = mutableSetOf<Int>()
        for (c in 0 until sn.nchains) {
            val classesInChain = mutableListOf<Int>()
            for (k in 0 until numClasses) {
                if (sn.chains.get(c, k) > 0) {
                    classesInChain.add(k)
                }
            }
            if (classesInChain.size > 1) {
                // Check if any class in this chain has jobs
                var chainHasJobs = false
                for (k in classesInChain) {
                    if (sn.njobs.get(k) > 0) {
                        chainHasJobs = true
                        break
                    }
                }
                if (chainHasJobs) {
                    // All classes in this chain can receive jobs via class-switching
                    classSwitchClasses.addAll(classesInChain)
                }
            }
        }

        // Finalize time-weighted queue length statistics for all service nodes
        // This ensures the time from the last update to the end of simulation is included
        val finalTime = Sim.time()
        for (qIdx in 0 until numServiceNodes) {
            for (k in 0 until numClasses) {
                val elapsed = finalTime - lastQueueUpdateTime[qIdx][k]
                if (elapsed > 0) {
                    // Include blocked jobs in the final queue length calculation
                    val effectiveQueueLength = currentQueueLength[qIdx][k] +
                                               currentBlockedServers[qIdx][k] +
                                               basBlockedAtDest[qIdx][k] +
                                               bbsBlockedAtDest[qIdx][k] +
                                               fcrBlockedAtDest[qIdx][k]
                    totalQueueTime[qIdx][k] += effectiveQueueLength * elapsed
                    lastQueueUpdateTime[qIdx][k] = finalTime
                }
                // Also finalize busy time
                val busyElapsed = finalTime - lastBusyUpdateTime[qIdx][k]
                if (busyElapsed > 0) {
                    totalBusyTime[qIdx][k] += currentBusyServers[qIdx][k] * busyElapsed
                    lastBusyUpdateTime[qIdx][k] = finalTime
                }
            }
        }

        // DEBUG: Print BAS event counts and queue length details
        // if (basBlockCount > 0 || basUnblockCount > 0) {
        //     val simTime = getActualSimTime()
        //     println("\n=== DES BAS Debug ===")
        //     println("BAS blocks: $basBlockCount, BAS unblocks: $basUnblockCount")
        //     println("Sim time: $simTime")
        //     for (qIdx in 0 until numServiceNodes) {
        //         val stationName = sn.stations[serviceStations[qIdx]].name
        //         for (k in 0 until numClasses) {
        //             val qlen = totalQueueTime[qIdx][k] / simTime
        //             val basBlocked = basBlockedAtDest[qIdx][k]
        //             println("  $stationName class $k: qlen=$qlen, basBlockedAtDest=$basBlocked, currentQueueLength=${currentQueueLength[qIdx][k]}")
        //         }
        //     }
        // }

        // Extract results for each service node (Queue or Delay)
        for ((svcIdx, serviceStation) in serviceStations.withIndex()) {
            for (k in 0 until numClasses) {
                // Skip closed classes with 0 initial population unless they can receive jobs via class-switching
                val njobs = sn.njobs.get(k)
                val isClosedWithZeroPopulation = !isOpenClass[k] && java.lang.Double.isFinite(njobs) && njobs == 0.0
                val canReceiveViaClassSwitch = classSwitchClasses.contains(k)
                if ((!isClosedWithZeroPopulation || canReceiveViaClassSwitch) && mus[svcIdx][k] < Double.MAX_VALUE) {
                    QN[serviceStation, k] = getAvgQueueLength(svcIdx, k)
                    UN[serviceStation, k] = getUtilization(svcIdx, k)
                    RN[serviceStation, k] = getAvgResponseTime(svcIdx, k)
                    TN[serviceStation, k] = getThroughput(svcIdx, k)
                    AN[serviceStation, k] = getArrivalRate(svcIdx, k)
                    TardN[serviceStation, k] = if (tardinessTally[svcIdx][k].numberObs() > 0) {
                        tardinessTally[svcIdx][k].average()
                    } else {
                        0.0
                    }
                }
            }
        }

        // Extract results for Place nodes (Petri nets)
        val simTime = getActualSimTime()
        for ((placeListIdx, placeNodeIdx) in placeNodes.withIndex()) {
            val placeStationIdx = sn.nodeToStation[placeNodeIdx].toInt()
            if (placeStationIdx >= 0 && placeStationIdx < numStations) {
                // Finalize time-weighted statistics
                val currentTime = Sim.time()
                val elapsed = currentTime - lastPlaceUpdateTime[placeListIdx]
                if (elapsed > 0) {
                    for (k in 0 until numClasses) {
                        totalPlaceTokenTime[placeListIdx][k] += placeTokens[placeListIdx][k] * elapsed
                    }
                }

                // Finalize transit token time (tokens in transit that originated from this place)
                val transitElapsed = currentTime - lastPlaceTransitUpdateTime[placeListIdx]
                if (transitElapsed > 0) {
                    for (k in 0 until numClasses) {
                        placeTransitTokenTime[placeListIdx][k] += placeTokensInTransit[placeListIdx][k] * transitElapsed
                    }
                }

                for (k in 0 until numClasses) {
                    // Average tokens in place (queue length equivalent)
                    // Include tokens at the place AND tokens in transit that originated from this place
                    val qlen = if (simTime > 0) {
                        (totalPlaceTokenTime[placeListIdx][k] + placeTransitTokenTime[placeListIdx][k]) / simTime
                    } else {
                        0.0
                    }
                    QN[placeStationIdx, k] = qlen

                    // Throughput: Use MSER-5 truncated throughput if available
                    val tput = if (mserEnabled && ::placeCompletionObservations.isInitialized &&
                        placeListIdx < placeCompletionObservations.size &&
                        k < placeCompletionObservations[placeListIdx].size &&
                        ::observationTimes.isInitialized) {
                        val compObs = placeCompletionObservations[placeListIdx][k]
                        val truncationIdx = mserTruncationBatch * effectiveMserBatchSize
                        if (truncationIdx < compObs.size && truncationIdx < observationTimes.size) {
                            val startCompletions = compObs[truncationIdx]
                            val endCompletions = placeCompletions[placeListIdx][k]
                            val startTime = observationTimes[truncationIdx]
                            val endTime = Sim.time()
                            val elapsed = endTime - startTime
                            if (elapsed > 0) {
                                (endCompletions - startCompletions).toDouble() / elapsed
                            } else {
                                0.0
                            }
                        } else if (simTime > 0) {
                            // Fallback if truncation data not available
                            placeCompletions[placeListIdx][k].toDouble() / simTime
                        } else {
                            0.0
                        }
                    } else if (simTime > 0) {
                        // Non-MSER mode
                        placeCompletions[placeListIdx][k].toDouble() / simTime
                    } else {
                        0.0
                    }
                    TN[placeStationIdx, k] = tput
                }
            }
        }

        // Extract results for Join nodes
        for ((joinListIdx, joinStationIdx) in joinStations.withIndex()) {
            if (joinStationIdx >= 0 && joinStationIdx < numStations) {
                // Finalize time-weighted queue length statistics
                val currentTime = Sim.time()
                for (k in 0 until numClasses) {
                    val elapsed = currentTime - lastJoinUpdateTime[joinListIdx][k]
                    if (elapsed > 0) {
                        totalJoinQueueTime[joinListIdx][k] += currentJoinQueueLength[joinListIdx][k] * elapsed
                    }

                    // Average queue length at Join
                    val qlen = if (simTime > 0) {
                        totalJoinQueueTime[joinListIdx][k] / simTime
                    } else {
                        0.0
                    }
                    QN[joinStationIdx, k] = qlen

                    // Throughput: completed synchronizations / time
                    val tput = if (simTime > 0) {
                        joinCompletions[joinListIdx][k].toDouble() / simTime
                    } else {
                        0.0
                    }
                    TN[joinStationIdx, k] = tput

                    // Response time: average from Tally
                    val respTime = if (joinResponseTimeTally[joinListIdx][k].numberObs() > 0) {
                        joinResponseTimeTally[joinListIdx][k].average()
                    } else {
                        0.0
                    }
                    RN[joinStationIdx, k] = respTime

                    // Utilization: 0 for Join (no servers)
                    UN[joinStationIdx, k] = 0.0

                    // Arrival rate: arrived forked jobs / time
                    val arvR = if (simTime > 0) {
                        arrivedAtJoin[joinListIdx][k].toDouble() / simTime
                    } else {
                        0.0
                    }
                    AN[joinStationIdx, k] = arvR
                }
            }
        }

        // Extract results for each source
        for ((srcIdx, sourceStation) in sourceStations.withIndex()) {
            for (k in 0 until numClasses) {
                if (lambdas[srcIdx][k] > 0) {
                    TN[sourceStation, k] = lambdas[srcIdx][k]
                }
            }
        }

        // System metrics
        for (k in 0 until numClasses) {
            if (isOpenClass[k]) {
                // Open class: use arrival rates to check if class is active
                var totalArrival = 0.0
                for (srcIdx in sourceStations.indices) {
                    totalArrival += lambdas[srcIdx][k]
                }
                if (totalArrival > 0) {
                    XN[0, k] = getSystemThroughput(k)
                    CN[0, k] = getSystemResponseTime(k)
                    SysTardN[0, k] = if (systemTardinessTally[k].numberObs() > 0) {
                        systemTardinessTally[k].average()
                    } else {
                        0.0
                    }
                }
            } else {
                // Closed class: use throughput at any station (population is conserved)
                var foundTput = false
                for (svcIdx in 0 until numServiceNodes) {
                    val tput = getThroughput(svcIdx, k)
                    if (tput > 0) {
                        XN[0, k] = tput
                        CN[0, k] = getSystemResponseTime(k)
                        SysTardN[0, k] = if (systemTardinessTally[k].numberObs() > 0) {
                            systemTardinessTally[k].average()
                        } else {
                            0.0
                        }
                        foundTput = true
                        break
                    }
                }
                // If no service node throughput, check Place nodes (Petri net)
                if (!foundTput && placeNodes.isNotEmpty()) {
                    for (placeListIdx in placeNodes.indices) {
                        // Use MSER-5 truncated throughput for Place nodes
                        val tput = if (mserEnabled && ::placeCompletionObservations.isInitialized &&
                            placeListIdx < placeCompletionObservations.size &&
                            k < placeCompletionObservations[placeListIdx].size &&
                            ::observationTimes.isInitialized) {
                            val compObs = placeCompletionObservations[placeListIdx][k]
                            val truncationIdx = mserTruncationBatch * effectiveMserBatchSize
                            if (truncationIdx < compObs.size && truncationIdx < observationTimes.size) {
                                val startCompletions = compObs[truncationIdx]
                                val endCompletions = placeCompletions[placeListIdx][k]
                                val startTime = observationTimes[truncationIdx]
                                val endTime = Sim.time()
                                val elapsed = endTime - startTime
                                if (elapsed > 0) {
                                    (endCompletions - startCompletions).toDouble() / elapsed
                                } else {
                                    0.0
                                }
                            } else if (simTime > 0) {
                                placeCompletions[placeListIdx][k].toDouble() / simTime
                            } else {
                                0.0
                            }
                        } else if (simTime > 0) {
                            placeCompletions[placeListIdx][k].toDouble() / simTime
                        } else {
                            0.0
                        }
                        if (tput > 0) {
                            XN[0, k] = tput
                            break
                        }
                    }
                }
            }
        }

        // Apply control variates correction if enabled
        if (useControlVariates) {
            applyControlVariateCorrection(QN, UN, RN, TN)
        }

        // Compute arrival rates from throughputs using routing matrix (aligns with JMT semantics)
        // Formula: AN[ist, k] = sum over (jst, r) of TN[jst, r] * rt[(jst*R)+r, (ist*R)+k]
        val computedAN = Matrix(numStations, numClasses)
        for (ist in 0 until numStations) {
            val istNode = sn.stationToNode[ist].toInt()
            if (sn.nodetype[istNode] == NodeType.Source) continue  // Source has no arrivals
            for (k in 0 until numClasses) {
                var a = 0.0
                for (jst in 0 until numStations) {
                    for (r in 0 until numClasses) {
                        a += TN[jst, r] * sn.rt[(jst * numClasses) + r, (ist * numClasses) + k]
                    }
                }
                computedAN[ist, k] = a
            }
        }

        val result = DESResult()
        result.QN = QN
        result.UN = UN
        result.RN = RN
        result.TN = TN
        result.AN = computedAN
        result.TardN = TardN
        result.SysTardN = SysTardN
        result.CN = CN
        result.XN = XN
        result.sn = sn

        // Compute FCR metrics if there are any finite capacity regions
        if (numRegions > 0) {
            val simTime = getActualSimTime()

            // Finalize time-weighted statistics for all regions
            for (f in 0 until numRegions) {
                updateRegionTimeWeightedStats(f)
            }

            result.QNfcr = Matrix(numRegions, numClasses)
            result.UNfcr = Matrix(numRegions, numClasses)
            result.RNfcr = Matrix(numRegions, numClasses)
            result.TNfcr = Matrix(numRegions, numClasses)
            result.ANfcr = Matrix(numRegions, numClasses)
            result.WNfcr = Matrix(numRegions, numClasses)

            for (f in 0 until numRegions) {
                for (k in 0 until numClasses) {
                    // Average queue length in region (time-weighted)
                    val qlen = if (simTime > 0) {
                        totalRegionJobTime[f][k] / simTime
                    } else {
                        0.0
                    }
                    result.QNfcr!![f, k] = qlen

                    // Utilization: NaN (not meaningful for regions)
                    result.UNfcr!![f, k] = Double.NaN

                    // Throughput: departures per time
                    val tput = if (simTime > 0) {
                        regionCompletions[f][k].toDouble() / simTime
                    } else {
                        0.0
                    }
                    result.TNfcr!![f, k] = tput

                    // Response time: Q/T (Little's law)
                    val respTime = if (tput > 0) {
                        qlen / tput
                    } else {
                        0.0
                    }
                    result.RNfcr!![f, k] = respTime
                    result.WNfcr!![f, k] = respTime  // Residence time = response time for regions

                    // Arrival rate: computed from inter-arrival times (corrected implementation)
                    // In JMT, FCR_ARRIVAL_RATE uses InverseMeasure which samples currentTime (buggy)
                    // We correctly compute arrival rate = 1 / mean(inter-arrival time)
                    val arrivalRate = if (regionArrivalCount[f][k] > 1) {
                        val meanInterArrival = regionInterArrivalTimeSum[f][k] / (regionArrivalCount[f][k] - 1)
                        if (meanInterArrival > 0) 1.0 / meanInterArrival else 0.0
                    } else {
                        0.0
                    }
                    result.ANfcr!![f, k] = arrivalRate
                }
            }
        }

        // Set variance reduction metadata
        result.varianceReductionMethod = options.config?.variates ?: "none"

        // Transfer convergence data
        result.converged = hasConverged
        result.stoppingReason = stoppingReason
        result.convergenceBatches = if (::queueBatchMeans.isInitialized && queueBatchMeans.isNotEmpty() && queueBatchMeans[0].isNotEmpty()) {
            queueBatchMeans[0][0].size
        } else {
            0
        }

        // Transfer CI half-widths
        val QNCI = Matrix(numStations, numClasses)
        QNCI.fill(0.0)
        val UNCI = Matrix(numStations, numClasses)
        UNCI.fill(0.0)
        val RNCI = Matrix(numStations, numClasses)
        RNCI.fill(0.0)
        val TNCI = Matrix(numStations, numClasses)
        TNCI.fill(0.0)

        // Transfer relative precision
        val QNRelPrec = Matrix(numStations, numClasses)
        QNRelPrec.fill(0.0)
        val UNRelPrec = Matrix(numStations, numClasses)
        UNRelPrec.fill(0.0)
        val RNRelPrec = Matrix(numStations, numClasses)
        RNRelPrec.fill(0.0)
        val TNRelPrec = Matrix(numStations, numClasses)
        TNRelPrec.fill(0.0)

        if (convergenceEnabled && ::finalQNCI.isInitialized) {
            for ((svcIdx, serviceStation) in serviceStations.withIndex()) {
                for (k in 0 until numClasses) {
                    QNCI[serviceStation, k] = finalQNCI[svcIdx][k]
                    UNCI[serviceStation, k] = finalUNCI[svcIdx][k]
                    RNCI[serviceStation, k] = finalRNCI[svcIdx][k]
                    TNCI[serviceStation, k] = finalTNCI[svcIdx][k]

                    QNRelPrec[serviceStation, k] = finalQNRelPrec[svcIdx][k]
                    UNRelPrec[serviceStation, k] = finalUNRelPrec[svcIdx][k]
                    RNRelPrec[serviceStation, k] = finalRNRelPrec[svcIdx][k]
                    TNRelPrec[serviceStation, k] = finalTNRelPrec[svcIdx][k]
                }
            }
        }

        // Compute OBM confidence intervals if enabled
        if (options.confint > 0 && mserEnabled) {
            val ciResults = computeOBMConfidenceIntervals(options.confint)
            result.QNCI = ciResults.QNCI
            result.UNCI = ciResults.UNCI
            result.RNCI = ciResults.RNCI
            result.TNCI = ciResults.TNCI
            result.ANCI = ciResults.ANCI
            result.WNCI = ciResults.WNCI
        } else {
            result.QNCI = QNCI
            result.UNCI = UNCI
            result.RNCI = RNCI
            result.TNCI = TNCI
        }
        result.QNRelPrec = QNRelPrec
        result.UNRelPrec = UNRelPrec
        result.RNRelPrec = RNRelPrec
        result.TNRelPrec = TNRelPrec

        // Compute impatience statistics
        computeImpatienceStatistics(result, simTime)

        // Store response time samples for CDF computation
        setRespTimeSamples(result, numStations, numClasses)

        return result
    }

    /**
     * Helper to set response time samples array on DESResult, bypassing Kotlin platform type checks.
     */
    private fun setRespTimeSamples(result: DESResult, numStations: Int, numClasses: Int) {
        // Create outer array using reflection
        val outerArray = java.lang.reflect.Array.newInstance(
            java.lang.reflect.Array.newInstance(java.util.List::class.java, 0).javaClass,
            numStations
        )

        for (ist in 0 until numStations) {
            // Create inner array using reflection
            val innerArray = java.lang.reflect.Array.newInstance(java.util.List::class.java, numClasses)

            for (k in 0 until numClasses) {
                val svcIdx = serviceStations.indexOf(ist)
                val samples: Any = if (svcIdx >= 0 && svcIdx < responseTimeSamples.size && k < responseTimeSamples[svcIdx].size) {
                    val list = java.util.ArrayList<Any>()
                    for (sample in responseTimeSamples[svcIdx][k]) {
                        list.add(sample)
                    }
                    list
                } else {
                    java.util.ArrayList<Any>()
                }
                java.lang.reflect.Array.set(innerArray, k, samples)
            }
            java.lang.reflect.Array.set(outerArray, ist, innerArray)
        }

        // Set via reflection to bypass Kotlin type checks
        val field = DESResult::class.java.getField("respTimeSamples")
        field.set(result, outerArray)
    }

    /**
     * Compute impatience (reneging, balking, retrial) statistics and populate DESResult.
     */
    private fun computeImpatienceStatistics(result: DESResult, simTime: Double) {
        val renegedMatrix = Matrix(numStations, numClasses)
        val avgRenegingWaitMatrix = Matrix(numStations, numClasses)
        val renegingRateMatrix = Matrix(numStations, numClasses)
        val balkedMatrix = Matrix(numStations, numClasses)
        val balkingProbMatrix = Matrix(numStations, numClasses)
        val retriedMatrix = Matrix(numStations, numClasses)
        val retrialDroppedMatrix = Matrix(numStations, numClasses)
        val avgOrbitMatrix = Matrix(numStations, numClasses)

        for (svcIdx in 0 until numServiceNodes) {
            val ist = serviceStations[svcIdx]
            for (k in 0 until numClasses) {
                // Reneging statistics
                val reneged = renegedCustomers[svcIdx][k]
                renegedMatrix.set(ist, k, reneged.toDouble())

                val avgRenegingWait = if (reneged > 0) {
                    totalRenegingWaitTime[svcIdx][k] / reneged
                } else 0.0
                avgRenegingWaitMatrix.set(ist, k, avgRenegingWait)

                val completed = completedCustomers[svcIdx][k]
                val dropped = droppedCustomers[svcIdx][k]
                val total = completed + reneged + dropped
                val renegingRate = if (total > 0) reneged.toDouble() / total else 0.0
                renegingRateMatrix.set(ist, k, renegingRate)

                // Balking statistics
                val balked = balkedCustomers[svcIdx][k]
                balkedMatrix.set(ist, k, balked.toDouble())

                val arrived = arrivedCustomers[svcIdx][k]
                val balkProb = if (arrived > 0) balked.toDouble() / arrived else 0.0
                balkingProbMatrix.set(ist, k, balkProb)

                // Retrial statistics
                val retried = retriedCustomers[svcIdx][k]
                retriedMatrix.set(ist, k, retried.toDouble())

                val retrialDrop = maxRetriesExceeded[svcIdx][k]
                retrialDroppedMatrix.set(ist, k, retrialDrop.toDouble())

                // Average orbit size (time-weighted)
                updateOrbitTimeStats(svcIdx, k)
                val avgOrbit = if (simTime > 0) totalOrbitTime[svcIdx][k] / simTime else 0.0
                avgOrbitMatrix.set(ist, k, avgOrbit)
            }
        }

        result.renegedCustomers = renegedMatrix
        result.avgRenegingWaitTime = avgRenegingWaitMatrix
        result.renegingRate = renegingRateMatrix
        result.balkedCustomers = balkedMatrix
        result.balkingProbability = balkingProbMatrix
        result.retriedCustomers = retriedMatrix
        result.retrialDropped = retrialDroppedMatrix
        result.avgOrbitSize = avgOrbitMatrix
    }

    /**
     * Apply control variates correction to performance metrics.
     * Uses the deviation of sampled service times from theoretical means to adjust estimates.
     *
     * The control variate correction is based on the relationship:
     * U = λ * E[S], so if service times were biased, utilization is biased proportionally.
     * The correction factor is E[S]_theoretical / E[S]_observed.
     */
    private fun applyControlVariateCorrection(QN: Matrix, UN: Matrix, RN: Matrix, TN: Matrix) {
        for ((svcIdx, serviceStation) in serviceStations.withIndex()) {
            for (k in 0 until numClasses) {
                val count = serviceSampleCount[svcIdx][k]
                val expectedMean = serviceExpectedMean[svcIdx][k]

                if (count > 0 && expectedMean > 0) {
                    val actualMean = serviceSampleSum[svcIdx][k] / count.toDouble()

                    if (actualMean > 0 && actualMean.isFinite()) {
                        val correctionFactor = expectedMean / actualMean

                        // Only apply correction if factor is reasonable (within 50%)
                        if (correctionFactor > 0.5 && correctionFactor < 2.0) {
                            // Correct utilization: U = λ * E[S], so if E[S] biased, U biased
                            val correctedU = UN[serviceStation, k] * correctionFactor
                            UN[serviceStation, k] = correctedU.coerceIn(0.0, 1.0)

                            // Response time also scales with service time deviation
                            // R = W + E[S] where W is waiting time
                            // If E[S] was biased high, R was biased high proportionally
                            val originalR = RN[serviceStation, k]
                            if (originalR > 0) {
                                // Estimate waiting time contribution: W = R - E[S]_observed
                                val waitingTime = (originalR - actualMean).coerceAtLeast(0.0)
                                // Corrected R = W + E[S]_theoretical
                                RN[serviceStation, k] = waitingTime + expectedMean
                            }
                        }
                    }
                }
            }
        }
    }

    // Transient analysis data structures
    private var transientMode = false
    private lateinit var transientTimes: MutableList<Double>
    private lateinit var transientQueueLengths: Array<Array<MutableList<Double>>>
    private lateinit var transientUtilizations: Array<Array<MutableList<Double>>>
    private lateinit var transientThroughputs: Array<Array<MutableList<Double>>>
    private lateinit var transientCompletions: Array<Array<MutableList<Int>>>
    private var transientSamplingInterval = 0.0
    private var lastTransientSampleTime = 0.0
    private lateinit var lastTransientQueueTime: Array<DoubleArray>
    private lateinit var lastTransientBusyTime: Array<DoubleArray>
    private lateinit var lastTransientCompletions: Array<IntArray>

    /**
     * Run transient simulation for the specified time horizon.
     * No warmup period - collects metrics from time 0.
     * MSER-5 is disabled during transient analysis.
     */
    fun simulateTransient(timeHorizon: Double) {
        transientMode = true
        mserEnabled = false  // Disable MSER-5 for transient analysis
        Sim.init()

        // Warn that timespan is used instead of samples for transient analysis
        if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
            println("DES transient: using options.timespan=[${options.timespan[0]}, ${options.timespan[1]}] (options.samples ignored)")
        }

        // Disable event-count based stopping for transient mode
        // Transient simulation uses time-based stopping via EndOfSimulation event at timeHorizon
        maxEvents = Long.MAX_VALUE
        totalEventCount = 0

        // Initialize trace logging if DEBUG verbose level
        initializeTracing()

        // Initialize Logger nodes for CSV output
        initializeLoggers()

        // Initialize routing strategies
        initializeRouting()

        // Initialize random generators
        initializeGenerators()

        // Initialize cache samplers (requires routingRng to be initialized)
        initializeCacheSamplers()

        // Initialize queues and statistics
        initializeState()

        // Validate closed class routing (no routes to Sink)
        validateClosedClassRouting()

        // Initialize queue populations - use init_sol if provided, otherwise reference stations
        if (initSol != null && !initSol.isEmpty) {
            initFromInitSol()
        } else {
            initClosedClassPopulations()
        }

        // Initialize transient data collection
        initializeTransient(timeHorizon)

        // Initialize convergence checking structures (but disable for transient mode)
        // This ensures batch arrays are initialized even though convergence won't be checked
        initializeConvergence(Long.MAX_VALUE)
        convergenceEnabled = false  // Force disable for transient mode

        // Mark warmup as done immediately (no warmup for transient)
        warmupDone = true
        warmupEndTime = 0.0

        // Schedule periodic sampling for transient metrics
        TransientSampleEvent().schedule(transientSamplingInterval)

        // Schedule end of simulation
        EndOfSimulation().schedule(timeHorizon)

        // Schedule initial arrivals from each source for OPEN classes only
        // Closed classes have fixed population already initialized at reference stations
        for (srcIdx in 0 until numSources) {
            for (k in 0 until numClasses) {
                if (isOpenClass[k]) {
                    val interarrivalTime = generateInterarrivalTime(srcIdx, k)
                    if (interarrivalTime > 0) {
                        ExternalArrival(srcIdx, k).schedule(interarrivalTime)
                    }
                }
            }
        }

        // Run simulation
        Sim.start()
    }

    /**
     * Validate that closed classes do not route to Sink nodes.
     * Closed class jobs must circulate perpetually to conserve population.
     */
    private fun validateClosedClassRouting() {
        for (k in 0 until numClasses) {
            if (isClosedClass[k]) {
                // Check that no node routes this class to a Sink
                for (fromNode in 0 until numNodes) {
                    for (sinkNode in sinkNodes) {
                        val fromIdx = fromNode * numClasses + k
                        val toIdx = sinkNode * numClasses + k
                        if (fromIdx < sn.rtnodes.numRows && toIdx < sn.rtnodes.numCols) {
                            val prob = sn.rtnodes.get(fromIdx, toIdx)
                            if (prob > 0) {
                                throw RuntimeException(
                                    "solver_ssj: Closed class $k has routing to Sink " +
                                    "(violates population conservation)"
                                )
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Initialize closed class populations by injecting jobs at their reference stations.
     * Called at simulation start (time 0) before any events are processed.
     */
    private fun initClosedClassPopulations() {
        for (k in 0 until numClasses) {
            if (isClosedClass[k]) {
                val population = closedClassPopulation[k]
                val refStationIdx = referenceStation[k]

                // Find the service node index for this reference station
                val refNodeIdx = sn.stationToNode[refStationIdx].toInt()
                val nodeType = sn.nodetype[refNodeIdx]

                // Check if reference station is a Queue/Delay or a Place node
                if (nodeType == NodeType.Place) {
                    // For Petri nets: tokens are already initialized in placeTokens array
                    // via initializePlaceTransitionParams() which reads from sn.state
                    // Just check the transitions are enabled and schedule any immediate firings
                    if (placeNodes.isNotEmpty() && transitionNodes.isNotEmpty()) {
                        checkAndFireTransitions()
                    }
                    // Continue to next class - population already handled
                    continue
                }

                val queueIdx = serviceNodes.indexOf(refNodeIdx)

                if (queueIdx < 0) {
                    throw RuntimeException(
                        "solver_ssj: Reference station $refStationIdx for closed class $k " +
                        "is not a service node (Queue, Delay, or Place)"
                    )
                }

                // Inject all jobs at reference station at time 0
                // Assign job IDs to classes that expect replies (for synchronous call tracking)
                val needsJobId = synchCallReplyClass[k] >= 0
                for (jobIdx in 0 until population) {
                    val customer = Customer(
                        classId = k,
                        priority = classPrio[k],
                        systemArrivalTime = 0.0,
                        queueArrivalTime = 0.0,
                        randomRank = siroRng.nextDouble(),
                        jobId = if (needsJobId) nextJobId++ else -1L,
                        absoluteDeadline = 0.0 + classDeadline[k]
                    )
                    arriveAtQueue(queueIdx, customer)
                }
            }
        }
    }

    /**
     * Initialize queue populations from init_sol matrix.
     * Called instead of initClosedClassPopulations() when init_sol is provided.
     *
     * The init_sol matrix is [1 x (M*K)] in row-major order where:
     * - M = number of stations (numStations)
     * - K = number of classes (numClasses)
     * - Element [stationIdx * K + classIdx] = initial queue length at (station, class)
     *
     * For closed classes, validates that total population matches njobs[k].
     */
    private fun initFromInitSol() {
        if (initSol == null || initSol.isEmpty) {
            // Fallback to standard initialization
            initClosedClassPopulations()
            return
        }

        // Parse init_sol: row-major [station0_class0, station0_class1, ..., stationM-1_classK-1]
        var idx = 0
        for (stationIdx in 0 until numStations) {
            for (k in 0 until numClasses) {
                val initialQueueLength = if (idx < initSol.length()) {
                    FastMath.max(0.0, initSol.get(idx)).toInt()
                } else {
                    0
                }
                idx++

                if (initialQueueLength <= 0) continue

                // Find the service node index for this station
                val nodeIdx = sn.stationToNode[stationIdx].toInt()
                val nodeType = sn.nodetype[nodeIdx]

                // Skip non-service nodes (Source, Sink, Router, etc.)
                if (nodeType != NodeType.Queue && nodeType != NodeType.Delay) {
                    continue
                }

                val queueIdx = serviceNodes.indexOf(nodeIdx)
                if (queueIdx < 0) {
                    // This station is not a service node in our list
                    continue
                }

                // Inject jobs at this station at time 0
                // Assign job IDs to classes that expect replies (for synchronous call tracking)
                val needsJobId = synchCallReplyClass[k] >= 0
                for (jobIdx in 0 until initialQueueLength) {
                    val customer = Customer(
                        classId = k,
                        priority = classPrio[k],
                        systemArrivalTime = 0.0,
                        queueArrivalTime = 0.0,
                        randomRank = siroRng.nextDouble(),
                        jobId = if (needsJobId) nextJobId++ else -1L,
                        absoluteDeadline = 0.0 + classDeadline[k]
                    )
                    arriveAtQueue(queueIdx, customer)
                }
            }
        }

        // For closed classes, verify population conservation
        for (k in 0 until numClasses) {
            if (isClosedClass[k]) {
                var totalInSystem = 0
                for (queueIdx in 0 until numServiceNodes) {
                    totalInSystem += currentQueueLength[queueIdx][k]
                }
                if (totalInSystem != closedClassPopulation[k]) {
                    throw RuntimeException(
                        "init_sol population mismatch for closed class $k: " +
                        "expected ${closedClassPopulation[k]}, got $totalInSystem"
                    )
                }
            }
        }
    }

    /**
     * Initialize transient data collection structures.
     */
    private fun initializeTransient(timeHorizon: Double) {
        // Target ~100 time points for smooth curves
        val targetObservations = 100
        transientSamplingInterval = timeHorizon / targetObservations

        transientTimes = mutableListOf()
        transientQueueLengths = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Double>() }
        }
        transientUtilizations = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Double>() }
        }
        transientThroughputs = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Double>() }
        }
        transientCompletions = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Int>() }
        }

        lastTransientSampleTime = 0.0
        lastTransientQueueTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        lastTransientBusyTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        lastTransientCompletions = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        // Record initial state at time 0 - use actual queue lengths (may be non-zero from init_sol)
        transientTimes.add(0.0)
        for (qIdx in 0 until numServiceNodes) {
            for (k in 0 until numClasses) {
                // Use actual current queue length (may have been set by initFromInitSol)
                transientQueueLengths[qIdx][k].add(currentQueueLength[qIdx][k].toDouble())
                // Initial per-server utilization based on busy servers (λ × S / c convention, consistent with MVA)
                val busyServers = currentBusyServers[qIdx][k]
                val initialUtil = if (isDelayNode[qIdx]) {
                    currentQueueLength[qIdx][k].toDouble()  // For delay nodes (c=inf), utilization = queue length
                } else {
                    busyServers.toDouble() / numServers[qIdx]
                }
                transientUtilizations[qIdx][k].add(initialUtil)
                transientThroughputs[qIdx][k].add(0.0)  // No throughput at t=0
                transientCompletions[qIdx][k].add(0)
            }
        }
    }

    /**
     * Transient sampling event - records metrics at regular intervals.
     */
    private inner class TransientSampleEvent : Event() {
        override fun actions() {
            val currentTime = Sim.time()
            val intervalDuration = currentTime - lastTransientSampleTime

            // Update queue stats before sampling
            for (qIdx in 0 until numServiceNodes) {
                for (k in 0 until numClasses) {
                    updateQueueStats(qIdx, k)
                    updateBusyStats(qIdx, k)
                }
            }

            transientTimes.add(currentTime)

            // Record time-weighted metrics over this interval
            for (qIdx in 0 until numServiceNodes) {
                for (k in 0 until numClasses) {
                    // Queue length: time-weighted average over interval
                    if (intervalDuration > 0) {
                        val intervalQueueTime = totalQueueTime[qIdx][k] - lastTransientQueueTime[qIdx][k]
                        val avgQueueLength = intervalQueueTime / intervalDuration
                        transientQueueLengths[qIdx][k].add(avgQueueLength)
                        lastTransientQueueTime[qIdx][k] = totalQueueTime[qIdx][k]

                        // Utilization: time-weighted average over interval (per-server)
                        val intervalBusyTime = totalBusyTime[qIdx][k] - lastTransientBusyTime[qIdx][k]
                        val avgUtilization = if (isDelayNode[qIdx]) {
                            avgQueueLength  // For delay nodes, utilization = queue length
                        } else {
                            intervalBusyTime / (intervalDuration * numServers[qIdx])
                        }
                        transientUtilizations[qIdx][k].add(avgUtilization)
                        lastTransientBusyTime[qIdx][k] = totalBusyTime[qIdx][k]

                        // Throughput: rate over interval
                        val intervalCompletions = completedCustomers[qIdx][k] - lastTransientCompletions[qIdx][k]
                        val throughput = intervalCompletions.toDouble() / intervalDuration
                        transientThroughputs[qIdx][k].add(throughput)
                        lastTransientCompletions[qIdx][k] = completedCustomers[qIdx][k]
                    } else {
                        transientQueueLengths[qIdx][k].add(currentQueueLength[qIdx][k].toDouble())
                        transientUtilizations[qIdx][k].add(0.0)
                        transientThroughputs[qIdx][k].add(0.0)
                    }
                    transientCompletions[qIdx][k].add(completedCustomers[qIdx][k])
                }
            }
            lastTransientSampleTime = currentTime

            // Schedule next sample if not at end
            if (currentTime + transientSamplingInterval < Double.MAX_VALUE) {
                TransientSampleEvent().schedule(transientSamplingInterval)
            }
        }
    }

    /**
     * Get transient analysis results.
     */
    fun getTransientDESResult(): DESResult {
        val result = DESResult()
        result.sn = sn

        // Create time points matrix
        val numTimePoints = transientTimes.size
        result.t = Matrix(numTimePoints, 1)
        for (i in 0 until numTimePoints) {
            result.t[i, 0] = transientTimes[i]
        }

        // Create transient matrices [stations x classes] at each time point
        // Format: QNt[station][class] is a Matrix(numTimePoints, 2) with columns [value, time]
        result.QNt = Array(numStations) { Array(numClasses) { Matrix(0, 0) } }
        result.UNt = Array(numStations) { Array(numClasses) { Matrix(0, 0) } }
        result.TNt = Array(numStations) { Array(numClasses) { Matrix(0, 0) } }

        // Fill in results for service nodes
        for ((svcIdx, serviceStation) in serviceStations.withIndex()) {
            for (k in 0 until numClasses) {
                if (svcIdx < transientQueueLengths.size && k < transientQueueLengths[svcIdx].size) {
                    val qObs = transientQueueLengths[svcIdx][k]
                    val uObs = transientUtilizations[svcIdx][k]
                    val tObs = transientThroughputs[svcIdx][k]
                    val obsCount = minOf(qObs.size, transientTimes.size)

                    // Create matrices with [value, time] format
                    val qMatrix = Matrix(obsCount, 2)
                    val uMatrix = Matrix(obsCount, 2)
                    val tMatrix = Matrix(obsCount, 2)

                    for (i in 0 until obsCount) {
                        qMatrix[i, 0] = qObs[i]
                        qMatrix[i, 1] = transientTimes[i]
                        uMatrix[i, 0] = uObs[i]
                        uMatrix[i, 1] = transientTimes[i]
                        tMatrix[i, 0] = tObs[i]
                        tMatrix[i, 1] = transientTimes[i]
                    }

                    result.QNt[serviceStation][k] = qMatrix
                    result.UNt[serviceStation][k] = uMatrix
                    result.TNt[serviceStation][k] = tMatrix
                }
            }
        }

        // Also compute final steady-state estimates (last interval averages)
        val QN = Matrix(numStations, numClasses)
        val UN = Matrix(numStations, numClasses)
        val RN = Matrix(numStations, numClasses)
        val TN = Matrix(numStations, numClasses)
        val CN = Matrix(1, numClasses)
        val XN = Matrix(1, numClasses)

        // Build set of classes that can receive jobs via class-switching
        val classSwitchClasses2 = mutableSetOf<Int>()
        for (c in 0 until sn.nchains) {
            val classesInChain = mutableListOf<Int>()
            for (k in 0 until numClasses) {
                if (sn.chains.get(c, k) > 0) {
                    classesInChain.add(k)
                }
            }
            if (classesInChain.size > 1) {
                var chainHasJobs = false
                for (k in classesInChain) {
                    if (sn.njobs.get(k) > 0) {
                        chainHasJobs = true
                        break
                    }
                }
                if (chainHasJobs) {
                    classSwitchClasses2.addAll(classesInChain)
                }
            }
        }

        for ((svcIdx, serviceStation) in serviceStations.withIndex()) {
            for (k in 0 until numClasses) {
                // Skip closed classes with 0 initial population unless they can receive jobs via class-switching
                val njobs = sn.njobs.get(k)
                val isClosedWithZeroPopulation = !isOpenClass[k] && java.lang.Double.isFinite(njobs) && njobs == 0.0
                val canReceiveViaClassSwitch2 = classSwitchClasses2.contains(k)
                if ((!isClosedWithZeroPopulation || canReceiveViaClassSwitch2) && mus[svcIdx][k] < Double.MAX_VALUE) {
                    QN[serviceStation, k] = getAvgQueueLength(svcIdx, k)
                    UN[serviceStation, k] = getUtilization(svcIdx, k)
                    RN[serviceStation, k] = getAvgResponseTime(svcIdx, k)
                    TN[serviceStation, k] = getThroughput(svcIdx, k)
                }
            }
        }

        for ((srcIdx, sourceStation) in sourceStations.withIndex()) {
            for (k in 0 until numClasses) {
                if (lambdas[srcIdx][k] > 0) {
                    TN[sourceStation, k] = lambdas[srcIdx][k]
                }
            }
        }

        for (k in 0 until numClasses) {
            var totalArrival = 0.0
            for (srcIdx in sourceStations.indices) {
                totalArrival += lambdas[srcIdx][k]
            }
            if (totalArrival > 0) {
                XN[0, k] = getSystemThroughput(k)
                CN[0, k] = getSystemResponseTime(k)
            }
        }

        result.QN = QN
        result.UN = UN
        result.RN = RN
        result.TN = TN
        result.CN = CN
        result.XN = XN

        // Store response time samples for CDF computation
        setRespTimeSamples(result, numStations, numClasses)

        return result
    }

    /**
     * Initialize MSER-5 data collection structures.
     */
    private fun initializeMSER(timeHorizon: Double) {
        // Target ~1000 observations for good MSER-5 estimation (200 batches of 5)
        val targetObservations = 1000
        mserSamplingInterval = timeHorizon / targetObservations

        queueLengthObservations = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Double>() }
        }
        throughputObservations = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Int>() }
        }
        placeCompletionObservations = Array(placeNodes.size) {
            Array(numClasses) { mutableListOf<Int>() }
        }
        observationTimes = mutableListOf()

        // Initialize MSER tracking arrays
        lastMserQueueTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        lastMserSampleTime = 0.0
    }

    /**
     * Initialize MSER-5 data collection structures for event-based stopping.
     * @param maxEventCount Maximum number of events (service completions) to simulate
     */
    private fun initializeMSEREventBased(maxEventCount: Long) {
        // Target ~1000 observations for good MSER-5 estimation (200 batches of 5)
        val targetObservations = 1000
        mserEventInterval = maxEventCount / targetObservations
        if (mserEventInterval < 1) mserEventInterval = 1
        lastMserEventCount = 0

        queueLengthObservations = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Double>() }
        }
        throughputObservations = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Int>() }
        }
        placeCompletionObservations = Array(placeNodes.size) {
            Array(numClasses) { mutableListOf<Int>() }
        }
        observationTimes = mutableListOf()

        // Initialize MSER tracking arrays
        lastMserQueueTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        lastMserSampleTime = 0.0
    }

    /**
     * Initialize transient detection and CI configuration from DESOptions.
     * Called at the start of simulation.
     */
    private fun initializeTransientAndCIConfig() {
        val desOptions = options as? DESOptions

        // Transient detection configuration
        val tranfilter = desOptions?.tranfilter ?: "mser5"
        effectiveMserBatchSize = desOptions?.mserbatch ?: DESOptions.DEFAULT_MSER_BATCH
        effectiveWarmupFraction = desOptions?.warmupfrac ?: DESOptions.DEFAULT_WARMUP_FRAC

        // Derive mserEnabled from tranfilter
        mserEnabled = (tranfilter == "mser5")

        // CI configuration
        effectiveCiMethod = desOptions?.cimethod ?: "obm"
        effectiveObmOverlap = desOptions?.obmoverlap ?: DESOptions.DEFAULT_OBM_OVERLAP
        effectiveCiMinBatch = desOptions?.ciminbatch ?: DESOptions.DEFAULT_CI_MIN_BATCH
        effectiveCiMinObs = desOptions?.ciminobs ?: DESOptions.DEFAULT_CI_MIN_OBS

        if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
            println("DES: Transient filter = $tranfilter" +
                    if (tranfilter == "mser5") " (batch size = $effectiveMserBatchSize)"
                    else if (tranfilter == "fixed") " (warmup fraction = $effectiveWarmupFraction)"
                    else "")
            println("DES: CI method = $effectiveCiMethod" +
                    if (effectiveCiMethod == "obm") " (overlap = $effectiveObmOverlap)" else "")
        }
    }

    /**
     * Initialize convergence checking data structures.
     * @param maxEventCount Maximum number of events (service completions) to simulate
     */
    private fun initializeConvergence(maxEventCount: Long) {
        // Read convergence options (cast to DESOptions if available)
        val desOptions = options as? DESOptions
        convergenceEnabled = desOptions?.cnvgon ?: false
        convergenceTolerance = desOptions?.cnvgtol ?: 0.05
        convergenceMinBatches = desOptions?.cnvgbatch ?: 20

        // Auto-calculate check interval if not specified
        val configuredInterval = desOptions?.cnvgchk ?: 0
        convergenceCheckInterval = if (configuredInterval > 0) {
            configuredInterval.toLong()
        } else {
            maxOf(1L, maxEventCount / 50)
        }
        lastConvergenceCheckEventCount = 0

        // Initialize batch means arrays
        queueBatchMeans = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Double>() }
        }
        utilBatchMeans = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Double>() }
        }
        respTimeBatchMeans = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Double>() }
        }
        throughputBatchMeans = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Double>() }
        }

        // Initialize batch accumulators
        batchStartQueueTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        batchStartBusyTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        batchStartCompletions = Array(numServiceNodes) { IntArray(numClasses) { 0 } }
        currentBatchRespTimeSum = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        currentBatchRespTimeCount = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        currentBatchObservations = 0
        batchStartTime = 0.0
        hasConverged = false
        stoppingReason = "max_events"

        // Initialize final CI arrays
        finalQNCI = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        finalUNCI = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        finalRNCI = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        finalTNCI = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        finalQNRelPrec = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        finalUNRelPrec = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        finalRNRelPrec = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        finalTNRelPrec = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
    }

    /**
     * Event counter for progress reporting
     */
    private var lastProgressEventCount: Long = 0

    /**
     * Check if event count thresholds have been reached and take appropriate action.
     * Called after each service completion.
     * - Handles warmup completion (for non-MSER mode)
     * - Handles MSER sampling (for MSER mode)
     * - Handles convergence sampling and checking
     * - Handles progress reporting
     * - Handles simulation stop when max events reached or converged
     */
    private fun checkEventCountStop() {
        totalEventCount++

        // Handle warmup completion (non-MSER mode)
        if (!mserEnabled && !warmupDone && totalEventCount >= warmupEventThreshold) {
            resetStatistics()
            warmupDone = true
        }

        // Handle MSER sampling
        if (mserEnabled && (totalEventCount - lastMserEventCount) >= mserEventInterval) {
            collectMSERSample()
            lastMserEventCount = totalEventCount
        }

        // Handle convergence sampling and checking at configured intervals
        if (convergenceEnabled &&
            (totalEventCount - lastConvergenceCheckEventCount) >= convergenceCheckInterval) {
            // Collect convergence sample
            collectConvergenceSample()
            lastConvergenceCheckEventCount = totalEventCount

            // Check for convergence
            if (checkConvergence()) {
                hasConverged = true
                stoppingReason = "convergence"
                if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                    println("\nDES: Convergence detected at event $totalEventCount")
                }
                finishSimulation()
                return
            }
        }

        // Handle progress reporting
        val progressInterval = maxEvents / 50
        if (progressInterval > 0 && (totalEventCount - lastProgressEventCount) >= progressInterval) {
            if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                if (lastProgressEventCount == 0L) {
                    System.out.printf("DES events: %6d ", totalEventCount)
                    System.out.flush()
                } else {
                    System.out.printf("\b\b\b\b\b\b\b %6d", totalEventCount)
                    System.out.flush()
                }
            }
            lastProgressEventCount = totalEventCount
        }

        // Handle streaming if collector is active
        if (stream != null) {
            val streamOpts = stream.getOptions()
            val currentTime = Sim.time()
            val shouldStream = when (streamOpts.mode) {
                jline.streaming.StreamingOptions.StreamMode.SAMPLED ->
                    (totalEventCount - lastStreamEventCount) >= streamOpts.sampleFrequency
                jline.streaming.StreamingOptions.StreamMode.TIME_WINDOW ->
                    (currentTime - lastStreamTime) >= streamOpts.timeWindowSeconds
            }
            if (shouldStream) {
                pushStreamingMetrics(currentTime)
                lastStreamEventCount = totalEventCount
                lastStreamTime = currentTime
            }
        }

        // Check if max events reached
        if (totalEventCount >= maxEvents) {
            stoppingReason = "max_events"
            finishSimulation()
        }
    }

    /**
     * Push current queue state metrics to the streaming collector.
     * Builds a queue length matrix [station x class] from current state.
     */
    private fun pushStreamingMetrics(currentTime: Double) {
        if (stream == null) return

        // Build queue length matrix from current state
        val nir = Matrix(numServiceNodes, numClasses)
        for (qIdx in 0 until numServiceNodes) {
            for (k in 0 until numClasses) {
                nir.set(qIdx, k, currentQueueLength[qIdx][k].toDouble())
            }
        }

        // Calculate time delta since last push (or from simulation start)
        val dt = if (lastStreamTime > 0.0) currentTime - lastStreamTime else currentTime

        // Record state with the streaming collector
        stream.recordState(currentTime, dt, nir, null, null)
    }

    /**
     * Collect an MSER sample based on current state.
     */
    private fun collectMSERSample() {
        val currentTime = Sim.time()
        val intervalDuration = currentTime - lastMserSampleTime

        // Update queue stats before sampling
        for (qIdx in 0 until numServiceNodes) {
            for (k in 0 until numClasses) {
                updateQueueStats(qIdx, k)
            }
        }

        observationTimes.add(currentTime)

        // Record time-weighted average queue length over this interval
        for (qIdx in 0 until numServiceNodes) {
            for (k in 0 until numClasses) {
                if (intervalDuration > 0) {
                    val intervalQueueTime = totalQueueTime[qIdx][k] - lastMserQueueTime[qIdx][k]
                    val avgQueueLength = intervalQueueTime / intervalDuration
                    queueLengthObservations[qIdx][k].add(avgQueueLength)
                    lastMserQueueTime[qIdx][k] = totalQueueTime[qIdx][k]
                } else {
                    // Include blocked jobs (BAS, BBS, and FCR) trying to enter this queue
                    val effectiveQueueLength = currentQueueLength[qIdx][k] + basBlockedAtDest[qIdx][k] + bbsBlockedAtDest[qIdx][k] + fcrBlockedAtDest[qIdx][k]
                    queueLengthObservations[qIdx][k].add(effectiveQueueLength.toDouble())
                }
                throughputObservations[qIdx][k].add(completedCustomers[qIdx][k])
            }
        }

        // Record Place completion counts for MSER truncation
        for ((placeListIdx, _) in placeNodes.withIndex()) {
            for (k in 0 until numClasses) {
                placeCompletionObservations[placeListIdx][k].add(placeCompletions[placeListIdx][k])
            }
        }

        lastMserSampleTime = currentTime
    }

    /**
     * Collect a convergence sample based on current state.
     * Called at the same interval as MSER sampling.
     * Aggregates observations into batches and computes batch means.
     */
    private fun collectConvergenceSample() {
        if (!convergenceEnabled) return

        val currentTime = Sim.time()

        // Update stats before sampling
        for (qIdx in 0 until numServiceNodes) {
            if (isPSScheduling(schedStrategies[qIdx])) {
                updatePSBusyStats(qIdx)
            }
            for (k in 0 until numClasses) {
                updateQueueStats(qIdx, k)
                if (!isPSScheduling(schedStrategies[qIdx])) {
                    updateBusyStats(qIdx, k)
                }
            }
        }

        currentBatchObservations++

        // Check if batch is complete (use same batch size as MSER)
        if (currentBatchObservations >= effectiveMserBatchSize) {
            finalizeBatch(currentTime)
        }
    }

    /**
     * Finalize a batch - compute batch means for all metrics and add to lists.
     */
    private fun finalizeBatch(currentTime: Double) {
        val batchDuration = currentTime - batchStartTime
        if (batchDuration <= 0) {
            // Reset batch without recording
            resetBatchAccumulators(currentTime)
            return
        }

        for (qIdx in 0 until numServiceNodes) {
            for (k in 0 until numClasses) {
                // Queue length batch mean (time-weighted)
                val queueTimeDelta = totalQueueTime[qIdx][k] - batchStartQueueTime[qIdx][k]
                val avgQueueLength = queueTimeDelta / batchDuration
                queueBatchMeans[qIdx][k].add(avgQueueLength)

                // Per-server Utilization batch mean (time-weighted, λ × S / c convention consistent with MVA)
                val busyTimeDelta = totalBusyTime[qIdx][k] - batchStartBusyTime[qIdx][k]
                val avgUtil = if (isDelayNode[qIdx]) {
                    avgQueueLength  // For delay nodes (c=inf), utilization = queue length per LINE convention
                } else {
                    busyTimeDelta / (batchDuration * numServers[qIdx])
                }
                utilBatchMeans[qIdx][k].add(avgUtil)

                // Response time batch mean
                val respTimeCount = currentBatchRespTimeCount[qIdx][k]
                val avgRespTime = if (respTimeCount > 0) {
                    currentBatchRespTimeSum[qIdx][k] / respTimeCount
                } else {
                    // Fallback: use overall average if available
                    if (responseTimeTally[qIdx][k].numberObs() > 0) {
                        responseTimeTally[qIdx][k].average()
                    } else {
                        0.0
                    }
                }
                respTimeBatchMeans[qIdx][k].add(avgRespTime)

                // Throughput batch mean
                val completionsDelta = completedCustomers[qIdx][k] - batchStartCompletions[qIdx][k]
                val avgThroughput = completionsDelta.toDouble() / batchDuration
                throughputBatchMeans[qIdx][k].add(avgThroughput)
            }
        }

        // Reset batch accumulators for next batch
        resetBatchAccumulators(currentTime)
    }

    /**
     * Reset batch accumulators for the start of a new batch.
     */
    private fun resetBatchAccumulators(currentTime: Double) {
        batchStartTime = currentTime
        currentBatchObservations = 0

        for (qIdx in 0 until numServiceNodes) {
            for (k in 0 until numClasses) {
                batchStartQueueTime[qIdx][k] = totalQueueTime[qIdx][k]
                batchStartBusyTime[qIdx][k] = totalBusyTime[qIdx][k]
                batchStartCompletions[qIdx][k] = completedCustomers[qIdx][k]
                currentBatchRespTimeSum[qIdx][k] = 0.0
                currentBatchRespTimeCount[qIdx][k] = 0
            }
        }
    }

    /**
     * Record a response time observation for the current batch.
     * Called when a customer completes service.
     */
    private fun recordResponseTimeForBatch(qIdx: Int, k: Int, responseTime: Double) {
        if (!convergenceEnabled) return
        if (qIdx < numServiceNodes && k < numClasses) {
            currentBatchRespTimeSum[qIdx][k] += responseTime
            currentBatchRespTimeCount[qIdx][k]++
        }
    }

    /**
     * Finish the simulation - called when max events reached.
     */
    private fun finishSimulation() {
        // Final update of queue and busy statistics
        for (qIdx in 0 until numServiceNodes) {
            // For PS queues, use the PS-specific busy stats update
            if (isPSScheduling(schedStrategies[qIdx])) {
                updatePSBusyStats(qIdx)
                for (k in 0 until numClasses) {
                    updateQueueStats(qIdx, k)
                }
            } else {
                for (k in 0 until numClasses) {
                    updateQueueStats(qIdx, k)
                    updateBusyStats(qIdx, k)
                }
            }
        }

        // Apply MSER-5 truncation to determine warmup period
        if (mserEnabled) {
            applyMSER5Truncation()
        }

        // Compute final CI matrices for convergence results
        computeFinalCIMatrices()

        // Close trace writer
        closeTracing()

        // Close Logger file writers
        closeLoggers()

        if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
            // Print final event count before newline
            System.out.printf("\b\b\b\b\b\b\b %6d", totalEventCount)
            println()
        }

        Sim.stop()
    }

    /**
     * Check if all metrics have converged.
     * @return true if all metrics for all active station/class combinations have converged
     */
    private fun checkConvergence(): Boolean {
        if (!convergenceEnabled) return false

        // Check if we have enough batches
        val numBatches = if (queueBatchMeans.isNotEmpty() && queueBatchMeans[0].isNotEmpty()) {
            queueBatchMeans[0][0].size
        } else {
            0
        }

        if (numBatches < convergenceMinBatches) return false

        // Check all metrics for all active station/class combinations
        for (qIdx in 0 until numServiceNodes) {
            for (k in 0 until numClasses) {
                // Skip inactive combinations (disabled service rate)
                if (mus[qIdx][k] >= Double.MAX_VALUE) continue

                // Check queue length convergence
                if (!isMetricConverged(queueBatchMeans[qIdx][k], "Q[$qIdx][$k]")) return false

                // Check utilization convergence
                if (!isMetricConverged(utilBatchMeans[qIdx][k], "U[$qIdx][$k]")) return false

                // Check response time convergence (skip if no data)
                if (respTimeBatchMeans[qIdx][k].isNotEmpty() && respTimeBatchMeans[qIdx][k].any { it > 0 }) {
                    if (!isMetricConverged(respTimeBatchMeans[qIdx][k], "R[$qIdx][$k]")) return false
                }

                // Check throughput convergence
                if (!isMetricConverged(throughputBatchMeans[qIdx][k], "T[$qIdx][$k]")) return false
            }
        }

        return true
    }

    /**
     * Check if a single metric has converged based on batch means.
     * @param batchMeans list of batch mean values
     * @param metricName name of the metric (for debugging)
     * @return true if the metric has converged (CI half-width / mean < tolerance)
     */
    private fun isMetricConverged(batchMeans: List<Double>, metricName: String): Boolean {
        if (batchMeans.size < 2) return false

        // Compute mean
        val n = batchMeans.size
        var sum = 0.0
        for (x in batchMeans) sum += x
        val mean = sum / n

        // Skip if mean is effectively zero (treat as converged)
        if (kotlin.math.abs(mean) < 1e-12) return true

        // Compute sample variance
        var variance = 0.0
        for (x in batchMeans) {
            val diff = x - mean
            variance += diff * diff
        }
        variance /= (n - 1)

        // Standard error
        val stdErr = kotlin.math.sqrt(variance / n)

        // Get confidence level (default to 0.95 if not configured)
        val confintLevel = if (options.confint > 0) options.confint else 0.95
        val alpha = 1.0 - confintLevel

        // t-critical value
        val tCrit = getTCriticalValue(alpha, n - 1)

        // CI half-width
        val ciHalfWidth = tCrit * stdErr

        // Relative precision
        val relPrec = ciHalfWidth / kotlin.math.abs(mean)

        if (options.verbose == VerboseLevel.DEBUG) {
            System.out.println("Convergence check $metricName: mean=$mean, relPrec=$relPrec, tol=$convergenceTolerance")
        }

        return relPrec <= convergenceTolerance
    }

    /**
     * Get t-distribution critical value for given alpha and degrees of freedom.
     * Uses table lookup for common values, approximation otherwise.
     */
    private fun getTCriticalValue(alpha: Double, df: Int): Double {
        // Common t-critical values for two-tailed test
        val halfAlpha = alpha / 2.0

        // For large df (>30), t approaches normal distribution
        if (df > 30) {
            return when {
                halfAlpha <= 0.005 -> 2.576  // 99% CI
                halfAlpha <= 0.025 -> 1.96   // 95% CI
                halfAlpha <= 0.05 -> 1.645   // 90% CI
                else -> 1.28                  // 80% CI
            }
        }

        // Table values for smaller df (95% CI values)
        val t95 = mapOf(
            1 to 12.71, 2 to 4.30, 3 to 3.18, 4 to 2.78, 5 to 2.57,
            6 to 2.45, 7 to 2.36, 8 to 2.31, 9 to 2.26, 10 to 2.23,
            11 to 2.20, 12 to 2.18, 13 to 2.16, 14 to 2.14, 15 to 2.13,
            16 to 2.12, 17 to 2.11, 18 to 2.10, 19 to 2.09, 20 to 2.09,
            25 to 2.06, 30 to 2.04
        )

        // For 95% CI (most common case)
        if (halfAlpha in 0.02..0.03) {
            return t95[df] ?: t95[minOf(df, 30)] ?: 2.0
        }

        // Rough scaling for other confidence levels
        val baseT = t95[df] ?: t95[minOf(df, 30)] ?: 2.0
        return when {
            halfAlpha <= 0.005 -> baseT * 1.32  // 99% CI
            halfAlpha <= 0.05 -> baseT * 0.84   // 90% CI
            else -> baseT * 0.65                 // 80% CI
        }
    }

    /**
     * Compute final CI matrices for all metrics.
     * Called at end of simulation to populate CI half-widths and relative precision.
     */
    private fun computeFinalCIMatrices() {
        if (!convergenceEnabled) return

        val confintLevel = if (options.confint > 0) options.confint else 0.95
        val alpha = 1.0 - confintLevel

        for (qIdx in 0 until numServiceNodes) {
            for (k in 0 until numClasses) {
                // Queue length CI
                val qResult = computeCIForMetric(queueBatchMeans[qIdx][k], alpha)
                finalQNCI[qIdx][k] = qResult.first
                finalQNRelPrec[qIdx][k] = qResult.second

                // Utilization CI
                val uResult = computeCIForMetric(utilBatchMeans[qIdx][k], alpha)
                finalUNCI[qIdx][k] = uResult.first
                finalUNRelPrec[qIdx][k] = uResult.second

                // Response time CI
                val rResult = computeCIForMetric(respTimeBatchMeans[qIdx][k], alpha)
                finalRNCI[qIdx][k] = rResult.first
                finalRNRelPrec[qIdx][k] = rResult.second

                // Throughput CI
                val tResult = computeCIForMetric(throughputBatchMeans[qIdx][k], alpha)
                finalTNCI[qIdx][k] = tResult.first
                finalTNRelPrec[qIdx][k] = tResult.second
            }
        }
    }

    /**
     * Compute CI half-width and relative precision for a single metric.
     * @param batchMeans list of batch mean values
     * @param alpha significance level (1 - confidence level)
     * @return Pair of (CI half-width, relative precision)
     */
    private fun computeCIForMetric(batchMeans: List<Double>, alpha: Double): Pair<Double, Double> {
        if (batchMeans.size < 2) return Pair(0.0, 0.0)

        val n = batchMeans.size
        var sum = 0.0
        for (x in batchMeans) sum += x
        val mean = sum / n

        if (kotlin.math.abs(mean) < 1e-12) return Pair(0.0, 0.0)

        var variance = 0.0
        for (x in batchMeans) {
            val diff = x - mean
            variance += diff * diff
        }
        variance /= (n - 1)

        val stdErr = kotlin.math.sqrt(variance / n)
        val tCrit = getTCriticalValue(alpha, n - 1)
        val ciHalfWidth = tCrit * stdErr
        val relPrec = ciHalfWidth / kotlin.math.abs(mean)

        return Pair(ciHalfWidth, relPrec)
    }

    /** Cumulative queue time at last MSER sample */
    private lateinit var lastMserQueueTime: Array<DoubleArray>
    /** Time of last MSER sample */
    private var lastMserSampleTime = 0.0

    /**
     * MSER-5 sampling event - records time-weighted average queue length over each interval.
     */
    private inner class MSERSampleEvent : Event() {
        override fun actions() {
            val currentTime = Sim.time()
            val intervalDuration = currentTime - lastMserSampleTime

            // Update queue stats before sampling
            for (qIdx in 0 until numServiceNodes) {
                for (k in 0 until numClasses) {
                    updateQueueStats(qIdx, k)
                }
            }

            observationTimes.add(currentTime)

            // Record time-weighted average queue length over this interval
            for (qIdx in 0 until numServiceNodes) {
                for (k in 0 until numClasses) {
                    if (intervalDuration > 0) {
                        val intervalQueueTime = totalQueueTime[qIdx][k] - lastMserQueueTime[qIdx][k]
                        val avgQueueLength = intervalQueueTime / intervalDuration
                        queueLengthObservations[qIdx][k].add(avgQueueLength)
                        lastMserQueueTime[qIdx][k] = totalQueueTime[qIdx][k]
                    } else {
                        // Include blocked jobs (BAS, BBS, and FCR) trying to enter this queue
                        val effectiveQueueLength = currentQueueLength[qIdx][k] + basBlockedAtDest[qIdx][k] + bbsBlockedAtDest[qIdx][k] + fcrBlockedAtDest[qIdx][k]
                        queueLengthObservations[qIdx][k].add(effectiveQueueLength.toDouble())
                    }
                    throughputObservations[qIdx][k].add(completedCustomers[qIdx][k])
                }
            }
            lastMserSampleTime = currentTime

            // Schedule next sample
            if (currentTime + mserSamplingInterval < Double.MAX_VALUE) {
                MSERSampleEvent().schedule(mserSamplingInterval)
            }
        }
    }

    /**
     * Compute MSER-5 truncation point for a series of observations.
     * Returns the optimal batch index d* that minimizes MSER(d) = Var(Z[d:n]) / (n-d)^2
     * where Z are batch means of size 5.
     */
    private fun computeMSER5TruncationPoint(observations: List<Double>): Int {
        val n = observations.size
        if (n < effectiveMserBatchSize * 4) {
            // Not enough data for MSER-5, no truncation
            return 0
        }

        // Form batch means
        val numBatches = n / effectiveMserBatchSize
        val batchMeans = DoubleArray(numBatches)
        for (j in 0 until numBatches) {
            var sum = 0.0
            for (i in 0 until effectiveMserBatchSize) {
                sum += observations[j * effectiveMserBatchSize + i]
            }
            batchMeans[j] = sum / effectiveMserBatchSize
        }

        // Find truncation point that minimizes MSER
        var minMSER = Double.MAX_VALUE
        var optimalD = 0
        val maxD = numBatches / 2  // Search up to half the batches

        for (d in 0 until maxD) {
            val remainingBatches = numBatches - d
            if (remainingBatches < 2) break

            // Compute mean and variance of batches from d to end
            var sum = 0.0
            for (j in d until numBatches) {
                sum += batchMeans[j]
            }
            val mean = sum / remainingBatches

            var variance = 0.0
            for (j in d until numBatches) {
                val diff = batchMeans[j] - mean
                variance += diff * diff
            }
            variance /= (remainingBatches - 1)  // Sample variance

            // MSER(d) = variance / (n - d)^2
            val mser = variance / (remainingBatches.toDouble() * remainingBatches)

            if (mser < minMSER) {
                minMSER = mser
                optimalD = d
            }
        }

        return optimalD
    }

    /**
     * Apply MSER-5 truncation to compute statistics.
     * Called at end of simulation to determine truncation point and recompute statistics.
     */
    private fun applyMSER5Truncation() {
        if (!mserEnabled || queueLengthObservations.isEmpty()) return

        // Find truncation point using aggregate queue length across all queues/classes
        val aggregateObs = mutableListOf<Double>()
        for (i in 0 until observationTimes.size) {
            var total = 0.0
            for (qIdx in 0 until numServiceNodes) {
                for (k in 0 until numClasses) {
                    if (i < queueLengthObservations[qIdx][k].size) {
                        total += queueLengthObservations[qIdx][k][i]
                    }
                }
            }
            aggregateObs.add(total)
        }

        mserTruncationBatch = computeMSER5TruncationPoint(aggregateObs)
        val truncationObsIdx = mserTruncationBatch * effectiveMserBatchSize

        if (truncationObsIdx > 0 && truncationObsIdx < observationTimes.size) {
            warmupEndTime = observationTimes[truncationObsIdx]
            if (options.verbose == VerboseLevel.DEBUG) {
                System.out.println("MSER-5: Truncation at batch $mserTruncationBatch (t=${warmupEndTime})")
            }
        } else {
            warmupEndTime = 0.0
        }
    }

    private fun initializeGenerators() {
        routingRng = MRG32k3a()
        if (seed > 0) {
            routingRng.setSeed(longArrayOf(seed, seed + 1, seed + 2, seed + 3, seed + 4, seed + 5))
            siroRng = Random(seed + 99999)
        } else {
            siroRng = Random()
        }

        // Determine arrival process types from sn data structure for each source and class
        arrivalProcessType = Array(numSources) { srcIdx ->
            val istStation = sourceStations[srcIdx]
            Array(numClasses) { k ->
                val station = sn.stations[istStation]
                val jobClass = sn.jobclasses[k]
                sn.procid[station]?.get(jobClass) ?: ProcessType.DISABLED
            }
        }

        // Initialize PH process matrices for arrivals from sn.proc
        arrivalProc = Array(numSources) { srcIdx ->
            val istStation = sourceStations[srcIdx]
            Array(numClasses) { k ->
                val station = sn.stations[istStation]
                val jobClass = sn.jobclasses[k]
                sn.proc[station]?.get(jobClass)
            }
        }

        // Initialize random generators for PH/MAP/BMAP/MMAP arrival sampling
        arrivalRng = Array(numSources) { srcIdx ->
            Array(numClasses) { k ->
                val procType = arrivalProcessType[srcIdx][k]
                if (procType == ProcessType.PH || procType == ProcessType.APH ||
                    procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN ||
                    procType == ProcessType.COX2 || procType == ProcessType.MAP ||
                    procType == ProcessType.MMPP2 || procType == ProcessType.BMAP ||
                    procType == ProcessType.MMAP || procType == ProcessType.ME || procType == ProcessType.RAP) {
                    if (seed > 0) {
                        val offset = (srcIdx * numClasses + k).toLong() * 10 + 2000
                        Random(seed + offset)
                    } else {
                        Random()
                    }
                } else null
            }
        }

        // Initialize batch size cache for BMAP arrivals (default batch size = 1)
        arrivalBatchSize = Array(numSources) { IntArray(numClasses) { 1 } }

        // Arrival generators for each source and class (Generic RandomVariateGen)
        arrivalGens = Array(numSources) { srcIdx ->
            Array(numClasses) { k ->
                val procType = arrivalProcessType[srcIdx][k]
                val proc = arrivalProc[srcIdx][k]

                var gen: RandomVariateGen? = null
                val stream = MRG32k3a()
                if (seed > 0) {
                    val offset = (srcIdx * numClasses + k).toLong() * 10
                    stream.setSeed(longArrayOf(seed + offset, seed + offset + 1, seed + offset + 2,
                        seed + offset + 3, seed + offset + 4, seed + offset + 5))
                }

                if (procType == ProcessType.DET) {
                    // Det distribution: use arrival rate to compute mean inter-arrival time (1/rate)
                    if (lambdas[srcIdx][k] > 0) {
                        val mean = 1.0 / lambdas[srcIdx][k]
                        gen = ConstantGen(stream, mean)
                    } else {
                        throw RuntimeException("DES: Deterministic arrival distribution for source $srcIdx, class $k has invalid rate ${lambdas[srcIdx][k]}")
                    }
                } else if (procType == ProcessType.EXP && lambdas[srcIdx][k] > 0) {
                    gen = ExponentialGen(stream, lambdas[srcIdx][k])
                } else {
                    // For non-Markovian distributions, read parameters directly from sn.proc
                    // sn.proc now contains actual distribution parameters: {[param1], [param2], ...}
                    when (procType) {
                        ProcessType.UNIFORM -> {
                            // Uniform: {[min], [max]}
                            if (proc != null && proc.size() >= 2) {
                                val minVal = proc.get(0).get(0, 0)
                                val maxVal = proc.get(1).get(0, 0)
                                gen = UniformGen(stream, minVal, maxVal)
                            }
                        }
                        ProcessType.GAMMA -> {
                            // Gamma: {[shape], [scale]}
                            if (proc != null && proc.size() >= 2) {
                                val shape = proc.get(0).get(0, 0)
                                val scale = proc.get(1).get(0, 0)
                                gen = GammaGen(stream, shape, 1.0/scale) // SSJ GammaGen takes alpha (shape) and lambda (inverse scale)
                            }
                        }
                        // TODO: Weibull support disabled - parameter mapping needs investigation
                        // ProcessType.WEIBULL -> {
                        //     // Weibull: {[shape(r)], [scale(alpha)]}
                        //     if (proc != null && proc.size() >= 2) {
                        //         val r = proc.get(0).get(0, 0)      // shape
                        //         val alpha = proc.get(1).get(0, 0)  // scale
                        //         gen = WeibullGen(stream, r, 1.0/alpha, 0.0)
                        //     }
                        // }
                        ProcessType.LOGNORMAL -> {
                            // Lognormal: {[mu], [sigma]}
                            if (proc != null && proc.size() >= 2) {
                                val mu = proc.get(0).get(0, 0)
                                val sigma = proc.get(1).get(0, 0)
                                gen = LognormalGen(stream, mu, sigma)
                            }
                        }
                        ProcessType.PARETO -> {
                            // Pareto: {[shape], [scale]}
                            if (proc != null && proc.size() >= 2) {
                                val shape = proc.get(0).get(0, 0)
                                val scale = proc.get(1).get(0, 0)
                                gen = ParetoGen(stream, shape, scale)
                            }
                        }
                        ProcessType.ERLANG -> {
                            // Erlang stores D0 in index 0 (as PH) if it came from PH representation
                            // Check if it is PH-like (multi-phase or negative diagonal)
                            if (proc != null && proc.size() >= 2 &&
                                (proc.get(0).getNumRows() > 1 || proc.get(0).get(0,0) < 0)) {
                                // It's PH representation
                                val D0 = proc.get(0)
                                val phases = D0.getNumRows()
                                val lambda = -D0.get(0, 0)
                                if (phases > 0 && lambda > 0.0) {
                                    gen = ErlangGen(stream, phases, lambda)
                                }
                            } else if (proc != null && proc.size() >= 2) {
                                // Mean/SCV representation (fallback for backwards compatibility)
                                val mean = proc.get(0).get(0, 0)
                                val scv = proc.get(1).get(0, 0)
                                if (mean > 0 && scv > 0) {
                                    val shape = (1.0 / scv).toInt().coerceAtLeast(1)
                                    val scale = mean / shape
                                    gen = ErlangGen(stream, shape, 1.0/scale)
                                }
                            }
                        }
                        ProcessType.POISSON -> {
                            if (proc != null && proc.size() >= 1) {
                                val mean = proc.get(0).get(0, 0)
                                gen = PoissonGen(stream, mean)
                            }
                        }
                        ProcessType.BERNOULLI -> {
                            if (proc != null && proc.size() >= 1) {
                                val mean = proc.get(0).get(0, 0)
                                gen = BernoulliGen(stream, mean)
                            }
                        }
                        ProcessType.BINOMIAL -> {
                            if (proc != null && proc.size() >= 2) {
                                val mean = proc.get(0).get(0, 0)
                                val scv = proc.get(1).get(0, 0)
                                val p = 1.0 - scv * mean
                                if (p in 0.0..1.0) {
                                    val n = (mean / p).toInt()
                                    gen = BinomialGen(stream, n, p)
                                }
                            }
                        }
                        ProcessType.GEOMETRIC -> {
                            if (proc != null && proc.size() >= 1) {
                                val mean = proc.get(0).get(0, 0)
                                val p = 1.0 / (mean + 1.0)
                                gen = GeometricGen(stream, p)
                            }
                        }
                        ProcessType.IMMEDIATE -> {
                            // Immediate: constant near-zero time (1/Immediate rate)
                            val immTime = 1.0 / GlobalConstants.Immediate
                            gen = ConstantGen(immTime)
                        }
                        else -> {}
                    }
                }
                // Handle IMMEDIATE even when proc is null
                if (gen == null && procType == ProcessType.IMMEDIATE) {
                    gen = ConstantGen(1.0 / GlobalConstants.Immediate)
                }
                gen
            }
        }

        // Determine process types from sn data structure for each service node and class
        serviceProcessType = Array(numServiceNodes) { svcIdx ->
            val istStation = serviceStations[svcIdx]
            Array(numClasses) { k ->
                val station = sn.stations[istStation]
                val jobClass = sn.jobclasses[k]
                sn.procid[station]?.get(jobClass) ?: ProcessType.DISABLED
            }
        }

        // Initialize PH process matrices from sn.proc
        serviceProc = Array(numServiceNodes) { svcIdx ->
            val istStation = serviceStations[svcIdx]
            Array(numClasses) { k ->
                val station = sn.stations[istStation]
                val jobClass = sn.jobclasses[k]
                sn.proc[station]?.get(jobClass)
            }
        }

        // Initialize random generators for PH/MAP/MMAP sampling
        serviceRng = Array(numServiceNodes) { svcIdx ->
            Array(numClasses) { k ->
                val procType = serviceProcessType[svcIdx][k]
                if (procType == ProcessType.PH || procType == ProcessType.APH ||
                    procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN ||
                    procType == ProcessType.COX2 || procType == ProcessType.MAP || procType == ProcessType.MMPP2 ||
                    procType == ProcessType.MMAP || procType == ProcessType.ME || procType == ProcessType.RAP) {
                    if (seed > 0) {
                        val offset = ((numSources + svcIdx) * numClasses + k).toLong() * 10 + 1000
                        Random(seed + offset)
                    } else {
                        Random()
                    }
                } else null
            }
        }

        // Service generators for each service node and class (Generic RandomVariateGen)
        serviceGens = Array(numServiceNodes) { svcIdx ->
            Array(numClasses) { k ->
                val procType = serviceProcessType[svcIdx][k]
                val proc = serviceProc[svcIdx][k]

                var gen: RandomVariateGen? = null
                val stream = MRG32k3a()
                if (seed > 0) {
                    val offset = ((numSources + svcIdx) * numClasses + k).toLong() * 10 + 1000
                    stream.setSeed(longArrayOf(seed + offset, seed + offset + 1, seed + offset + 2,
                        seed + offset + 3, seed + offset + 4, seed + offset + 5))
                }

                if (procType == ProcessType.DET) {
                    // Det distribution: use service rate to compute mean (1/rate)
                    if (mus[svcIdx][k] > 0 && mus[svcIdx][k] < Double.MAX_VALUE) {
                        val mean = 1.0 / mus[svcIdx][k]
                        gen = ConstantGen(stream, mean)
                    } else {
                        throw RuntimeException("DES: Deterministic service distribution for station $svcIdx, class $k has invalid rate ${mus[svcIdx][k]}")
                    }
                } else if (procType == ProcessType.EXP && mus[svcIdx][k] > 0 && mus[svcIdx][k] < Double.MAX_VALUE) {
                    gen = ExponentialGen(stream, mus[svcIdx][k])
                } else {
                    // For non-Markovian distributions, read parameters directly from sn.proc
                    // sn.proc now contains actual distribution parameters: {[param1], [param2], ...}
                    when (procType) {
                        ProcessType.UNIFORM -> {
                            // Uniform: {[min], [max]}
                            if (proc != null && proc.size() >= 2) {
                                val minVal = proc.get(0).get(0, 0)
                                val maxVal = proc.get(1).get(0, 0)
                                gen = UniformGen(stream, minVal, maxVal)
                            }
                        }
                        ProcessType.GAMMA -> {
                            // Gamma: {[shape], [scale]}
                            if (proc != null && proc.size() >= 2) {
                                val shape = proc.get(0).get(0, 0)
                                val scale = proc.get(1).get(0, 0)
                                gen = GammaGen(stream, shape, 1.0/scale) // SSJ GammaGen takes alpha (shape) and lambda (inverse scale)
                            }
                        }
                        // TODO: Weibull support disabled - parameter mapping needs investigation
                        // ProcessType.WEIBULL -> {
                        //     // Weibull: {[shape(r)], [scale(alpha)]}
                        //     if (proc != null && proc.size() >= 2) {
                        //         val r = proc.get(0).get(0, 0)      // shape
                        //         val alpha = proc.get(1).get(0, 0)  // scale
                        //         gen = WeibullGen(stream, r, 1.0/alpha, 0.0)
                        //     }
                        // }
                        ProcessType.LOGNORMAL -> {
                            // Lognormal: {[mu], [sigma]}
                            if (proc != null && proc.size() >= 2) {
                                val mu = proc.get(0).get(0, 0)
                                val sigma = proc.get(1).get(0, 0)
                                gen = LognormalGen(stream, mu, sigma)
                            }
                        }
                        ProcessType.PARETO -> {
                            // Pareto: {[shape], [scale]}
                            if (proc != null && proc.size() >= 2) {
                                val shape = proc.get(0).get(0, 0)
                                val scale = proc.get(1).get(0, 0)
                                gen = ParetoGen(stream, shape, scale)
                            }
                        }
                        ProcessType.ERLANG -> {
                            // Erlang stores D0 in index 0 (as PH) if it came from PH representation
                            // Check if it is PH-like (multi-phase or negative diagonal)
                            if (proc != null && proc.size() >= 2 &&
                                (proc.get(0).getNumRows() > 1 || proc.get(0).get(0,0) < 0)) {
                                // It's PH representation
                                val D0 = proc.get(0)
                                val phases = D0.getNumRows()
                                val lambda = -D0.get(0, 0)
                                if (phases > 0 && lambda > 0.0) {
                                    gen = ErlangGen(stream, phases, lambda)
                                }
                            } else if (proc != null && proc.size() >= 2) {
                                // Mean/SCV representation (fallback for backwards compatibility)
                                val mean = proc.get(0).get(0, 0)
                                val scv = proc.get(1).get(0, 0)
                                if (mean > 0 && scv > 0) {
                                    val shape = (1.0 / scv).toInt().coerceAtLeast(1)
                                    val scale = mean / shape
                                    gen = ErlangGen(stream, shape, 1.0/scale)
                                }
                            }
                        }
                        ProcessType.POISSON -> {
                            if (proc != null && proc.size() >= 1) {
                                val mean = proc.get(0).get(0, 0)
                                gen = PoissonGen(stream, mean)
                            }
                        }
                        ProcessType.BERNOULLI -> {
                            if (proc != null && proc.size() >= 1) {
                                val mean = proc.get(0).get(0, 0)
                                gen = BernoulliGen(stream, mean)
                            }
                        }
                        ProcessType.BINOMIAL -> {
                            if (proc != null && proc.size() >= 2) {
                                val mean = proc.get(0).get(0, 0)
                                val scv = proc.get(1).get(0, 0)
                                val p = 1.0 - scv * mean
                                if (p in 0.0..1.0) {
                                    val n = (mean / p).toInt()
                                    gen = BinomialGen(stream, n, p)
                                }
                            }
                        }
                        ProcessType.GEOMETRIC -> {
                            if (proc != null && proc.size() >= 1) {
                                val mean = proc.get(0).get(0, 0)
                                val p = 1.0 / (mean + 1.0)
                                gen = GeometricGen(stream, p)
                            }
                        }
                        ProcessType.IMMEDIATE -> {
                            // Immediate: constant near-zero time (1/Immediate rate)
                            val immTime = 1.0 / GlobalConstants.Immediate
                            gen = ConstantGen(immTime)
                        }
                        else -> {}
                    }
                }
                // Handle IMMEDIATE even when proc is null
                if (gen == null && procType == ProcessType.IMMEDIATE) {
                    gen = ConstantGen(1.0 / GlobalConstants.Immediate)
                }
                gen
            }
        }

        // Initialize trace samplers for Replayer distributions (arrivals)
        arrivalTraceSamplers = Array(numSources) { srcIdx ->
            val istStation = sourceStations[srcIdx]
            val station = sn.stations[istStation]
            Array(numClasses) { k ->
                val procType = arrivalProcessType[srcIdx][k]
                if (procType == ProcessType.REPLAYER && station is Source) {
                    val jobClass = sn.jobclasses[k]
                    val distr = station.getArrivalDistribution(jobClass)
                    if (distr is Replayer) {
                        // Access the data field using reflection since it's package-private
                        val dataField = Replayer::class.java.getDeclaredField("data")
                        dataField.isAccessible = true
                        val data = dataField.get(distr) as DoubleArray
                        TraceSampler(data)
                    } else null
                } else null
            }
        }

        // Initialize trace samplers for Replayer distributions (service)
        serviceTraceSamplers = Array(numServiceNodes) { svcIdx ->
            val istStation = serviceStations[svcIdx]
            val station = sn.stations[istStation]
            Array(numClasses) { k ->
                val procType = serviceProcessType[svcIdx][k]
                if (procType == ProcessType.REPLAYER) {
                    val jobClass = sn.jobclasses[k]
                    val distr = station.server?.getServiceDistribution(jobClass)
                    if (distr is Replayer) {
                        // Access the data field using reflection since it's package-private
                        val dataField = Replayer::class.java.getDeclaredField("data")
                        dataField.isAccessible = true
                        val data = dataField.get(distr) as DoubleArray
                        TraceSampler(data)
                    } else null
                } else null
            }
        }

        // Initialize setup and delayoff generators for queues with setup/delayoff enabled
        // First, initialize the hasSetupDelayoff array
        hasSetupDelayoff = BooleanArray(numServiceNodes) { false }

        setupGens = Array(numServiceNodes) { svcIdx ->
            val istStation = serviceStations[svcIdx]
            val station = sn.stations[istStation]

            Array(numClasses) { k ->
                if (station is Queue && station.isDelayOffEnabled()) {
                    hasSetupDelayoff[svcIdx] = true
                    val jobClass = sn.jobclasses[k]
                    val setupDist = station.getSetupTime(jobClass)

                    if (setupDist != null && setupDist.getMean() > 0) {
                        // Create exponential generator with rate = 1/mean
                        // TODO: Support other distribution types (Erlang, HyperExp, PH, etc.)
                        val stream = MRG32k3a()
                        if (seed > 0) {
                            val offset = ((numSources + numServiceNodes + svcIdx) * numClasses + k).toLong() * 10 + 3000
                            stream.setSeed(longArrayOf(seed + offset, seed + offset + 1, seed + offset + 2,
                                seed + offset + 3, seed + offset + 4, seed + offset + 5))
                        }
                        ExponentialGen(stream, 1.0 / setupDist.getMean())
                    } else null
                } else null
            }
        }

        delayoffGens = Array(numServiceNodes) { svcIdx ->
            val istStation = serviceStations[svcIdx]
            val station = sn.stations[istStation]

            Array(numClasses) { k ->
                if (station is Queue && station.isDelayOffEnabled()) {
                    val jobClass = sn.jobclasses[k]
                    val delayoffDist = station.getDelayOffTime(jobClass)

                    if (delayoffDist != null && delayoffDist.getMean() > 0) {
                        // Create exponential generator with rate = 1/mean
                        // TODO: Support other distribution types (Erlang, HyperExp, PH, etc.)
                        val stream = MRG32k3a()
                        if (seed > 0) {
                            val offset = ((numSources + numServiceNodes + svcIdx) * numClasses + k).toLong() * 10 + 4000
                            stream.setSeed(longArrayOf(seed + offset, seed + offset + 1, seed + offset + 2,
                                seed + offset + 3, seed + offset + 4, seed + offset + 5))
                        }
                        ExponentialGen(stream, 1.0 / delayoffDist.getMean())
                    } else null
                } else null
            }
        }

        // Initialize antithetic variates generators if enabled
        initializeAntitheticGenerators()

        // Initialize control variates tracking if enabled
        initializeControlVariates()
    }

    /**
     * Initialize antithetic variates generators.
     * Creates a second set of generators using AntitheticStream wrappers for synchronized 1-U sampling.
     */
    private fun initializeAntitheticGenerators() {
        if (!useAntitheticVariates) {
            // Initialize empty arrays to avoid lateinit issues
            antitheticArrivalGens = Array(numSources) { Array(numClasses) { null as RandomVariateGen? } }
            antitheticServiceGens = Array(numServiceNodes) { Array(numClasses) { null as RandomVariateGen? } }
            antitheticArrivalRng = Array(numSources) { Array(numClasses) { null as Random? } }
            antitheticServiceRng = Array(numServiceNodes) { Array(numClasses) { null as Random? } }
            return
        }

        // Offset for antithetic seed to ensure different but deterministic streams
        val antitheticSeedOffset = 50000L

        // Initialize antithetic arrival generators
        antitheticArrivalGens = Array(numSources) { srcIdx ->
            Array(numClasses) { k ->
                val procType = arrivalProcessType[srcIdx][k]

                var gen: RandomVariateGen? = null
                val stream = MRG32k3a()
                if (seed > 0) {
                    val offset = (srcIdx * numClasses + k).toLong() * 10 + antitheticSeedOffset
                    stream.setSeed(longArrayOf(seed + offset, seed + offset + 1, seed + offset + 2,
                        seed + offset + 3, seed + offset + 4, seed + offset + 5))
                }
                // Wrap the stream with AntitheticStream for 1-U transformation
                val antiStream = umontreal.ssj.rng.AntitheticStream(stream)

                if (procType == ProcessType.EXP && lambdas[srcIdx][k] > 0) {
                    gen = ExponentialGen(antiStream, lambdas[srcIdx][k])
                }
                gen
            }
        }

        // Initialize antithetic service generators
        antitheticServiceGens = Array(numServiceNodes) { svcIdx ->
            Array(numClasses) { k ->
                val procType = serviceProcessType[svcIdx][k]

                var gen: RandomVariateGen? = null
                val stream = MRG32k3a()
                if (seed > 0) {
                    val offset = (numSources * numClasses + svcIdx * numClasses + k).toLong() * 10 + antitheticSeedOffset
                    stream.setSeed(longArrayOf(seed + offset, seed + offset + 1, seed + offset + 2,
                        seed + offset + 3, seed + offset + 4, seed + offset + 5))
                }
                val antiStream = umontreal.ssj.rng.AntitheticStream(stream)

                if (procType == ProcessType.EXP && mus[svcIdx][k] > 0 && mus[svcIdx][k] < Double.MAX_VALUE) {
                    gen = ExponentialGen(antiStream, mus[svcIdx][k])
                }
                gen
            }
        }

        // Initialize antithetic RNG for PH/MAP/MMAP distributions
        antitheticArrivalRng = Array(numSources) { srcIdx ->
            Array(numClasses) { k ->
                val procType = arrivalProcessType[srcIdx][k]
                if (procType == ProcessType.PH || procType == ProcessType.APH ||
                    procType == ProcessType.MAP || procType == ProcessType.MMPP2 ||
                    procType == ProcessType.MMAP) {
                    val rng = Random()
                    if (seed > 0) {
                        rng.setSeed(seed + antitheticSeedOffset + (srcIdx * numClasses + k) * 100)
                    }
                    rng
                } else null
            }
        }

        antitheticServiceRng = Array(numServiceNodes) { svcIdx ->
            Array(numClasses) { k ->
                val procType = serviceProcessType[svcIdx][k]
                if (procType == ProcessType.PH || procType == ProcessType.APH ||
                    procType == ProcessType.MAP || procType == ProcessType.MMPP2 ||
                    procType == ProcessType.MMAP) {
                    val rng = Random()
                    if (seed > 0) {
                        rng.setSeed(seed + antitheticSeedOffset + (numSources * numClasses + svcIdx * numClasses + k) * 100)
                    }
                    rng
                } else null
            }
        }
    }

    /**
     * Initialize control variates tracking structures.
     */
    private fun initializeControlVariates() {
        // Always initialize arrays to avoid lateinit issues
        arrivalSampleSum = Array(numSources) { DoubleArray(numClasses) { 0.0 } }
        arrivalSampleCount = Array(numSources) { LongArray(numClasses) { 0L } }
        serviceSampleSum = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        serviceSampleCount = Array(numServiceNodes) { LongArray(numClasses) { 0L } }

        // Initialize expected means (1/rate)
        arrivalExpectedMean = Array(numSources) { srcIdx ->
            DoubleArray(numClasses) { k ->
                if (lambdas[srcIdx][k] > 0) 1.0 / lambdas[srcIdx][k] else 0.0
            }
        }

        serviceExpectedMean = Array(numServiceNodes) { svcIdx ->
            DoubleArray(numClasses) { k ->
                if (mus[svcIdx][k] > 0 && mus[svcIdx][k] < Double.MAX_VALUE)
                    1.0 / mus[svcIdx][k] else 0.0
            }
        }

        // Initialize heterogeneous server support
        initializeHeterogeneousServers()
    }

    /**
     * Initialize heterogeneous server support.
     * Reads server type configuration from NetworkStruct and sets up
     * per-type tracking arrays and service generators.
     */
    private fun initializeHeterogeneousServers() {
        // Initialize arrays with default values (homogeneous case)
        numServerTypes = IntArray(numServiceNodes) { 0 }
        serversPerType = Array(numServiceNodes) { IntArray(0) }
        serverCompat = Array(numServiceNodes) { Array(0) { BooleanArray(0) } }
        busyCountPerType = Array(numServiceNodes) { IntArray(0) }
        serverToType = Array(numServiceNodes) { IntArray(0) }
        heteroSchedPolicy = Array(numServiceNodes) { null }
        heteroServiceGens = Array(numServiceNodes) { Array(0) { Array(0) { null as RandomVariateGen? } } }
        heteroMus = Array(numServiceNodes) { Array(0) { DoubleArray(0) } }
        heteroServiceProcType = Array(numServiceNodes) { Array(0) { Array(0) { ProcessType.EXP } } }
        heteroServiceProc = Array(numServiceNodes) { Array(0) { Array(0) { null as MatrixCell? } } }
        heteroServiceRng = Array(numServiceNodes) { Array(0) { Array(0) { null as Random? } } }
        serverTypeOrder = Array(numServiceNodes) { mutableListOf() }
        alfsOrder = Array(numServiceNodes) { IntArray(0) }

        // Check if any station has heterogeneous servers
        if (sn.nservertypes == null) return

        for (svcIdx in 0 until numServiceNodes) {
            if (isDelayNode[svcIdx]) continue // Skip delay nodes

            val stationIdx = serviceStations[svcIdx]
            if (stationIdx >= sn.nservertypes.length()) continue

            val nTypes = sn.nservertypes.get(stationIdx).toInt()
            if (nTypes <= 0) continue

            numServerTypes[svcIdx] = nTypes
            val station = sn.stations[stationIdx]

            // Get servers per type
            val serversMatrix = sn.serverspertype?.get(station)
            serversPerType[svcIdx] = IntArray(nTypes) { typeId ->
                serversMatrix?.get(typeId)?.toInt() ?: 1
            }

            // Get compatibility matrix: servercompat[station] is (nTypes x K)
            val compatMatrix = sn.servercompat?.get(station)
            serverCompat[svcIdx] = Array(nTypes) { typeId ->
                BooleanArray(numClasses) { classId ->
                    val compatible = compatMatrix?.get(typeId, classId) ?: 0.0
                    compatible > 0.5
                }
            }

            // Initialize busy counts per type
            busyCountPerType[svcIdx] = IntArray(nTypes) { 0 }

            // Build server ID to type mapping
            val totalServers = serversPerType[svcIdx].sum()
            serverToType[svcIdx] = IntArray(totalServers)
            var globalServerId = 0
            for (typeId in 0 until nTypes) {
                repeat(serversPerType[svcIdx][typeId]) {
                    if (globalServerId < totalServers) {
                        serverToType[svcIdx][globalServerId] = typeId
                        globalServerId++
                    }
                }
            }

            // Get heterogeneous scheduling policy
            heteroSchedPolicy[svcIdx] = sn.heteroschedpolicy?.get(station) ?: HeteroSchedPolicy.ORDER

            // Initialize server type order for ALIS/FAIRNESS (round-robin)
            serverTypeOrder[svcIdx] = (0 until nTypes).toMutableList()

            // Initialize ALFS order (sorted by number of compatible classes, ascending)
            alfsOrder[svcIdx] = (0 until nTypes).sortedBy { typeId ->
                serverCompat[svcIdx][typeId].count { it }
            }.toIntArray()

            // Get heterogeneous service rates
            val stationRates = sn.heterorates?.get(station)
            heteroMus[svcIdx] = Array(nTypes) { typeId ->
                DoubleArray(numClasses) { classId ->
                    stationRates?.get(typeId)?.get(classId) ?: Double.MAX_VALUE
                }
            }

            // Get heterogeneous process types
            val stationProcId = sn.heteroprocid?.get(station)
            heteroServiceProcType[svcIdx] = Array(nTypes) { typeId ->
                Array(numClasses) { classId ->
                    stationProcId?.get(typeId)?.get(classId) ?: ProcessType.EXP
                }
            }

            // Get heterogeneous PH processes
            val stationProc = sn.heteroproc?.get(station)
            heteroServiceProc[svcIdx] = Array(nTypes) { typeId ->
                Array(numClasses) { classId ->
                    stationProc?.get(typeId)?.get(classId)
                }
            }

            // Initialize RNG for PH sampling
            heteroServiceRng[svcIdx] = Array(nTypes) { typeId ->
                Array(numClasses) { classId ->
                    if (seed > 0) {
                        val offset = (svcIdx * 1000L + typeId * 100L + classId) + 8000
                        Random(seed + offset)
                    } else {
                        Random()
                    }
                }
            }

            // Create heterogeneous service generators
            heteroServiceGens[svcIdx] = Array(nTypes) { typeId ->
                Array(numClasses) { classId ->
                    createHeteroServiceGenerator(svcIdx, typeId, classId)
                }
            }
        }
    }

    /**
     * Creates a service time generator for a specific server type and job class.
     */
    private fun createHeteroServiceGenerator(svcIdx: Int, typeId: Int, classId: Int): RandomVariateGen? {
        val rate = heteroMus[svcIdx][typeId][classId]
        if (rate <= 0 || rate >= Double.MAX_VALUE) return null

        val procType = heteroServiceProcType[svcIdx][typeId][classId]
        val proc = heteroServiceProc[svcIdx][typeId][classId]

        val stream = MRG32k3a()
        if (seed > 0) {
            val offset = (svcIdx * 1000L + typeId * 100L + classId) + 7000
            stream.setSeed(longArrayOf(seed + offset, seed + offset + 1, seed + offset + 2,
                seed + offset + 3, seed + offset + 4, seed + offset + 5))
        }

        return when (procType) {
            ProcessType.DET -> ConstantGen(stream, 1.0 / rate)
            ProcessType.EXP -> ExponentialGen(stream, rate)
            ProcessType.ERLANG -> {
                // Erlang: proc = {[D0], [D1]} where D0 has shape phases
                if (proc != null && proc.size() >= 2) {
                    val D0 = proc.get(0)
                    val nPhases = D0.getNumCols()
                    val phaseRate = rate * nPhases
                    ErlangGen(stream, nPhases, phaseRate)
                } else null
            }
            ProcessType.HYPEREXP -> {
                // HyperExp will use map_sample at runtime
                null
            }
            ProcessType.PH, ProcessType.APH, ProcessType.COXIAN, ProcessType.COX2, ProcessType.MAP -> {
                // PH distributions will use map_sample at runtime
                null
            }
            else -> {
                // Try exponential as fallback
                if (rate > 0 && rate < Double.MAX_VALUE) {
                    ExponentialGen(stream, rate)
                } else null
            }
        }
    }

    private fun initializeState() {
        // Create wait queues with appropriate comparator based on scheduling strategy
        waitQueues = Array(numServiceNodes) { svcIdx ->
            PriorityQueue<Customer>(getComparatorForStrategy(schedStrategies[svcIdx], svcIdx))
        }
        // For Delay nodes (infinite servers), use empty array since we don't track individual servers
        serverBusy = Array(numServiceNodes) { svcIdx ->
            if (isDelayNode[svcIdx]) BooleanArray(0) else BooleanArray(numServers[svcIdx]) { false }
        }
        // Server blocked state for synchronous calls (waiting for REPLY signal)
        serverBlocked = Array(numServiceNodes) { svcIdx ->
            if (isDelayNode[svcIdx]) BooleanArray(0) else BooleanArray(numServers[svcIdx]) { false }
        }
        customersInService = IntArray(numServiceNodes) { 0 }

        responseTimeTally = Array(numServiceNodes) { svcIdx ->
            Array(numClasses) { k -> Tally("Response time S$svcIdx C$k") }
        }
        responseTimeSamples = Array(numServiceNodes) {
            Array(numClasses) { mutableListOf<Double>() }
        }
        completedCustomers = Array(numServiceNodes) { IntArray(numClasses) { 0 } }
        totalQueueTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        lastQueueUpdateTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        currentQueueLength = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        // Initialize utilization tracking
        totalBusyTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        lastBusyUpdateTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        currentBusyServers = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        // Initialize blocking time tracking for REPLY signals
        totalBlockingTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        currentBlockedServers = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        // Initialize blocked-at-destination tracking (blocked jobs count at destination for queue length)
        basBlockedAtDest = Array(numServiceNodes) { IntArray(numClasses) { 0 } }
        bbsBlockedAtDest = Array(numServiceNodes) { IntArray(numClasses) { 0 } }
        fcrBlockedAtDest = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        // Initialize setup/delayoff state tracking (all servers start ACTIVE)
        serverState = Array(numServiceNodes) { svcIdx ->
            if (isDelayNode[svcIdx]) {
                emptyArray()  // Delay nodes don't track individual servers
            } else {
                Array(numServers[svcIdx]) { ServerState.ACTIVE }
            }
        }
        // Note: hasSetupDelayoff is initialized in initializeGenerators() before use
        serverLastClass = Array(numServiceNodes) { svcIdx ->
            if (isDelayNode[svcIdx]) IntArray(0) else IntArray(numServers[svcIdx]) { -1 }
        }
        pendingDelayoffEvents = Array(numServiceNodes) { svcIdx ->
            if (isDelayNode[svcIdx]) emptyArray() else Array<Event?>(numServers[svcIdx]) { null }
        }

        // Initialize setup time tracking
        totalSetupTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        lastSetupUpdateTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        currentServersInSetup = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        // Initialize delayoff time tracking
        totalDelayoffTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        lastDelayoffUpdateTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        currentServersInDelayoff = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        systemResponseTimeTally = Array(numClasses) { k -> Tally("System response time C$k") }
        systemTardinessTally = Array(numClasses) { k -> Tally("System tardiness C$k") }
        tardinessTally = Array(numServiceNodes) { svcIdx ->
            Array(numClasses) { k -> Tally("Tardiness node${svcIdx} C$k") }
        }
        systemCompletedCustomers = IntArray(numClasses) { 0 }

        // Initialize dropped customers tracking
        droppedCustomers = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        // Initialize arrived customers tracking (for arrival rate calculation)
        arrivedCustomers = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        // Initialize finite capacity region tracking
        currentJobsInRegion = Array(numRegions) { IntArray(numClasses) { 0 } }
        droppedByRegion = Array(numRegions) { IntArray(numClasses) { 0 } }

        // Initialize FCR time-weighted metrics tracking
        totalRegionJobTime = Array(numRegions) { DoubleArray(numClasses) { 0.0 } }
        lastRegionUpdateTime = DoubleArray(numRegions) { 0.0 }
        regionCompletions = Array(numRegions) { IntArray(numClasses) { 0 } }
        regionResponseTimeTally = Array(numRegions) { f ->
            Array(numClasses) { k -> Tally("FCR$f response time C$k") }
        }

        // Initialize FCR blocking (waiting queue) support
        fcRegionBlockedQueue = Array(numRegions) { java.util.LinkedList<BlockedCustomer>() }
        blockedInRegion = Array(numRegions) { IntArray(numClasses) { 0 } }

        // Initialize FCR class weights from NetworkStruct (JMT parity: each job class can have different weight for capacity)
        fcRegionClassWeights = Array(numRegions) { f ->
            DoubleArray(numClasses) { k ->
                if (sn.regionweight != null) sn.regionweight.get(f, k) else 1.0
            }
        }

        // Initialize FCR arrival rate tracking (for computing ANfcr using inter-arrival times)
        lastRegionArrivalTime = Array(numRegions) { DoubleArray(numClasses) { 0.0 } }
        regionArrivalCount = Array(numRegions) { IntArray(numClasses) { 0 } }
        regionInterArrivalTimeSum = Array(numRegions) { DoubleArray(numClasses) { 0.0 } }

        // Initialize impatience (reneging, balking, retrial) statistics and configuration
        initializeImpatienceSupport()

        // Initialize PS (Processor Sharing) scheduling state
        psJobsInService = Array(numServiceNodes) { ArrayList<PSCustomer>() }
        psLastUpdateTime = DoubleArray(numServiceNodes) { 0.0 }
        psLastBusyUpdateTime = DoubleArray(numServiceNodes) { 0.0 }

        // Initialize preemptive LCFS and SRPT scheduling state
        preemptiveJobsInService = Array(numServiceNodes) { ArrayList<PreemptiveCustomer>() }
        isPreemptiveScheduling = BooleanArray(numServiceNodes) { idx ->
            isPreemptiveLCFSScheduling(schedStrategies[idx]) || isSizeBasedPreemptiveScheduling(schedStrategies[idx])
        }
        preemptedJobHistory = mutableMapOf()

        // Cache DPS/GPS weights from sn.schedparam
        schedWeights = Array(numServiceNodes) { svcIdx ->
            DoubleArray(numClasses) { k ->
                val stationIdx = serviceStations[svcIdx]
                val weight = sn.schedparam.get(stationIdx, k)
                if (weight > 0 && !weight.isNaN()) weight else 1.0
            }
        }

        // Initialize load-dependent service support
        initializeLoadDependentService()

        // Initialize limited joint-dependent (LJD) scaling
        initializeLjdScaling()

        // Initialize limited joint class-dependent (LJCD) scaling
        initializeLjcdScaling()

        // Initialize polling scheduling support
        initializePollingState()
    }

    /**
     * Initialize impatience support (reneging, balking, retrial).
     * Sets up statistics arrays, configuration caches, and random generators for patience/retrial times.
     */
    private fun initializeImpatienceSupport() {
        // Initialize reneging statistics
        renegedCustomers = Array(numServiceNodes) { IntArray(numClasses) { 0 } }
        totalRenegingWaitTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }

        // Initialize balking statistics
        balkedCustomers = Array(numServiceNodes) { IntArray(numClasses) { 0 } }

        // Initialize retrial/orbit statistics
        orbitJobs = Array(numServiceNodes) { mutableListOf<OrbitJob>() }
        retriedCustomers = Array(numServiceNodes) { IntArray(numClasses) { 0 } }
        maxRetriesExceeded = Array(numServiceNodes) { IntArray(numClasses) { 0 } }
        currentOrbitSize = Array(numServiceNodes) { IntArray(numClasses) { 0 } }
        totalOrbitTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }
        lastOrbitUpdateTime = Array(numServiceNodes) { DoubleArray(numClasses) { 0.0 } }

        // Initialize configuration caches
        hasPatienceConfig = Array(numServiceNodes) { BooleanArray(numClasses) { false } }
        patienceGens = Array(numServiceNodes) { arrayOfNulls<RandomVariateGen>(numClasses) }
        hasBalkingConfig = Array(numServiceNodes) { BooleanArray(numClasses) { false } }
        hasRetrialConfig = Array(numServiceNodes) { BooleanArray(numClasses) { false } }
        retrialGens = Array(numServiceNodes) { arrayOfNulls<RandomVariateGen>(numClasses) }
        retrialMaxAttemptsConfig = Array(numServiceNodes) { IntArray(numClasses) { -1 } }

        // Read impatience configuration from NetworkStruct
        for (svcIdx in 0 until numServiceNodes) {
            val istStation = serviceStations[svcIdx]
            val station = sn.stations[istStation]

            for (k in 0 until numClasses) {
                val jobClass = sn.jobclasses[k]

                // Check for patience/reneging configuration
                val impatienceType = sn.impatienceType?.get(station)?.get(jobClass)
                if (impatienceType != null && impatienceType != ProcessType.DISABLED) {
                    hasPatienceConfig[svcIdx][k] = true

                    // Create patience time generator
                    val stream = MRG32k3a()
                    if (seed > 0) {
                        val offset = (svcIdx * numClasses + k).toLong() * 10 + 50000
                        stream.setSeed(longArrayOf(seed + offset, seed + offset + 1, seed + offset + 2,
                            seed + offset + 3, seed + offset + 4, seed + offset + 5))
                    }

                    patienceGens[svcIdx][k] = createImpatienceGenerator(
                        impatienceType,
                        sn.impatienceMu?.get(station)?.get(jobClass),
                        sn.impatiencePhi?.get(station)?.get(jobClass),
                        sn.impatienceProc?.get(station)?.get(jobClass),
                        stream,
                        svcIdx, k, "patience"
                    )
                }

                // Check for balking configuration
                val balkStrat = sn.balkingStrategy?.get(station)?.get(jobClass)
                if (balkStrat != null) {
                    hasBalkingConfig[svcIdx][k] = true
                }

                // Check for retrial configuration
                val retrialType = sn.retrialType?.get(station)?.get(jobClass)
                if (retrialType != null && retrialType != ProcessType.DISABLED) {
                    hasRetrialConfig[svcIdx][k] = true

                    // Create retrial delay generator
                    val stream = MRG32k3a()
                    if (seed > 0) {
                        val offset = (svcIdx * numClasses + k).toLong() * 10 + 60000
                        stream.setSeed(longArrayOf(seed + offset, seed + offset + 1, seed + offset + 2,
                            seed + offset + 3, seed + offset + 4, seed + offset + 5))
                    }

                    retrialGens[svcIdx][k] = createImpatienceGenerator(
                        retrialType,
                        sn.retrialMu?.get(station)?.get(jobClass),
                        sn.retrialPhi?.get(station)?.get(jobClass),
                        sn.retrialProc?.get(station)?.get(jobClass),
                        stream,
                        svcIdx, k, "retrial"
                    )

                    // Get max retrial attempts (-1 = unlimited)
                    retrialMaxAttemptsConfig[svcIdx][k] = sn.retrialMaxAttempts?.get(station)?.get(jobClass) ?: -1
                }
            }
        }
    }

    /**
     * Create a random variate generator for impatience (patience/retrial) distributions.
     * Supports EXP, ERLANG, HYPEREXP, PH, DET, UNIFORM, GAMMA, PARETO, WEIBULL, LOGNORMAL
     * and other standard distributions.
     */
    private fun createImpatienceGenerator(
        procType: ProcessType,
        muMatrix: jline.util.matrix.Matrix?,
        phiMatrix: jline.util.matrix.Matrix?,
        proc: jline.util.matrix.MatrixCell?,
        stream: MRG32k3a,
        svcIdx: Int,
        classIdx: Int,
        label: String
    ): RandomVariateGen? {
        if (procType == ProcessType.DISABLED) return null

        val rate = muMatrix?.get(0, 0) ?: 0.0
        val scv = phiMatrix?.get(0, 0) ?: 1.0  // Squared coefficient of variation (default 1.0 for exponential)

        return when (procType) {
            ProcessType.EXP -> {
                if (rate > 0 && rate < Double.MAX_VALUE) {
                    ExponentialGen(stream, rate)
                } else null
            }
            ProcessType.DET -> {
                if (rate > 0 && rate < Double.MAX_VALUE) {
                    ConstantGen(stream, 1.0 / rate)
                } else null
            }
            ProcessType.ERLANG -> {
                if (proc != null && proc.size() >= 2) {
                    val k = proc.get(0).get(0, 0).toInt()
                    val lambda = proc.get(1).get(0, 0)
                    if (k > 0 && lambda > 0) {
                        ErlangGen(stream, k, lambda)
                    } else null
                } else null
            }
            ProcessType.HYPEREXP -> {
                if (proc != null && proc.size() >= 2) {
                    val d0 = proc.get(0)
                    val d1 = proc.get(1)
                    if (d0 != null && d1 != null) {
                        // HyperExp: D0 has negative rates on diagonal, D1 has absorption rates
                        val lambda1 = -d0.get(0, 0)
                        val lambda2 = if (d0.getNumRows() > 1) -d0.get(1, 1) else lambda1
                        val p = if (lambda1 > 0) d1.get(0, 0) / lambda1 else 0.5
                        HyperExponentialDistGen(stream, doubleArrayOf(p, 1.0 - p), doubleArrayOf(lambda1, lambda2))
                    } else null
                } else null
            }
            ProcessType.UNIFORM -> {
                if (proc != null && proc.size() >= 2) {
                    val minVal = proc.get(0).get(0, 0)
                    val maxVal = proc.get(1).get(0, 0)
                    if (maxVal > minVal) {
                        UniformGen(stream, minVal, maxVal)
                    } else null
                } else null
            }
            ProcessType.GAMMA -> {
                // Gamma distribution: shape = 1/scv, scale = scv/rate
                // SSJ GammaGen takes alpha (shape) and lambda (1/scale = rate/scv)
                if (rate > 0 && rate < Double.MAX_VALUE && scv > 0) {
                    val shape = 1.0 / scv
                    val lambda = rate / scv  // 1/scale where scale = scv/rate
                    GammaGen(stream, shape, lambda)
                } else null
            }
            ProcessType.PARETO -> {
                // Pareto distribution: shape = sqrt(1 + 1/scv) + 1, scale computed from mean
                // Mean = scale * shape / (shape - 1) => scale = mean * (shape - 1) / shape
                if (rate > 0 && rate < Double.MAX_VALUE && scv > 0) {
                    val mean = 1.0 / rate
                    val shape = FastMath.sqrt(1.0 + 1.0 / scv) + 1.0
                    val scale = mean * (shape - 1.0) / shape
                    if (shape > 1.0 && scale > 0) {
                        ParetoGen(stream, shape, scale)
                    } else null
                } else null
            }
            ProcessType.WEIBULL -> {
                // Weibull distribution using Justus approximation (1976) for shape parameter
                // c = sqrt(scv), r (shape) = c^(-1.086), alpha (scale) = mean / Gamma(1 + 1/r)
                if (rate > 0 && rate < Double.MAX_VALUE && scv > 0) {
                    val mean = 1.0 / rate
                    val c = FastMath.sqrt(scv)
                    val r = FastMath.pow(c, -1.086)  // shape parameter
                    val alpha = mean / org.apache.commons.math3.special.Gamma.gamma(1.0 + 1.0 / r)  // scale
                    if (r > 0 && alpha > 0) {
                        // SSJ WeibullGen(stream, alpha, lambda, delta) where:
                        // alpha = shape, lambda = 1/scale, delta = location (0 for standard)
                        WeibullGen(stream, r, 1.0 / alpha, 0.0)
                    } else null
                } else null
            }
            ProcessType.LOGNORMAL -> {
                // Lognormal distribution: mu and sigma from mean and scv
                // c = sqrt(scv), mu = log(mean / sqrt(c^2 + 1)), sigma = sqrt(log(c^2 + 1))
                if (rate > 0 && rate < Double.MAX_VALUE && scv > 0) {
                    val mean = 1.0 / rate
                    val c = FastMath.sqrt(scv)
                    val c2plus1 = c * c + 1.0
                    val mu = FastMath.log(mean / FastMath.sqrt(c2plus1))
                    val sigma = FastMath.sqrt(FastMath.log(c2plus1))
                    if (sigma > 0) {
                        LognormalGen(stream, mu, sigma)
                    } else null
                } else null
            }
            ProcessType.PH, ProcessType.APH, ProcessType.COXIAN, ProcessType.COX2, ProcessType.MAP, ProcessType.MMPP2, ProcessType.ME, ProcessType.RAP -> {
                // Phase-type and MAP distributions: use map_sample directly via wrapper
                if (proc != null && proc.size() >= 2) {
                    val rng = if (seed > 0) {
                        val offset = (svcIdx * numClasses + classIdx).toLong() * 10 + 70000
                        Random(seed + offset)
                    } else Random()
                    MapSampleGen(stream, proc, rng)
                } else null
            }
            else -> {
                // Default: if we have a rate, use exponential
                if (rate > 0 && rate < Double.MAX_VALUE) {
                    ExponentialGen(stream, rate)
                } else null
            }
        }
    }

    /**
     * Initialize polling scheduling support.
     * Sets up per-class queues, switchover generators, and polling state for stations with POLLING strategy.
     */
    private fun initializePollingState() {
        isPollingStation = BooleanArray(numServiceNodes) { svcIdx ->
            schedStrategies[svcIdx] == SchedStrategy.POLLING
        }

        pollingType = Array(numServiceNodes) { null }
        pollingK = IntArray(numServiceNodes) { 1 }
        pollingCurrentClass = IntArray(numServiceNodes) { 0 }
        pollingJobsServedInRound = IntArray(numServiceNodes) { 0 }
        pollingGateSize = IntArray(numServiceNodes) { 0 }
        pollingInSwitchover = BooleanArray(numServiceNodes) { false }

        // Initialize per-class queues for polling stations
        pollingQueues = Array(numServiceNodes) { svcIdx ->
            Array(numClasses) { LinkedList<Customer>() }
        }

        // Initialize switchover generators
        pollingSwitchoverGens = Array(numServiceNodes) { svcIdx ->
            Array(numClasses) { null as RandomVariateGen? }
        }

        // Extract polling parameters from sn.nodeparam
        for (svcIdx in 0 until numServiceNodes) {
            if (!isPollingStation[svcIdx]) continue

            val nodeIdx = serviceNodes[svcIdx]
            val stationIdx = serviceStations[svcIdx]
            val node = sn.nodes[nodeIdx]

            // Get NodeParam for this station
            val nodeParam = sn.nodeparam?.get(node)
            if (nodeParam is QueueNodeParam) {
                // Set polling type
                pollingType[svcIdx] = nodeParam.pollingType ?: PollingType.EXHAUSTIVE

                // Set K for K-LIMITED
                pollingK[svcIdx] = nodeParam.pollingPar ?: 1

                // Initialize switchover time generators
                for (k in 0 until numClasses) {
                    val switchoverDist = nodeParam.switchoverTime?.get(sn.jobclasses[k])
                    if (switchoverDist != null) {
                        val rate = switchoverDist.rate
                        if (rate > 0 && rate < Double.MAX_VALUE) {
                            val stream = MRG32k3a()
                            stream.setSeed(longArrayOf(seed + 5000 + svcIdx * 100L + k, seed + 5001, seed + 5002, seed + 5003, seed + 5004, seed + 5005))
                            pollingSwitchoverGens[svcIdx][k] = ExponentialGen(stream, rate)
                        }
                    }
                }
            } else {
                // Default to EXHAUSTIVE if no param specified
                pollingType[svcIdx] = PollingType.EXHAUSTIVE
            }
        }
    }

    /**
     * Initialize load-dependent service support.
     * Reads sn.lldscaling matrix and sets up per-station scaling factors.
     * lldscaling[station, n-1] gives the scaling factor when n jobs are at that station.
     */
    private fun initializeLoadDependentService() {
        val hasLld = sn.lldscaling != null && !sn.lldscaling.isEmpty && sn.lldscaling.numCols > 0

        lldScaling = arrayOfNulls(numServiceNodes)
        isLoadDependent = BooleanArray(numServiceNodes) { false }

        if (!hasLld) {
            return  // No load dependence in this model
        }

        for (svcIdx in 0 until numServiceNodes) {
            val stationIdx = serviceStations[svcIdx]
            val maxN = sn.lldscaling.numCols

            // Check if this station has any non-trivial scaling (not all 1s)
            var hasScaling = false
            for (n in 0 until maxN) {
                val scale = sn.lldscaling.get(stationIdx, n)
                if (scale != 1.0 && scale > 0) {
                    hasScaling = true
                    break
                }
            }

            if (hasScaling) {
                isLoadDependent[svcIdx] = true
                lldScaling[svcIdx] = DoubleArray(maxN) { n ->
                    val scale = sn.lldscaling.get(stationIdx, n)
                    if (scale > 0) scale else 1.0
                }
            }
        }
    }

    /**
     * Initialize limited joint-dependent (LJD) scaling tables.
     * LJD uses a single scaling table per station indexed by population vector.
     */
    private fun initializeLjdScaling() {
        // Check if LJD is used
        if (sn.ljdscaling == null || sn.ljdscaling.isEmpty()) {
            hasLjd = false
            return
        }

        ljdScaling = arrayOfNulls(numServiceNodes)
        ljdCutoffs = arrayOfNulls(numServiceNodes)

        for (svcIdx in 0 until numServiceNodes) {
            val stationIdx = serviceStations[svcIdx]
            val station = sn.stations[stationIdx]
            val scaling = sn.ljdscaling[station]
            val cutoffs = sn.ljdcutoffs?.get(station)

            if (scaling != null && !scaling.isEmpty && cutoffs != null) {
                hasLjd = true
                // Convert Matrix to DoubleArray
                ljdScaling!![svcIdx] = DoubleArray(scaling.length()) { i -> scaling.get(i) }
                // Convert cutoffs Matrix to IntArray
                ljdCutoffs!![svcIdx] = IntArray(cutoffs.length()) { i -> cutoffs.get(i).toInt() }
            }
        }
    }

    /**
     * Initialize limited joint class-dependent (LJCD) scaling tables.
     * LJCD uses per-class scaling tables indexed by population vector.
     */
    private fun initializeLjcdScaling() {
        // Check if LJCD is used
        if (sn.ljcdscaling == null || sn.ljcdscaling.isEmpty()) {
            hasLjcd = false
            return
        }

        ljcdScaling = Array(numServiceNodes) { arrayOfNulls<DoubleArray>(numClasses) }
        ljcdCutoffs = arrayOfNulls(numServiceNodes)

        for (svcIdx in 0 until numServiceNodes) {
            val stationIdx = serviceStations[svcIdx]
            val station = sn.stations[stationIdx]
            val stationScaling = sn.ljcdscaling[station]
            val cutoffs = sn.ljcdcutoffs?.get(station)

            if (stationScaling != null && cutoffs != null) {
                hasLjcd = true
                // Convert cutoffs Matrix to IntArray
                ljcdCutoffs!![svcIdx] = IntArray(cutoffs.length()) { i -> cutoffs.get(i).toInt() }

                // Convert per-class scaling tables
                for (classIdx in 0 until numClasses) {
                    val jobClass = sn.jobclasses[classIdx]
                    val classScaling = stationScaling[jobClass]
                    if (classScaling != null && !classScaling.isEmpty) {
                        ljcdScaling!![svcIdx]!![classIdx] = DoubleArray(classScaling.length()) { i -> classScaling.get(i) }
                    }
                }
            }
        }
    }

    /**
     * Compute linearized index from population vector (0-indexed).
     * This is the DES version of ljd_linearize from the MVA solver.
     *
     * @param nvec Per-class population counts
     * @param cutoffs Per-class cutoffs
     * @return Linearized index
     */
    private fun ljdLinearize(nvec: IntArray, cutoffs: IntArray): Int {
        var idx = 0
        var multiplier = 1
        for (k in nvec.indices) {
            val nk = minOf(nvec[k], cutoffs[k])
            idx += nk * multiplier
            multiplier *= (cutoffs[k] + 1)
        }
        return idx
    }

    /**
     * Initialize trace logging if verbose level is DEBUG.
     * Creates a CSV file with detailed simulation event trace in a temporary directory.
     */
    private fun initializeTracing() {
        if (options.verbose == VerboseLevel.DEBUG) {
            traceEnabled = true
            val timestamp = System.currentTimeMillis()
            val tempDir = lineTempName("ssj")
            val traceFile = File(tempDir, "ssj_trace_$timestamp.csv")
            traceWriter = PrintWriter(traceFile)
            // Write CSV header
            traceWriter!!.println("time,event_type,station_idx,class_id,queue_length,busy_servers,phase")
        }
    }

    /**
     * Log a simulation event to the trace file (only if DEBUG verbose level).
     *
     * @param eventType Type of event (ARRIVAL, DEPARTURE, DROP, etc.)
     * @param stationIdx Station index where event occurred
     * @param classId Job class ID
     * @param queueLength Current queue length at station for class
     * @param busyServers Current number of busy servers for class
     */
    private fun logEvent(eventType: String, stationIdx: Int, classId: Int, queueLength: Int, busyServers: Int) {
        if (traceEnabled && traceWriter != null) {
            val phase = if (warmupDone) "STEADY_STATE" else "WARMUP"
            traceWriter!!.println("${Sim.time()},$eventType,$stationIdx,$classId,$queueLength,$busyServers,$phase")
        }
    }

    /**
     * Close the trace writer and flush remaining data.
     */
    private fun closeTracing() {
        if (traceWriter != null) {
            traceWriter!!.flush()
            traceWriter!!.close()
            traceWriter = null
        }
    }

    // ==================== Logger Node Support ====================

    /**
     * Initialize Logger nodes by reading configuration from network structure
     * and creating CSV file writers with JMT-compatible format.
     */
    private fun initializeLoggers() {
        if (loggerNodes.isEmpty()) return

        // Record simulation start time
        val dateFormat = java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss")
        simulationStartTime = dateFormat.format(java.util.Date())

        for (loggerNodeIdx in loggerNodes) {
            // Get Logger node from model
            val loggerNode = sn.nodes[loggerNodeIdx]
            if (loggerNode !is jline.lang.nodes.Logger) continue

            val fileName = loggerNode.fileName ?: "logger_${loggerNode.name}.csv"
            val filePath = loggerNode.filePath ?: ""
            val loggerName = loggerNode.name ?: "Logger"

            val config = LoggerConfig(
                nodeIdx = loggerNodeIdx,
                fileName = fileName,
                filePath = filePath,
                logLoggerName = loggerNode.loggerName,
                logTimestamp = loggerNode.timestamp,
                logJobID = loggerNode.jobID,
                logJobClass = loggerNode.jobClass,
                logTimeSameClass = loggerNode.timeSameClass,
                logTimeAnyClass = loggerNode.timeAnyClass,
                logStartTime = loggerNode.startTime,
                loggerName = loggerName
            )
            loggerConfigs[loggerNodeIdx] = config

            // Initialize last job time tracking
            loggerLastJobTimePerClass[loggerNodeIdx] = DoubleArray(numClasses) { 0.0 }
            loggerLastJobTimeAny[loggerNodeIdx] = 0.0

            // Create output file
            try {
                val outputPath = if (filePath.isNotEmpty()) {
                    File(filePath, fileName)
                } else {
                    File(fileName)
                }
                // Ensure parent directory exists
                outputPath.parentFile?.mkdirs()

                val writer = java.io.BufferedWriter(java.io.FileWriter(outputPath, false))
                loggerWriters[loggerNodeIdx] = writer

                // Write CSV header (JMT format - fixed column order)
                writer.write("LOGGERNAME,TIMESTAMP,JOB_ID,CLASS_ID,INTERARRIVAL_SAMECLASS,INTERARRIVAL_ANYCLASS,SIMUL_START_TIME")
                writer.newLine()
            } catch (e: Exception) {
                line_warning("solver_ssj", "Failed to create logger file for %s: %s", loggerName, e.message)
            }
        }
    }

    /**
     * Initialize routing strategies for each node and class.
     * Sets up state needed for RROBIN, WRROBIN, JSQ, and KCHOICES routing.
     * Also initializes class switch matrices for ClassSwitch nodes.
     */
    private fun initializeRouting() {
        val I = numNodes
        val R = numClasses

        // Initialize routing strategies from sn.routing
        nodeRoutingStrategies = Array(I) { nodeIdx ->
            Array(R) { classIdx ->
                val node = sn.nodes.getOrNull(nodeIdx)
                if (node != null && sn.routing.containsKey(node)) {
                    val classRouting = sn.routing[node]
                    val jobClass = sn.jobclasses.getOrNull(classIdx)
                    if (jobClass != null && classRouting != null) {
                        classRouting[jobClass]
                    } else null
                } else null
            }
        }

        // Initialize round-robin counters for each node
        roundRobinCounters = IntArray(I) { 0 }

        // Initialize WRROBIN weights and RROBIN/WRROBIN outlinks from nodeparam
        wrrobinWeights = Array(I) { Array(R) { null as DoubleArray? } }
        rroutlinks = Array(I) { Array(R) { null as IntArray? } }

        for (nodeIdx in 0 until I) {
            val node = sn.nodes.getOrNull(nodeIdx) ?: continue
            val param = sn.nodeparam?.get(node) ?: continue

            for (classIdx in 0 until R) {
                val jobClass = sn.jobclasses.getOrNull(classIdx) ?: continue
                val routingStrategy = nodeRoutingStrategies.getOrNull(nodeIdx)?.getOrNull(classIdx)

                // Get outlinks for RROBIN/WRROBIN (ordered list of destination node indices)
                if (routingStrategy == RoutingStrategy.RROBIN || routingStrategy == RoutingStrategy.WRROBIN) {
                    val outlinkMatrix = param.outlinks?.get(jobClass)
                    if (outlinkMatrix != null && outlinkMatrix.length() > 0) {
                        val outlinks = IntArray(outlinkMatrix.length().toInt()) { i ->
                            outlinkMatrix.get(i).toInt()
                        }
                        rroutlinks[nodeIdx][classIdx] = outlinks
                    }
                }

                // Get weights for WRROBIN
                if (routingStrategy == RoutingStrategy.WRROBIN) {
                    val weightMatrix = param.weights?.get(jobClass)
                    if (weightMatrix != null) {
                        // weights is a 1xN matrix where weights[0, destNodeIdx] is the weight
                        val weights = DoubleArray(I) { destIdx ->
                            if (destIdx < weightMatrix.getNumCols()) {
                                weightMatrix.get(0, destIdx)
                            } else {
                                0.0
                            }
                        }
                        wrrobinWeights[nodeIdx][classIdx] = weights
                    }
                }
            }
        }

        // Initialize class switch matrices (for potential future use, though class switching
        // is already encoded in the routing matrix by LINE)
        classSwitchMatrices = arrayOfNulls<Array<DoubleArray>?>(I)
        for (csNodeIdx in classSwitchNodes) {
            val node = sn.nodes.getOrNull(csNodeIdx) ?: continue
            val csNode = node as? jline.lang.nodes.ClassSwitch ?: continue
            val csServer = csNode.server as? jline.lang.sections.ClassSwitcher ?: continue

            // Build class switch probability matrix
            val csMatrix = Array(R) { fromClass ->
                DoubleArray(R) { toClass ->
                    csServer.applyCsFun(fromClass, toClass)
                }
            }
            classSwitchMatrices[csNodeIdx] = csMatrix
        }
    }

    /**
     * Log a job passage through a Logger node.
     * Writes a CSV row in JMT-compatible format.
     *
     * @param loggerNodeIdx Logger node index
     * @param classId Job class ID
     * @param jobId Unique job identifier
     */
    private fun logJobPassage(loggerNodeIdx: Int, classId: Int, jobId: Long) {
        val config = loggerConfigs[loggerNodeIdx] ?: return
        val writer = loggerWriters[loggerNodeIdx] ?: return

        val currentTime = Sim.time()
        val delimiter = config.delimiter

        // Calculate interarrival times
        val lastTimePerClass = loggerLastJobTimePerClass[loggerNodeIdx]
        val lastTimeAny = loggerLastJobTimeAny[loggerNodeIdx] ?: 0.0

        val interarrivalSameClass = currentTime - (lastTimePerClass?.get(classId) ?: 0.0)
        val interarrivalAnyClass = currentTime - lastTimeAny

        // Update last job times
        lastTimePerClass?.set(classId, currentTime)
        loggerLastJobTimeAny[loggerNodeIdx] = currentTime

        // Get class name
        val className = if (classId >= 0 && classId < sn.jobclasses.size) {
            sn.jobclasses[classId].name
        } else {
            classId.toString()
        }

        // Build CSV row in JMT format (fixed column order):
        // LOGGERNAME, TIMESTAMP, JOB_ID, CLASS_ID, INTERARRIVAL_SAMECLASS, INTERARRIVAL_ANYCLASS, SIMUL_START_TIME
        val row = StringBuilder()

        // LOGGERNAME
        if (config.logLoggerName) row.append(config.loggerName)
        row.append(delimiter)

        // TIMESTAMP
        if (config.logTimestamp) row.append(formatNumber(currentTime, config.decimalSeparator))
        row.append(delimiter)

        // JOB_ID
        if (config.logJobID) row.append(jobId)
        row.append(delimiter)

        // CLASS_ID
        if (config.logJobClass) row.append(className)
        row.append(delimiter)

        // INTERARRIVAL_SAMECLASS
        if (config.logTimeSameClass) row.append(formatNumber(interarrivalSameClass, config.decimalSeparator))
        row.append(delimiter)

        // INTERARRIVAL_ANYCLASS
        if (config.logTimeAnyClass) row.append(formatNumber(interarrivalAnyClass, config.decimalSeparator))
        row.append(delimiter)

        // SIMUL_START_TIME
        if (config.logStartTime) row.append(simulationStartTime)

        try {
            writer.write(row.toString())
            writer.newLine()
        } catch (e: Exception) {
            // Silently ignore write errors to avoid disrupting simulation
        }
    }

    /**
     * Format a number with the specified decimal separator.
     */
    private fun formatNumber(value: Double, decimalSeparator: String): String {
        val formatted = value.toString()
        return if (decimalSeparator != ".") {
            formatted.replace(".", decimalSeparator)
        } else {
            formatted
        }
    }

    /**
     * Close all Logger file writers.
     */
    private fun closeLoggers() {
        for ((_, writer) in loggerWriters) {
            try {
                writer.flush()
                writer.close()
            } catch (e: Exception) {
                // Ignore close errors
            }
        }
        loggerWriters.clear()
    }

    /**
     * Route a job through pass-through nodes (Logger, Router, ClassSwitch, Cache) until reaching a service node or sink.
     * Logs job passage at each Logger node encountered.
     * Router and ClassSwitch nodes are traversed instantly.
     * ClassSwitch nodes may change the job's class.
     * Cache nodes simulate cache access and determine hit/miss class switching.
     *
     * @param fromNode Starting node
     * @param classId Job class ID
     * @param jobId Unique job identifier
     * @return RoutingResult with final destination node and (possibly switched) class
     */
    private fun routeThroughPassthroughNodes(fromNode: Int, classId: Int, jobId: Long): RoutingResult {
        var currentNode = fromNode
        var currentClass = classId
        var maxIterations = 100  // Prevent infinite loops

        while (maxIterations > 0) {
            val result = selectDestinationWithClassSwitch(currentNode, currentClass)

            if (result.destNode < 0) {
                return RoutingResult(-1, currentClass)  // No valid destination
            }

            if (loggerNodes.contains(result.destNode)) {
                // Log job passage at this Logger node
                logJobPassage(result.destNode, result.destClassId, jobId)
                currentNode = result.destNode
                currentClass = result.destClassId
                maxIterations--
            } else if (routerNodes.contains(result.destNode)) {
                // Router is pass-through (no logging, instant routing)
                currentNode = result.destNode
                currentClass = result.destClassId
                maxIterations--
            } else if (classSwitchNodes.contains(result.destNode)) {
                // ClassSwitch node: the rtnodes matrix encodes class switching
                // Entry format: rtnodes[CS*K+inputClass, dest*K+outputClass]
                // Keep inputClass so next iteration's lookup finds correct entries
                currentNode = result.destNode
                currentClass = result.destClassId
                maxIterations--
            } else if (cacheNodes.contains(result.destNode)) {
                // Cache node: process cache access and determine hit/miss class
                val newClass = processCacheAccess(result.destNode, result.destClassId)
                currentNode = result.destNode
                currentClass = newClass
                maxIterations--
            } else {
                // Reached non-pass-through destination (Queue, Delay, Sink)
                return result
            }
        }

        return RoutingResult(-1, currentClass)  // Max iterations exceeded
    }

    /**
     * Apply class switch at a ClassSwitch node.
     * Uses the csMatrix to determine the new class based on probabilities.
     *
     * @param csNodeIdx ClassSwitch node index
     * @param inClass Input class ID
     * @return Output class ID after switching
     */
    private fun applyClassSwitch(csNodeIdx: Int, inClass: Int): Int {
        val csMatrix = classSwitchMatrices[csNodeIdx]
        if (csMatrix == null) return inClass  // No switch matrix, keep same class

        val probs = csMatrix[inClass]
        val cumProb = probs.sum()
        if (cumProb <= 0) return inClass  // No valid switch probabilities

        val rand = routingRng.nextDouble() * cumProb
        var cumulative = 0.0

        for (toClass in probs.indices) {
            cumulative += probs[toClass]
            if (rand <= cumulative) {
                return toClass
            }
        }

        return inClass  // Fallback: keep same class
    }

    /**
     * Select destination node and class, supporting class switching via routing matrix.
     * For self-looping classes, always returns the reference station node with same class.
     * Supports PROB, RAND, RROBIN, WRROBIN, JSQ, and KCHOICES routing strategies.
     */
    private fun selectDestinationWithClassSwitch(fromNode: Int, classId: Int): RoutingResult {
        // Check if this is a self-looping class - always route to reference station
        if (sn.isslc != null && sn.isslc.get(classId) == 1.0) {
            val refStationIdx = referenceStation[classId]
            return RoutingResult(sn.stationToNode[refStationIdx].toInt(), classId)
        }

        val R = numClasses
        val I = numNodes

        // Get valid destinations considering class switching
        // destination = (nodeIdx, classIdx) pair encoded as RoutingResult
        val destinations = mutableListOf<RoutingResult>()
        val destProbs = mutableListOf<Double>()

        for (toNode in 0 until I) {
            for (toClass in 0 until R) {
                val prob = sn.rtnodes[fromNode * R + classId, toNode * R + toClass]
                if (prob > 0) {
                    destinations.add(RoutingResult(toNode, toClass))
                    destProbs.add(prob)
                }
            }
        }

        if (destinations.isEmpty()) {
            return RoutingResult(-1, classId) // No valid destination
        }

        // If only one destination, return it directly
        if (destinations.size == 1) {
            return destinations[0]
        }

        // Get routing strategy for this node and class
        val routingStrategy = nodeRoutingStrategies.getOrNull(fromNode)?.getOrNull(classId)

        // For JSQ and KCHOICES, we need to compare queue lengths (use node only, ignore class switch for selection)
        return when (routingStrategy) {
            RoutingStrategy.RAND -> selectRandomDestinationWithClass(destinations)
            RoutingStrategy.RROBIN -> selectRoundRobinDestinationWithClass(fromNode, destinations)
            RoutingStrategy.WRROBIN -> selectWeightedRoundRobinDestinationWithClass(fromNode, classId, destinations)
            RoutingStrategy.JSQ -> selectJSQDestinationWithClass(destinations)
            RoutingStrategy.KCHOICES -> selectKChoicesDestinationWithClass(destinations)
            else -> selectProbabilisticDestinationWithClass(destinations, destProbs)  // PROB or default
        }
    }

    /**
     * Select destination node based on routing strategy (direct, without pass-through handling).
     * For self-looping classes, always returns the reference station node.
     * Supports PROB, RAND, RROBIN, WRROBIN, JSQ, and KCHOICES routing strategies.
     * Does NOT support class switching - use selectDestinationWithClassSwitch for that.
     */
    private fun selectDestinationDirect(fromNode: Int, classId: Int): Int {
        // Check if this is a self-looping class - always route to reference station
        if (sn.isslc != null && sn.isslc.get(classId) == 1.0) {
            val refStationIdx = referenceStation[classId]
            return sn.stationToNode[refStationIdx].toInt()
        }

        val R = numClasses
        val I = numNodes

        // Get valid destinations (nodes with positive routing probability) - same class only
        val destinations = mutableListOf<Int>()
        val destProbs = mutableListOf<Double>()

        for (toNode in 0 until I) {
            val prob = sn.rtnodes[fromNode * R + classId, toNode * R + classId]
            if (prob > 0) {
                destinations.add(toNode)
                destProbs.add(prob)
            }
        }

        if (destinations.isEmpty()) {
            return -1 // No valid destination
        }

        // If only one destination, return it directly
        if (destinations.size == 1) {
            return destinations[0]
        }

        // Get routing strategy for this node and class
        val routingStrategy = nodeRoutingStrategies.getOrNull(fromNode)?.getOrNull(classId)

        return when (routingStrategy) {
            RoutingStrategy.RAND -> selectRandomDestination(destinations)
            RoutingStrategy.RROBIN -> selectRoundRobinDestination(fromNode, destinations)
            RoutingStrategy.WRROBIN -> selectWeightedRoundRobinDestination(fromNode, classId, destinations)
            RoutingStrategy.JSQ -> selectJSQDestination(destinations)
            RoutingStrategy.KCHOICES -> selectKChoicesDestination(destinations)
            else -> selectProbabilisticDestination(destinations, destProbs)  // PROB or default
        }
    }

    /**
     * Select destination with class using probabilistic routing (PROB strategy).
     */
    private fun selectProbabilisticDestinationWithClass(destinations: List<RoutingResult>, probs: List<Double>): RoutingResult {
        val cumProb = probs.sum()
        if (cumProb <= 0) return destinations.first()

        val rand = routingRng.nextDouble() * cumProb
        var cumulative = 0.0

        for (i in destinations.indices) {
            cumulative += probs[i]
            if (rand <= cumulative) {
                return destinations[i]
            }
        }

        return destinations.last()
    }

    /**
     * Select destination with class using random routing.
     */
    private fun selectRandomDestinationWithClass(destinations: List<RoutingResult>): RoutingResult {
        val idx = (routingRng.nextDouble() * destinations.size).toInt().coerceIn(0, destinations.size - 1)
        return destinations[idx]
    }

    /**
     * Select destination with class using round-robin routing.
     */
    private fun selectRoundRobinDestinationWithClass(fromNode: Int, destinations: List<RoutingResult>): RoutingResult {
        val idx = roundRobinCounters[fromNode]++ % destinations.size
        return destinations[idx]
    }

    /**
     * Select destination with class using weighted round-robin routing.
     * Uses actual weights from nodeparam instead of rtnodes probabilities.
     */
    private fun selectWeightedRoundRobinDestinationWithClass(fromNode: Int, classId: Int, destinations: List<RoutingResult>): RoutingResult {
        // Get actual weights from wrrobinWeights
        val weights = wrrobinWeights.getOrNull(fromNode)?.getOrNull(classId)

        if (weights == null || destinations.size <= 1) {
            // Fallback to simple round-robin if no weights available
            val idx = roundRobinCounters[fromNode]++ % destinations.size
            return destinations[idx]
        }

        // Build integer weights for each destination
        val scale = 100
        val intWeights = destinations.map { dest ->
            val w = weights.getOrNull(dest.destNode) ?: 0.0
            (w * scale).toInt().coerceAtLeast(if (w > 0) 1 else 0)
        }

        val totalWeight = intWeights.sum()
        if (totalWeight <= 0) {
            // No valid weights, use simple round-robin
            val idx = roundRobinCounters[fromNode]++ % destinations.size
            return destinations[idx]
        }

        // Use weighted round-robin: cycle through based on weights
        val counter = roundRobinCounters[fromNode] % totalWeight
        roundRobinCounters[fromNode]++

        var cumulative = 0
        for (i in destinations.indices) {
            cumulative += intWeights[i]
            if (counter < cumulative) {
                return destinations[i]
            }
        }

        return destinations.last()
    }

    /**
     * Select destination with class using Join Shortest Queue routing.
     */
    private fun selectJSQDestinationWithClass(destinations: List<RoutingResult>): RoutingResult {
        var minQueueLength = Int.MAX_VALUE
        var bestDest = destinations.first()

        for (dest in destinations) {
            val destNode = dest.destNode
            if (destNode < 0 || destNode >= sn.nodeToStation.length()) continue
            val stationIdx = sn.nodeToStation.get(destNode).toInt()
            if (stationIdx < 0 || stationIdx >= currentQueueLength.size) continue

            val queueLen = currentQueueLength[stationIdx].sum()
            if (queueLen < minQueueLength) {
                minQueueLength = queueLen
                bestDest = dest
            }
        }

        return bestDest
    }

    /**
     * Select destination with class using Power of K Choices routing.
     */
    private fun selectKChoicesDestinationWithClass(destinations: List<RoutingResult>): RoutingResult {
        val k = minOf(kchoicesK, destinations.size)
        val candidates = (0 until k).map {
            val idx = (routingRng.nextDouble() * destinations.size).toInt().coerceIn(0, destinations.size - 1)
            destinations[idx]
        }

        var minQueueLength = Int.MAX_VALUE
        var bestDest = candidates.first()

        for (dest in candidates) {
            val destNode = dest.destNode
            if (destNode < 0 || destNode >= sn.nodeToStation.length()) continue
            val stationIdx = sn.nodeToStation.get(destNode).toInt()
            if (stationIdx < 0 || stationIdx >= currentQueueLength.size) continue

            val queueLen = currentQueueLength[stationIdx].sum()
            if (queueLen < minQueueLength) {
                minQueueLength = queueLen
                bestDest = dest
            }
        }

        return bestDest
    }

    /**
     * Select destination using probabilistic routing (PROB strategy).
     * Uses routing probabilities from sn.rtnodes.
     */
    private fun selectProbabilisticDestination(destinations: List<Int>, probs: List<Double>): Int {
        val cumProb = probs.sum()
        if (cumProb <= 0) return destinations.first()

        val rand = routingRng.nextDouble() * cumProb
        var cumulative = 0.0

        for (i in destinations.indices) {
            cumulative += probs[i]
            if (rand <= cumulative) {
                return destinations[i]
            }
        }

        return destinations.last()
    }

    /**
     * Select destination using random uniform routing (RAND strategy).
     */
    private fun selectRandomDestination(destinations: List<Int>): Int {
        val idx = (routingRng.nextDouble() * destinations.size).toInt().coerceIn(0, destinations.size - 1)
        return destinations[idx]
    }

    /**
     * Select destination using round-robin routing (RROBIN strategy).
     * Cycles through destinations in order.
     */
    private fun selectRoundRobinDestination(fromNode: Int, destinations: List<Int>): Int {
        val idx = roundRobinCounters[fromNode] % destinations.size
        roundRobinCounters[fromNode]++
        return destinations[idx]
    }

    /**
     * Select destination using weighted round-robin routing (WRROBIN strategy).
     * Uses actual weights from nodeparam instead of rtnodes probabilities.
     */
    private fun selectWeightedRoundRobinDestination(fromNode: Int, classId: Int, destinations: List<Int>): Int {
        // Get actual weights from wrrobinWeights
        val weights = wrrobinWeights.getOrNull(fromNode)?.getOrNull(classId)

        if (weights == null || destinations.size <= 1) {
            // Fallback to simple round-robin if no weights available
            val idx = roundRobinCounters[fromNode] % destinations.size
            roundRobinCounters[fromNode]++
            return destinations[idx]
        }

        // Build integer weights for each destination
        val scale = 100
        val intWeights = destinations.map { destNode ->
            val w = weights.getOrNull(destNode) ?: 0.0
            (w * scale).toInt().coerceAtLeast(if (w > 0) 1 else 0)
        }

        val totalWeight = intWeights.sum()
        if (totalWeight <= 0) {
            // No valid weights, use simple round-robin
            val idx = roundRobinCounters[fromNode] % destinations.size
            roundRobinCounters[fromNode]++
            return destinations[idx]
        }

        // Use weighted round-robin: cycle through based on weights
        val counter = roundRobinCounters[fromNode] % totalWeight
        roundRobinCounters[fromNode]++

        var cumulative = 0
        for (i in destinations.indices) {
            cumulative += intWeights[i]
            if (counter < cumulative) {
                return destinations[i]
            }
        }

        return destinations.last()
    }

    /**
     * Select destination using Join Shortest Queue routing (JSQ strategy).
     * Routes to the destination with the fewest jobs (queue + in service).
     */
    private fun selectJSQDestination(destinations: List<Int>): Int {
        var minQueue = Int.MAX_VALUE
        var minDest = destinations.first()

        for (destNode in destinations) {
            // Check if destination is a service node (Queue or Delay)
            if (destNode < 0 || destNode >= sn.nodeToStation.length()) continue
            val stationIdx = sn.nodeToStation.get(destNode).toInt()
            if (stationIdx < 0) continue

            val svcIdx = serviceStations.indexOf(stationIdx)

            if (svcIdx >= 0) {
                // Sum queue lengths across all classes at this station
                val totalQueue = currentQueueLength[svcIdx].sum()
                if (totalQueue < minQueue) {
                    minQueue = totalQueue
                    minDest = destNode
                }
            }
        }

        return minDest
    }

    /**
     * Select destination using Power of K Choices routing (KCHOICES strategy).
     * Randomly selects K candidates and routes to the one with shortest queue.
     */
    private fun selectKChoicesDestination(destinations: List<Int>): Int {
        if (destinations.size <= kchoicesK) {
            // If K >= number of destinations, use JSQ
            return selectJSQDestination(destinations)
        }

        // Randomly select K candidates
        val candidates = mutableListOf<Int>()
        val available = destinations.toMutableList()

        repeat(kchoicesK) {
            if (available.isEmpty()) return@repeat
            val idx = (routingRng.nextDouble() * available.size).toInt().coerceIn(0, available.size - 1)
            candidates.add(available.removeAt(idx))
        }

        // Select the one with shortest queue among candidates
        return selectJSQDestination(candidates)
    }

    // ==================== PS (Processor Sharing) Rate Calculations ====================

    /**
     * Calculate effective service rates for all jobs at a PS queue.
     * Dispatches to appropriate rate calculation based on scheduling strategy.
     *
     * @param queueIdx Service node index
     * @param strategy Scheduling strategy
     * @param jobs List of PSCustomer jobs currently in service
     * @param c Number of servers
     * @return Array of effective rates corresponding to each job
     */
    private fun calculatePSRates(
        queueIdx: Int,
        strategy: SchedStrategy,
        jobs: List<PSCustomer>,
        c: Double
    ): DoubleArray {
        val n = jobs.size
        if (n == 0) return DoubleArray(0)

        return when (strategy) {
            SchedStrategy.PS -> calculatePSBasicRates(queueIdx, jobs, c)
            SchedStrategy.DPS -> calculateDPSRates(queueIdx, jobs, c)
            SchedStrategy.GPS -> calculateGPSRates(queueIdx, jobs, c)
            SchedStrategy.PSPRIO -> calculatePSPRIORates(queueIdx, jobs, c)
            SchedStrategy.DPSPRIO -> calculateDPSPRIORates(queueIdx, jobs, c)
            SchedStrategy.GPSPRIO -> calculateGPSPRIORates(queueIdx, jobs, c)
            else -> calculatePSBasicRates(queueIdx, jobs, c)  // Fallback
        }
    }

    /**
     * PS: Equal sharing - each job gets min(1, c/n) of the server capacity.
     * When n <= c, each job gets full rate (1.0); when n > c, jobs share c servers.
     *
     * Returns the "speed factor" for each job: 1.0 means full speed, 0.5 means half speed.
     * The actual time to complete = remainingServiceWork / speedFactor.
     */
    private fun calculatePSBasicRates(
        queueIdx: Int,
        jobs: List<PSCustomer>,
        c: Double
    ): DoubleArray {
        val n = jobs.size
        val rates = DoubleArray(n)
        val sharePerJob = if (n > 0) minOf(1.0, c / n) else 0.0

        for ((idx, job) in jobs.withIndex()) {
            // Rate is the fraction of full speed this job gets (not mu*share)
            rates[idx] = sharePerJob
        }
        return rates
    }

    /**
     * DPS: Discriminatory Processor Sharing - rate proportional to class weight.
     * Weight for class k at station i is in sn.schedparam.
     * Speed factor for class k: w_k * c / sum_j(w_j) where sum is over all jobs.
     *
     * Returns speed factors (not mu*speed), so time to complete = remainingServiceWork / speedFactor.
     */
    private fun calculateDPSRates(
        queueIdx: Int,
        jobs: List<PSCustomer>,
        c: Double
    ): DoubleArray {
        val n = jobs.size
        val rates = DoubleArray(n)
        if (n == 0) return rates

        // Calculate total weighted population: sum(weight_j) for each job
        var totalWeightedPop = 0.0
        for (job in jobs) {
            totalWeightedPop += schedWeights[queueIdx][job.classId]
        }

        if (totalWeightedPop <= 0) {
            // Fallback to equal sharing if weights are zero
            return calculatePSBasicRates(queueIdx, jobs, c)
        }

        // Each job k gets speed factor: w_k * c / totalWeightedPop (capped at 1.0)
        for ((idx, job) in jobs.withIndex()) {
            val weight = schedWeights[queueIdx][job.classId]
            val share = weight * c / totalWeightedPop
            rates[idx] = minOf(1.0, share)
        }
        return rates
    }

    /**
     * GPS: Generalized Processor Sharing.
     * Unlike DPS, GPS divides capacity among CLASS weights, then distributes equally to jobs within each class.
     * Formula: slice(class k) = w_k / Σ(w_i) / n_k
     * where the sum is over all classes with jobs present.
     *
     * Returns speed factors (not mu*speed).
     */
    private fun calculateGPSRates(
        queueIdx: Int,
        jobs: List<PSCustomer>,
        c: Double
    ): DoubleArray {
        val n = jobs.size
        val rates = DoubleArray(n)
        if (n == 0) return rates

        // Count jobs per class
        val jobsPerClass = IntArray(numClasses) { 0 }
        for (job in jobs) {
            jobsPerClass[job.classId]++
        }

        // Calculate total weight sum (over classes with jobs, not jobs)
        var totalWeight = 0.0
        for (k in 0 until numClasses) {
            if (jobsPerClass[k] > 0) {
                totalWeight += schedWeights[queueIdx][k]
            }
        }

        if (totalWeight <= 0) {
            return calculatePSBasicRates(queueIdx, jobs, c)
        }

        // GPS: each job in class k gets w_k / totalWeight / n_k * c
        for ((idx, job) in jobs.withIndex()) {
            val classId = job.classId
            val weight = schedWeights[queueIdx][classId]
            val jobsInClass = jobsPerClass[classId]
            val share = (weight / totalWeight / jobsInClass) * c
            rates[idx] = minOf(1.0, share)
        }
        return rates
    }

    /**
     * PSPRIO: PS with strict priorities - highest priority jobs get service first.
     * Lower priority jobs get rate 0 only when there are more jobs than servers.
     * When n <= c (fewer jobs than servers), all jobs get full rate (like regular PS).
     * Among jobs at the highest priority level, use equal PS sharing.
     *
     * Returns speed factors (not mu*speed).
     */
    private fun calculatePSPRIORates(
        queueIdx: Int,
        jobs: List<PSCustomer>,
        c: Double
    ): DoubleArray {
        val n = jobs.size
        val rates = DoubleArray(n)
        if (n == 0) return rates

        // JMT behavior: when n <= c (fewer jobs than servers), all jobs get full rate
        // Priority only matters when there's contention (n > c)
        if (n <= c) {
            // All jobs get full rate (each job gets its own "virtual server")
            for (idx in 0 until n) {
                rates[idx] = 1.0
            }
            return rates
        }

        // n > c: Apply priority - only highest priority jobs get service
        // Find highest priority among all jobs (LINE/JMT: lower value = higher priority)
        var highestPriority = Int.MAX_VALUE
        for (job in jobs) {
            if (job.priority < highestPriority) {
                highestPriority = job.priority
            }
        }

        // Count jobs at highest priority level
        var highPrioJobCount = 0
        for (job in jobs) {
            if (job.priority == highestPriority) {
                highPrioJobCount++
            }
        }

        // Only highest priority jobs get service
        val sharePerJob = minOf(1.0, c / highPrioJobCount)
        for ((idx, job) in jobs.withIndex()) {
            if (job.priority == highestPriority) {
                rates[idx] = sharePerJob
            } else {
                rates[idx] = 0.0  // Lower priority jobs get nothing
            }
        }
        return rates
    }

    /**
     * DPSPRIO: DPS with strict priorities - highest priority jobs get service first.
     * Lower priority jobs get rate 0 only when there are more jobs than servers.
     * When n <= c (fewer jobs than servers), all jobs get full rate (like regular DPS).
     * Among jobs at the highest priority level, use DPS weighted sharing.
     *
     * Returns speed factors (not mu*speed).
     */
    private fun calculateDPSPRIORates(
        queueIdx: Int,
        jobs: List<PSCustomer>,
        c: Double
    ): DoubleArray {
        val n = jobs.size
        val rates = DoubleArray(n)
        if (n == 0) return rates

        // JMT behavior: when n <= c (fewer jobs than servers), all jobs get full rate
        // Priority only matters when there's contention (n > c)
        if (n <= c) {
            // All jobs get full rate (each job gets its own "virtual server")
            for (idx in 0 until n) {
                rates[idx] = 1.0
            }
            return rates
        }

        // n > c: Apply priority - only highest priority jobs get service
        // Find highest priority among all jobs (LINE/JMT: lower value = higher priority)
        var highestPriority = Int.MAX_VALUE
        for (job in jobs) {
            if (job.priority < highestPriority) {
                highestPriority = job.priority
            }
        }

        // Calculate weighted population for highest priority jobs only
        var totalWeightedPop = 0.0
        for (job in jobs) {
            if (job.priority == highestPriority) {
                totalWeightedPop += schedWeights[queueIdx][job.classId]
            }
        }

        if (totalWeightedPop <= 0) {
            // Fallback to PSPRIO if weights are zero
            return calculatePSPRIORates(queueIdx, jobs, c)
        }

        // Only highest priority jobs get service, using DPS weights
        for ((idx, job) in jobs.withIndex()) {
            if (job.priority == highestPriority) {
                val weight = schedWeights[queueIdx][job.classId]
                val share = weight * c / totalWeightedPop
                rates[idx] = minOf(1.0, share)
            } else {
                rates[idx] = 0.0  // Lower priority jobs get nothing
            }
        }
        return rates
    }

    /**
     * GPSPRIO: GPS with strict priorities - highest priority jobs get service first.
     * Lower priority jobs get rate 0 only when there are more jobs than servers.
     * When n <= c (fewer jobs than servers), all jobs get full rate (like regular GPS).
     * Among jobs at the highest priority level, use GPS weighted sharing (by class, not job count).
     *
     * Returns speed factors (not mu*speed).
     */
    private fun calculateGPSPRIORates(
        queueIdx: Int,
        jobs: List<PSCustomer>,
        c: Double
    ): DoubleArray {
        val n = jobs.size
        val rates = DoubleArray(n)
        if (n == 0) return rates

        // JMT behavior: when n <= c (fewer jobs than servers), all jobs get full rate
        // Priority only matters when there's contention (n > c)
        if (n <= c) {
            // All jobs get full rate (each job gets its own "virtual server")
            for (idx in 0 until n) {
                rates[idx] = 1.0
            }
            return rates
        }

        // n > c: Apply priority - only highest priority jobs get service
        // Find highest priority among all jobs (LINE/JMT: lower value = higher priority)
        var highestPriority = Int.MAX_VALUE
        for (job in jobs) {
            if (job.priority < highestPriority) {
                highestPriority = job.priority
            }
        }

        // Count jobs per class at highest priority level
        val jobsPerClassAtHighPrio = IntArray(numClasses) { 0 }
        for (job in jobs) {
            if (job.priority == highestPriority) {
                jobsPerClassAtHighPrio[job.classId]++
            }
        }

        // Calculate weight sum for highest priority classes only
        var totalWeight = 0.0
        for (k in 0 until numClasses) {
            if (jobsPerClassAtHighPrio[k] > 0) {
                totalWeight += schedWeights[queueIdx][k]
            }
        }

        if (totalWeight <= 0) {
            // Fallback to PSPRIO if weights are zero
            return calculatePSPRIORates(queueIdx, jobs, c)
        }

        // Only highest priority jobs get service, using GPS weights
        for ((idx, job) in jobs.withIndex()) {
            if (job.priority == highestPriority) {
                val classId = job.classId
                val weight = schedWeights[queueIdx][classId]
                val jobsInClass = jobsPerClassAtHighPrio[classId]
                // GPS: class share / number of jobs in class
                val share = (weight / totalWeight / jobsInClass) * c
                rates[idx] = minOf(1.0, share)
            } else {
                rates[idx] = 0.0  // Lower priority jobs get nothing
            }
        }
        return rates
    }

    // ==================== PS Core Functions ====================

    /**
     * Update remaining service work for all jobs at a PS queue based on elapsed time.
     * Must be called before any state change (arrival or departure).
     *
     * @param queueIdx Service node index
     * @param currentTime Current simulation time
     */
    private fun updatePSRemainingWork(queueIdx: Int, currentTime: Double) {
        val elapsedTime = currentTime - psLastUpdateTime[queueIdx]
        if (elapsedTime <= 0) {
            psLastUpdateTime[queueIdx] = currentTime
            return
        }

        val jobs = psJobsInService[queueIdx]
        if (jobs.isEmpty()) {
            psLastUpdateTime[queueIdx] = currentTime
            return
        }

        val strategy = schedStrategies[queueIdx]
        // For load-dependent PS queues, use the effective server count based on current population
        val c = getEffectivePSServerCount(queueIdx)

        // Calculate effective rates for each job
        val rates = calculatePSRates(queueIdx, strategy, jobs, c)

        // Reduce remaining work for each job based on elapsed time and rate
        for ((idx, job) in jobs.withIndex()) {
            val effectiveRate = rates[idx]
            if (effectiveRate > 0) {
                val workDone = elapsedTime * effectiveRate
                job.remainingServiceWork = maxOf(0.0, job.remainingServiceWork - workDone)
            }
        }

        psLastUpdateTime[queueIdx] = currentTime
    }

    /**
     * Cancel all existing departure events and reschedule based on current rates.
     * Called after every arrival and departure at a PS queue.
     */
    private fun rescheduleAllPSDepartures(queueIdx: Int) {
        val jobs = psJobsInService[queueIdx]
        if (jobs.isEmpty()) return

        // Cancel all existing departure events
        for (job in jobs) {
            job.scheduledDepartureEvent?.cancel()
            job.scheduledDepartureEvent = null
        }

        val strategy = schedStrategies[queueIdx]
        // For load-dependent PS queues, use the effective server count based on current population
        val c = getEffectivePSServerCount(queueIdx)
        val rates = calculatePSRates(queueIdx, strategy, jobs, c)

        // Schedule new departures based on remaining work and current rates
        for ((idx, job) in jobs.withIndex()) {
            val rate = rates[idx]

            // Handle numerical edge cases
            if (job.remainingServiceWork <= 1e-12) {
                // Job essentially complete - schedule immediate departure
                val departureEvent = PSDeparture(queueIdx, job)
                departureEvent.schedule(1e-12)  // Tiny delay to maintain event order
                job.scheduledDepartureEvent = departureEvent
            } else if (rate > 0) {
                // Normal case - schedule based on time to complete
                val timeToComplete = job.remainingServiceWork / rate
                val departureEvent = PSDeparture(queueIdx, job)
                departureEvent.schedule(timeToComplete)
                job.scheduledDepartureEvent = departureEvent
            }
            // If rate is 0 (blocked by higher priority), don't schedule departure
        }
    }

    // ==================== PS Arrival Handler ====================

    /**
     * Handle customer arrival at a PS-family queue (PS, DPS, GPS, PSPRIO, DPSPRIO, GPSPRIO).
     * In PS, all customers receive service simultaneously with shared capacity.
     *
     * @param queueIdx Service node index
     * @param customer The arriving customer
     * @param forkedJob Optional forked job info for fork-join tracking
     * @return true if customer was accepted, false if dropped
     */
    private fun arriveAtPSQueue(queueIdx: Int, customer: Customer, forkedJob: ForkedJob? = null): Boolean {
        val classId = customer.classId
        val currentTime = Sim.time()

        // Update remaining work for all current jobs based on elapsed time
        // Must be done BEFORE creating new job, so the new job's remaining work is not decremented
        updatePSRemainingWork(queueIdx, currentTime)

        // Generate service requirement for new customer
        val serviceRequirement = generateServiceTime(queueIdx, classId)

        // Create PS customer
        val psCustomer = PSCustomer(
            classId = classId,
            priority = customer.priority,
            systemArrivalTime = customer.systemArrivalTime,
            queueArrivalTime = customer.queueArrivalTime,
            totalServiceRequirement = serviceRequirement,
            remainingServiceWork = serviceRequirement,
            scheduledDepartureEvent = null,
            forkedJob = forkedJob
        )

        // Update busy time statistics for PS BEFORE adding the new job
        // This accounts for the busy time during the previous state
        updatePSBusyStats(queueIdx)

        // Add to jobs in service (in PS, all jobs are always "in service")
        psJobsInService[queueIdx].add(psCustomer)

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]++

        // Track total customers in service (for internal use, not for busy time)
        customersInService[queueIdx]++

        // Update region job counts if in a region
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]++
            updateRegionArrivalTracking(regionIdx, classId)
        }

        // Reschedule all departures with new rates
        rescheduleAllPSDepartures(queueIdx)

        // Log arrival event
        logEvent("PS_ARRIVAL", serviceStations[queueIdx], classId,
            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])

        return true
    }

    // ==================== PS Departure Event ====================

    /**
     * Departure event for PS-family queues.
     * When a PS customer completes service, we reschedule all remaining jobs.
     */
    private inner class PSDeparture(
        private val queueIdx: Int,
        private val customer: PSCustomer
    ) : Event() {
        override fun actions() {
            val currentTime = Sim.time()
            val classId = customer.classId

            // Update remaining work for all jobs (including this one, should be ~0)
            updatePSRemainingWork(queueIdx, currentTime)

            // Update busy time statistics for PS BEFORE removing the job
            // This accounts for the busy time during the state with this job present
            updatePSBusyStats(queueIdx)

            // Remove departing customer from jobs in service
            psJobsInService[queueIdx].remove(customer)

            // Record queue response time
            val queueResponseTime = currentTime - customer.queueArrivalTime
            responseTimeTally[queueIdx][classId].add(queueResponseTime)
            responseTimeSamples[queueIdx][classId].add(queueResponseTime)
            completedCustomers[queueIdx][classId]++

            // Check event count for stopping/warmup/MSER sampling
            checkEventCountStop()

            // Update queue length statistics
            updateQueueStats(queueIdx, classId)
            currentQueueLength[queueIdx][classId]--

            // Update region job counts if in a region and try releasing blocked customers
            val stationIdx = serviceStations[queueIdx]
            if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
                val regionIdx = fcRegionIndices[stationIdx]
                updateRegionTimeWeightedStats(regionIdx)
                currentJobsInRegion[regionIdx][classId]--
                if (warmupDone) {
                    regionCompletions[regionIdx][classId]++
                }
                // Try to release blocked customers now that space is available
                tryReleaseBlockedCustomers(regionIdx)
            }

            customersInService[queueIdx]--

            // Log departure event
            logEvent("PS_DEPARTURE", stationIdx, classId,
                currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])

            // Reschedule remaining departures with new rates (fewer jobs = faster rates)
            rescheduleAllPSDepartures(queueIdx)

            // Route customer to next destination
            val currentNode = serviceNodes[queueIdx]
            val routingResult = selectDestination(currentNode, classId)
            val destNode = routingResult.destNode
            val destClassId = routingResult.destClassId

            if (destNode >= 0 && sinkNodes.contains(destNode)) {
                // Customer leaves system (count in original class for consistency)
                val systemResponseTime = currentTime - customer.systemArrivalTime
                systemResponseTimeTally[classId].add(systemResponseTime)
                systemCompletedCustomers[classId]++
            } else if (destNode >= 0 && forkNodes.contains(destNode)) {
                // Destination is a Fork node - create forked children
                val parentJobId = nextJobId++
                handleForkArrival(destNode, parentJobId, destClassId, customer.systemArrivalTime)
            } else if (destNode >= 0 && joinNodes.contains(destNode)) {
                // Destination is a Join node - check if this is a forked job
                val forkedJob = customer.forkedJob
                if (forkedJob != null) {
                    val updatedForkedJob = ForkedJob(
                        forkJobId = forkedJob.forkJobId,
                        parentJobId = forkedJob.parentJobId,
                        classId = destClassId,
                        priority = forkedJob.priority,
                        systemArrivalTime = forkedJob.systemArrivalTime,
                        queueArrivalTime = currentTime,
                        randomRank = forkedJob.randomRank
                    )
                    handleJoinArrival(destNode, updatedForkedJob)
                } else {
                    handleUnknownJoinArrival(destNode, destClassId, customer.systemArrivalTime)
                }
            } else if (destNode >= 0 && placeNodes.contains(destNode)) {
                // Destination is a Place node - add token
                handlePlaceArrival(destNode, destClassId, customer.systemArrivalTime)
            } else if (destNode >= 0 && transitionNodes.contains(destNode)) {
                // Destination is a Transition node - check and fire
                checkAndFireTransitions()
            } else if (destNode >= 0) {
                // Route to another queue (possibly with switched class)
                val nextQueueIdx = serviceNodes.indexOf(destNode)
                if (nextQueueIdx >= 0) {
                    val forkedJob = customer.forkedJob
                    val nextCustomer = Customer(destClassId, classPrio[destClassId],
                        customer.systemArrivalTime, currentTime, siroRng.nextDouble())
                    if (forkedJob != null) {
                        // Propagate forked job tracking
                        val nextForkedJob = ForkedJob(
                            forkJobId = forkedJob.forkJobId,
                            parentJobId = forkedJob.parentJobId,
                            classId = destClassId,
                            priority = forkedJob.priority,
                            systemArrivalTime = forkedJob.systemArrivalTime,
                            queueArrivalTime = currentTime,
                            randomRank = forkedJob.randomRank
                        )
                        arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob)
                    } else {
                        arriveAtQueue(nextQueueIdx, nextCustomer)
                    }
                }
            }
        }
    }

    /**
     * External arrival from a source node.
     * For BMAP (Batch Markovian Arrival Process), multiple customers may arrive simultaneously.
     */
    private inner class ExternalArrival(private val sourceIdx: Int, private val classId: Int) : Event() {
        override fun actions() {
            // Get batch size before generating next arrival (BMAP sets this during generateInterarrivalTime)
            val batchSize = arrivalBatchSize[sourceIdx][classId]
            // Reset batch size to 1 (will be updated by next generateInterarrivalTime if BMAP)
            arrivalBatchSize[sourceIdx][classId] = 1

            // Schedule next arrival using appropriate generator (EXP or PH/Erlang/HyperExp/BMAP)
            val interarrivalTime = generateInterarrivalTime(sourceIdx, classId)
            if (interarrivalTime > 0) {
                ExternalArrival(sourceIdx, classId).schedule(interarrivalTime)
            }

            // Process batch: route each customer to destination
            for (b in 0 until batchSize) {
                // Route customer to first destination
                val sourceNode = sourceNodes[sourceIdx]
                val routingResult = selectDestination(sourceNode, classId)
                val destNode = routingResult.destNode
                val destClassId = routingResult.destClassId

                if (destNode >= 0 && sinkNodes.contains(destNode)) {
                    // Customer leaves system immediately (e.g., Source -> Cache -> Sink)
                    // Track completion in original class for consistency with service node departures
                    systemResponseTimeTally[classId].add(0.0)  // Zero response time (immediate departure)
                    systemCompletedCustomers[classId]++
                    checkEventCountStop()
                } else if (destNode >= 0 && forkNodes.contains(destNode)) {
                    // Destination is a Fork node - create forked children
                    val parentJobId = nextJobId++
                    handleForkArrival(destNode, parentJobId, destClassId, Sim.time())
                } else if (destNode >= 0 && placeNodes.contains(destNode)) {
                    // Destination is a Place node - add token
                    handlePlaceArrival(destNode, destClassId, Sim.time())
                } else if (destNode >= 0 && transitionNodes.contains(destNode)) {
                    // Destination is a Transition node - check and fire
                    checkAndFireTransitions()
                } else if (destNode >= 0) {
                    // Find which queue this destination corresponds to
                    val queueIdx = serviceNodes.indexOf(destNode)
                    if (queueIdx >= 0) {
                        // Check if this is a negative signal or catastrophe
                        if ((hasNegativeSignals && isNegativeSignal[destClassId]) ||
                            (hasCatastropheSignals && isCatastropheSignal[destClassId])) {
                            handleNegativeSignalArrival(queueIdx, destClassId, Sim.time())
                        } else {
                            // Use destClassId (possibly switched from original classId)
                            // Assign unique job ID if this class expects a reply (for synchronous call tracking)
                            val jobId = if (synchCallReplyClass[destClassId] >= 0) nextJobId++ else -1L
                            val currentTime = Sim.time()
                            val customer = Customer(destClassId, classPrio[destClassId], currentTime, currentTime, siroRng.nextDouble(), jobId = jobId,
                                absoluteDeadline = currentTime + classDeadline[destClassId])
                            arriveAtQueue(queueIdx, customer)
                        }
                    }
                }
            }
        }
    }

    /**
     * Customer arrives at a service node (Queue or Delay).
     * Checks both total station capacity (cap) and per-class capacity (classcap).
     * Returns true if customer was accepted, false if dropped due to capacity exceeded.
     *
     * @param skipCapacityCheck If true, skip the capacity check (used for BAS-admitted jobs)
     */
    private fun arriveAtQueue(queueIdx: Int, customer: Customer, skipCapacityCheck: Boolean = false): Boolean {
        val classId = customer.classId

        // Track arrivals for arrival rate calculation (including jobs that will be dropped)
        if (warmupDone) {
            arrivedCustomers[queueIdx][classId]++
        }

        // Check for balking BEFORE joining queue
        if (hasBalkingConfig[queueIdx][classId] && shouldBalk(queueIdx, classId)) {
            if (warmupDone) {
                balkedCustomers[queueIdx][classId]++
            }
            logEvent("BALK", serviceStations[queueIdx], classId,
                currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
            return false  // Customer refuses to join - does not enter system
        }

        if (!skipCapacityCheck) {
            // Check total station capacity (Kendall K = queue + in service)
            val currentTotal = getTotalCustomersAtStation(queueIdx)
            if (currentTotal >= bufferCapacities[queueIdx]) {
                // Total capacity exceeded - check for retrial before dropping
                if (hasRetrialConfig[queueIdx][classId]) {
                    // Move to orbit for retrial
                    val maxAttempts = retrialMaxAttemptsConfig[queueIdx][classId]
                    scheduleRetrial(queueIdx, customer, 1, maxAttempts)
                    logEvent("RETRIAL_CAPACITY", serviceStations[queueIdx], classId,
                        currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
                    return false  // Customer goes to orbit, not to queue
                }
                // No retrial configured - drop customer
                if (warmupDone) {
                    droppedCustomers[queueIdx][classId]++
                }
                logEvent("DROP_CAPACITY", serviceStations[queueIdx], classId,
                    currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
                return false
            }

            // Check per-class capacity constraint (only if explicitly set)
            // classCapacities is Int.MAX_VALUE for classes without explicit per-class limits
            if (classCapacities[queueIdx][classId] < Int.MAX_VALUE) {
                val currentClassCount = currentQueueLength[queueIdx][classId]
                if (currentClassCount >= classCapacities[queueIdx][classId]) {
                    // Per-class capacity exceeded - check for retrial before dropping
                    if (hasRetrialConfig[queueIdx][classId]) {
                        // Move to orbit for retrial
                        val maxAttempts = retrialMaxAttemptsConfig[queueIdx][classId]
                        scheduleRetrial(queueIdx, customer, 1, maxAttempts)
                        logEvent("RETRIAL_CLASS_CAPACITY", serviceStations[queueIdx], classId,
                            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
                        return false
                    }
                    // No retrial configured - drop customer
                    if (warmupDone) {
                        droppedCustomers[queueIdx][classId]++
                    }
                    logEvent("DROP_CLASS_CAPACITY", serviceStations[queueIdx], classId,
                        currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
                    return false
                }
            }
        }

        // Check finite capacity region constraints
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]

            // Check global region capacity
            val currentRegionTotal = getTotalCustomersInRegion(regionIdx)
            val globalMax = fcRegionGlobalMax.getOrNull(regionIdx) ?: Int.MAX_VALUE
            if (currentRegionTotal >= globalMax) {
                // Region global capacity exceeded - check per-class dropRule
                val shouldDrop = fcRegionDropRule.getOrNull(regionIdx)?.getOrNull(classId) ?: true

                if (shouldDrop) {
                    if (warmupDone) {
                        droppedByRegion[regionIdx][classId]++
                    }
                    logEvent("DROP_REGION_CAPACITY", stationIdx, classId,
                        currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
                    return false
                } else {
                    // Blocking behavior: add customer to region's blocked queue
                    fcRegionBlockedQueue[regionIdx].add(BlockedCustomer(customer, queueIdx))

                    // JMT semantics: blocked jobs count in FCR QLen AND queue QLen
                    updateRegionTimeWeightedStats(regionIdx)
                    blockedInRegion[regionIdx][classId]++

                    // Track FCR blocked jobs at destination for queue length calculation
                    updateQueueStats(queueIdx, classId)
                    fcrBlockedAtDest[queueIdx][classId]++

                    // NOTE: Do NOT increment currentJobsInRegion or currentQueueLength here!
                    // Blocked customers are tracked separately via fcrBlockedAtDest for queue QLen.

                    logEvent("BLOCK_REGION_CAPACITY", stationIdx, classId,
                        currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
                    return true  // Customer is blocked but not dropped
                }
            }

            // Check per-class region capacity
            val regionClassMax = fcRegionClassMax.get(regionIdx, classId).toInt()
            if (regionClassMax < Int.MAX_VALUE) {
                val currentRegionClass = getClassJobsInRegion(regionIdx, classId)
                if (currentRegionClass >= regionClassMax) {
                    // Region per-class capacity exceeded - check per-class dropRule
                    val shouldDrop = fcRegionDropRule.getOrNull(regionIdx)?.getOrNull(classId) ?: true

                    if (shouldDrop) {
                        if (warmupDone) {
                            droppedByRegion[regionIdx][classId]++
                        }
                        logEvent("DROP_REGION_CLASS_CAPACITY", stationIdx, classId,
                            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
                        return false
                    } else {
                        // Blocking behavior: add customer to region's blocked queue
                        fcRegionBlockedQueue[regionIdx].add(BlockedCustomer(customer, queueIdx))

                        // JMT semantics: blocked jobs count in FCR QLen AND queue QLen
                        updateRegionTimeWeightedStats(regionIdx)
                        blockedInRegion[regionIdx][classId]++

                        // Track FCR blocked jobs at destination for queue length calculation
                        updateQueueStats(queueIdx, classId)
                        fcrBlockedAtDest[queueIdx][classId]++

                        // NOTE: Do NOT increment currentJobsInRegion or currentQueueLength here!
                        // Blocked customers are tracked separately via fcrBlockedAtDest for queue QLen.

                        logEvent("BLOCK_REGION_CLASS_CAPACITY", stationIdx, classId,
                            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
                        return true  // Customer is blocked but not dropped
                    }
                }
            }
        }

        // Dispatch to PS handler for Processor Sharing scheduling
        if (isPSScheduling(schedStrategies[queueIdx])) {
            return arriveAtPSQueue(queueIdx, customer)
        }

        // Dispatch to preemptive LCFS handler
        if (isPreemptiveScheduling[queueIdx]) {
            return arriveAtPreemptiveLCFSQueue(queueIdx, customer)
        }

        // Dispatch to polling handler
        if (isPollingStation[queueIdx]) {
            return arriveAtPollingQueue(queueIdx, customer)
        }

        // For SJF/LJF, we must generate service time upon arrival to sort the queue
        val strategy = schedStrategies[queueIdx]
        if (strategy == SchedStrategy.SJF || strategy == SchedStrategy.LJF) {
            customer.serviceTime = generateServiceTime(queueIdx, classId)
        }

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]++

        // Track max queue length reached
        val totalAtStation = getTotalCustomersAtStation(queueIdx)
        if (totalAtStation > maxQueueLengthReached) {
            maxQueueLengthReached = totalAtStation
        }

        // Update region job counts if in a region
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]++
            updateRegionArrivalTracking(regionIdx, classId)
        }

        if (isDelayNode[queueIdx]) {
            // Delay node (infinite server): always start service immediately
            customersInService[queueIdx]++

            // Update busy time statistics before incrementing
            updateBusyStats(queueIdx, classId)
            currentBusyServers[queueIdx][classId]++
            // Reset busy time tracking for this job - prevents including warmup time in first job of class
            lastBusyUpdateTime[queueIdx][classId] = Sim.time()

            val serviceTime = generateServiceTime(queueIdx, classId)
            val departureEvent = DelayDeparture(queueIdx, customer)
            departureEvent.schedule(serviceTime)
            // Track for signal-based removal (G-networks)
            if (hasNegativeSignals) {
                val jobId = nextDelayJobId++
                delayJobs[jobId] = DelayJob(queueIdx, customer, departureEvent)
            }
        } else {
            // Queue node: check if server is available (heterogeneous-aware)
            val serverSelection = findFreeServerForClass(queueIdx, classId)
            val freeServer = serverSelection.serverId
            val serverTypeId = serverSelection.serverTypeId
            if (freeServer >= 0) {
                // Check server state to determine action (for setup/delayoff feature)
                if (hasSetupDelayoff[queueIdx]) {
                    val state = serverState[queueIdx][freeServer]
                    when (state) {
                        ServerState.ACTIVE -> {
                            // Server is active and idle, start service immediately
                            startService(queueIdx, freeServer, customer, serverTypeId)
                        }
                        ServerState.DELAYOFF -> {
                            // Server is in delayoff, cancel it and start service immediately
                            cancelDelayoff(queueIdx, freeServer)
                            startService(queueIdx, freeServer, customer, serverTypeId)
                        }
                        ServerState.OFF -> {
                            // Server is off, initiate setup and queue the job
                            startServerSetup(queueIdx, freeServer, classId)
                            waitQueues[queueIdx].add(customer)
                            scheduleRenegingEvent(queueIdx, customer)
                        }
                        ServerState.SETUP -> {
                            // Server is already in setup, job must wait
                            waitQueues[queueIdx].add(customer)
                            scheduleRenegingEvent(queueIdx, customer)
                        }
                    }
                } else {
                    // No setup/delayoff, use original logic with heterogeneous tracking
                    markServerBusy(queueIdx, freeServer, serverTypeId)
                    customersInService[queueIdx]++
                    customer.assignedServerType = serverTypeId

                    // Update busy time statistics before incrementing
                    updateBusyStats(queueIdx, classId)
                    currentBusyServers[queueIdx][classId]++
                    // Reset busy time tracking for this job - prevents including warmup time in first job of class
                    lastBusyUpdateTime[queueIdx][classId] = Sim.time()

                    val serviceTime = if (customer.serviceTime > 0) {
                        customer.serviceTime
                    } else {
                        generateHeteroServiceTime(queueIdx, classId, serverTypeId)
                    }
                    if (serviceTime < 0) {
                        throw RuntimeException("DES: Service time is negative ($serviceTime) for station $queueIdx, class $classId at time ${Sim.time()}")
                    }
                    val departureEvent = Departure(queueIdx, freeServer, customer)
                    departureEvent.schedule(serviceTime)
                    // Track for signal-based removal (G-networks)
                    if (hasNegativeSignals) {
                        inServiceJobs[Pair(queueIdx, freeServer)] = InServiceJob(customer, departureEvent)
                    }
                }
            } else {
                // Join queue (priority queue sorts by priority, then FCFS)
                waitQueues[queueIdx].add(customer)

                // Schedule reneging event if patience is configured (customer abandons if wait too long)
                scheduleRenegingEvent(queueIdx, customer)
            }
        }

        // Log arrival event
        logEvent("ARRIVAL", serviceStations[queueIdx], classId,
            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])

        return true
    }

    /**
     * Handles arrival of a negative signal at a service node (G-network semantics).
     * A negative signal removes one job from the queue (if present) uniformly at random
     * from all jobs (waiting + in service), then continues to its destination.
     *
     * G-network behavior:
     * - If queue has jobs: remove one job uniformly at random
     * - If queue is empty: signal has no effect (passes through)
     * - Signal itself is NOT counted as a job (no queue length impact, no service)
     *
     * @param queueIdx The service node index
     * @param signalClassId The class of the negative signal
     * @param systemArrivalTime When the signal entered the system
     */
    private fun handleNegativeSignalArrival(queueIdx: Int, signalClassId: Int, systemArrivalTime: Double) {
        val stationIdx = serviceStations[queueIdx]

        // Get total jobs at this station (all classes)
        val totalJobs = getTotalCustomersAtStation(queueIdx)

        if (totalJobs == 0) {
            // Empty queue - signal has no effect
            logEvent("SIGNAL_MISS", stationIdx, signalClassId, 0, 0)
        } else {
            // Determine number of jobs to remove
            val numToRemove = if (isCatastropheSignal[signalClassId]) {
                // Catastrophe: remove ALL jobs
                totalJobs
            } else {
                // Sample from removal distribution (default: remove exactly 1)
                val dist = signalRemovalDist[signalClassId]
                if (dist != null) {
                    val sampled = dist.sample(1, siroRng)[0].toInt()
                    minOf(sampled, totalJobs)  // Cannot remove more than present
                } else {
                    1  // Default: remove exactly 1 job
                }
            }

            // Get removal policy (default: RANDOM)
            val policy = signalRemovalPolicy[signalClassId] ?: RemovalPolicy.RANDOM

            // Remove jobs according to policy
            var removedCount = 0
            repeat(numToRemove) {
                if (getTotalCustomersAtStation(queueIdx) > 0) {
                    val removedClassId = removeJobBySignalWithPolicy(queueIdx, policy)
                    if (removedClassId >= 0) {
                        removedCount++
                        logEvent("SIGNAL_KILL", stationIdx, signalClassId,
                            getTotalCustomersAtStation(queueIdx), removedClassId)
                    }
                }
            }

            // Log catastrophe event if it was one
            if (isCatastropheSignal[signalClassId] && removedCount > 0) {
                logEvent("CATASTROPHE", stationIdx, signalClassId, removedCount, 0)
            }
        }

        // Route signal to next destination
        routeSignalToNextDestination(queueIdx, signalClassId, systemArrivalTime)
    }

    /**
     * Handles arrival of a REPLY signal at a queue for synchronous call semantics.
     *
     * REPLY signal behavior (LQN SynchCall semantics):
     * - Looks up pending reply by job ID
     * - If found: unblocks the waiting server, emits a job of the ORIGINAL class
     *   that continues from the blocked server (class-switches REPLY back to original)
     * - The REPLY signal is consumed (disappears) - not routed further
     * - If not found: signal passes through (no effect)
     *
     * @param queueIdx The service node index where signal arrived
     * @param signalClassId The class of the REPLY signal
     * @param replyJobId The job ID this reply is responding to
     * @param systemArrivalTime When the signal entered the system
     */
    private fun handleReplySignalArrival(queueIdx: Int, signalClassId: Int, replyJobId: Long, systemArrivalTime: Double) {
        val stationIdx = serviceStations[queueIdx]

        // Count the Reply signal as a completion at this queue (for throughput tracking)
        // The Reply "visits" this queue to unblock a server, so count it as a throughput event
        if (warmupDone) {
            completedCustomers[queueIdx][signalClassId]++
        }

        // Look up pending reply by job ID
        val pendingReply = pendingReplyMap.remove(replyJobId)

        if (pendingReply == null) {
            // No pending reply for this job - signal passes through
            logEvent("REPLY_MISS", stationIdx, signalClassId, 0, replyJobId.toInt())
            routeSignalToNextDestination(queueIdx, signalClassId, systemArrivalTime)
            return
        }

        val blockedQueueIdx = pendingReply.queueIdx
        val blockedServerId = pendingReply.serverId
        val blockedStationIdx = serviceStations[blockedQueueIdx]
        val originalClassId = pendingReply.originalClassId

        // Calculate and accumulate blocking time
        val blockingDuration = Sim.time() - pendingReply.blockStartTime
        totalBlockingTime[blockedQueueIdx][originalClassId] += blockingDuration

        // Update queue stats before decrementing blocked count (to properly time-weight)
        updateQueueStats(blockedQueueIdx, originalClassId)

        // Unblock server
        serverBlocked[blockedQueueIdx][blockedServerId] = false
        currentBlockedServers[blockedQueueIdx][originalClassId]--

        logEvent("REPLY_UNBLOCK", blockedStationIdx, originalClassId, blockedServerId, replyJobId.toInt())

        // Unblock server and update stats (with heterogeneous tracking)
        markServerIdle(blockedQueueIdx, blockedServerId)
        customersInService[blockedQueueIdx]--

        // Start next waiting customer (if any) - heterogeneous-aware
        if (waitQueues[blockedQueueIdx].isNotEmpty()) {
            val nextCustomer = waitQueues[blockedQueueIdx].poll()
            val nextClassId = nextCustomer.classId

            // Cancel any scheduled reneging event for this customer
            cancelRenegingIfScheduled(blockedQueueIdx, nextCustomer)

            // Find compatible server for next customer
            val serverSelection = findFreeServerForClass(blockedQueueIdx, nextClassId)
            if (serverSelection.serverId >= 0) {
                markServerBusy(blockedQueueIdx, serverSelection.serverId, serverSelection.serverTypeId)
                customersInService[blockedQueueIdx]++
                nextCustomer.assignedServerType = serverSelection.serverTypeId

                // Update busy time for the new customer's class
                updateBusyStats(blockedQueueIdx, nextClassId)
                currentBusyServers[blockedQueueIdx][nextClassId]++

                val serviceTime = if (nextCustomer.serviceTime > 0) {
                    nextCustomer.serviceTime
                } else {
                    generateHeteroServiceTime(blockedQueueIdx, nextClassId, serverSelection.serverTypeId)
                }
                val departureEvent = Departure(blockedQueueIdx, serverSelection.serverId, nextCustomer)
                departureEvent.schedule(serviceTime)

                // Track for signal-based removal (G-networks)
                if (hasNegativeSignals) {
                    inServiceJobs[Pair(blockedQueueIdx, serverSelection.serverId)] = InServiceJob(nextCustomer, departureEvent)
                }

                logEvent("UNBLOCK_SERVICE", blockedStationIdx, nextClassId,
                    currentQueueLength[blockedQueueIdx][nextClassId], currentBusyServers[blockedQueueIdx][nextClassId])
            } else {
                // No compatible server - put customer back
                waitQueues[blockedQueueIdx].add(nextCustomer)
            }
        }

        // The REPLY signal continues routing to its destination (e.g., Delay then Sink).
        // This is different from the previous implementation which tried to reconstitute the
        // original class job - that was wrong because the job flow is:
        // Request@Q1 -> Q2 -> Reply@Q1 -> Delay -> Sink
        // The Reply signal itself carries the job forward after unblocking.
        routeSignalToNextDestination(queueIdx, signalClassId, systemArrivalTime)
    }

    /**
     * Removes one job uniformly at random from all jobs at the station (waiting + in service).
     * This implements G-network semantics where negative signals remove a random job.
     *
     * @param queueIdx The service node index
     * @return The classId of the removed job, or -1 if no job could be removed
     */
    private fun removeJobBySignal(queueIdx: Int): Int {
        val strategy = schedStrategies[queueIdx]

        return when {
            isDelayNode[queueIdx] -> removeJobFromDelayBySignal(queueIdx)
            isPSScheduling(strategy) -> removeJobFromPSBySignal(queueIdx)
            isPreemptiveScheduling[queueIdx] -> removeJobFromPreemptiveBySignal(queueIdx)
            isPollingStation[queueIdx] -> removeJobFromPollingBySignal(queueIdx)
            else -> removeJobFromStandardQueueBySignal(queueIdx)
        }
    }

    /**
     * Removes a job from the station using the specified removal policy.
     * This supports batch removal with RANDOM, FCFS, or LCFS selection policies.
     *
     * @param queueIdx The service node index
     * @param policy The removal policy (RANDOM, FCFS, or LCFS)
     * @return The classId of the removed job, or -1 if no job could be removed
     */
    private fun removeJobBySignalWithPolicy(queueIdx: Int, policy: RemovalPolicy): Int {
        val strategy = schedStrategies[queueIdx]

        // For non-standard queue types, fall back to random removal (policy doesn't apply)
        return when {
            isDelayNode[queueIdx] -> removeJobFromDelayBySignal(queueIdx)
            isPSScheduling(strategy) -> removeJobFromPSBySignal(queueIdx)
            isPreemptiveScheduling[queueIdx] -> removeJobFromPreemptiveBySignal(queueIdx)
            isPollingStation[queueIdx] -> removeJobFromPollingBySignal(queueIdx)
            else -> removeJobFromStandardQueueBySignalWithPolicy(queueIdx, policy)
        }
    }

    /**
     * Removes a job from a standard queue using the specified policy.
     * - RANDOM: Uniform random selection from waiting + in-service jobs
     * - FCFS: Remove oldest job first (first arrived)
     * - LCFS: Remove newest job first (last arrived)
     */
    private fun removeJobFromStandardQueueBySignalWithPolicy(queueIdx: Int, policy: RemovalPolicy): Int {
        val waitingCount = waitQueues[queueIdx].size
        val inServiceCount = customersInService[queueIdx]
        val totalJobs = waitingCount + inServiceCount

        if (totalJobs == 0) return -1

        return when (policy) {
            RemovalPolicy.RANDOM -> {
                // Current behavior: uniform random selection
                val victimIdx = (siroRng.nextDouble() * totalJobs).toInt().coerceIn(0, totalJobs - 1)
                if (victimIdx < waitingCount) {
                    removeFromWaitingQueue(queueIdx, victimIdx)
                } else {
                    removeFromInService(queueIdx, victimIdx - waitingCount)
                }
            }
            RemovalPolicy.FCFS -> {
                // Remove oldest job (first in waiting queue, or oldest in service)
                if (waitingCount > 0) {
                    // Waiting queue is ordered - remove from head for FCFS
                    removeFromWaitingQueueOldest(queueIdx)
                } else {
                    removeFromInServiceOldest(queueIdx)
                }
            }
            RemovalPolicy.LCFS -> {
                // Remove newest job (last in waiting queue, or newest in service)
                if (waitingCount > 0) {
                    removeFromWaitingQueueNewest(queueIdx)
                } else {
                    removeFromInServiceNewest(queueIdx)
                }
            }
        }
    }

    /**
     * Removes the oldest job from the waiting queue (for FCFS removal policy).
     */
    private fun removeFromWaitingQueueOldest(queueIdx: Int): Int {
        val waitList = waitQueues[queueIdx].toMutableList()
        if (waitList.isEmpty()) return -1

        // Find oldest by arrival time
        var oldestIdx = 0
        var oldestTime = waitList[0].queueArrivalTime
        for (i in 1 until waitList.size) {
            if (waitList[i].queueArrivalTime < oldestTime) {
                oldestTime = waitList[i].queueArrivalTime
                oldestIdx = i
            }
        }

        return removeFromWaitingQueue(queueIdx, oldestIdx)
    }

    /**
     * Removes the newest job from the waiting queue (for LCFS removal policy).
     */
    private fun removeFromWaitingQueueNewest(queueIdx: Int): Int {
        val waitList = waitQueues[queueIdx].toMutableList()
        if (waitList.isEmpty()) return -1

        // Find newest by arrival time
        var newestIdx = 0
        var newestTime = waitList[0].queueArrivalTime
        for (i in 1 until waitList.size) {
            if (waitList[i].queueArrivalTime > newestTime) {
                newestTime = waitList[i].queueArrivalTime
                newestIdx = i
            }
        }

        return removeFromWaitingQueue(queueIdx, newestIdx)
    }

    /**
     * Removes the oldest job from in-service (for FCFS removal policy).
     */
    private fun removeFromInServiceOldest(queueIdx: Int): Int {
        var oldestTime = Double.MAX_VALUE
        var oldestServerId = -1

        for (sid in 0 until numServers[queueIdx]) {
            if (serverBusy[queueIdx][sid]) {
                val inServiceJob = inServiceJobs[Pair(queueIdx, sid)]
                if (inServiceJob != null && inServiceJob.customer.systemArrivalTime < oldestTime) {
                    oldestTime = inServiceJob.customer.systemArrivalTime
                    oldestServerId = sid
                }
            }
        }

        return if (oldestServerId >= 0) {
            removeFromInServiceById(queueIdx, oldestServerId)
        } else -1
    }

    /**
     * Removes the newest job from in-service (for LCFS removal policy).
     */
    private fun removeFromInServiceNewest(queueIdx: Int): Int {
        var newestTime = Double.MIN_VALUE
        var newestServerId = -1

        for (sid in 0 until numServers[queueIdx]) {
            if (serverBusy[queueIdx][sid]) {
                val inServiceJob = inServiceJobs[Pair(queueIdx, sid)]
                if (inServiceJob != null && inServiceJob.customer.systemArrivalTime > newestTime) {
                    newestTime = inServiceJob.customer.systemArrivalTime
                    newestServerId = sid
                }
            }
        }

        return if (newestServerId >= 0) {
            removeFromInServiceById(queueIdx, newestServerId)
        } else -1
    }

    /**
     * Removes a job in-service at a specific server ID.
     */
    private fun removeFromInServiceById(queueIdx: Int, serverId: Int): Int {
        if (!serverBusy[queueIdx][serverId]) return -1

        // Get and cancel the departure event
        val inServiceJob = inServiceJobs.remove(Pair(queueIdx, serverId)) ?: return -1
        inServiceJob.departureEvent.cancel()

        val classId = inServiceJob.customer.classId

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]--

        // Update region counts if applicable
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]--
        }

        // Update busy time before decrementing
        updateBusyStats(queueIdx, classId)
        currentBusyServers[queueIdx][classId]--

        // Free server (with heterogeneous tracking using departing customer's assigned type)
        val departingServerType = inServiceJob.customer.assignedServerType
        markServerIdle(queueIdx, serverId)
        customersInService[queueIdx]--

        // Start service for next customer in queue (if any)
        if (waitQueues[queueIdx].isNotEmpty()) {
            val nextCustomer = waitQueues[queueIdx].poll()
            val nextClassId = nextCustomer.classId

            // Cancel any scheduled reneging event for this customer
            cancelRenegingIfScheduled(queueIdx, nextCustomer)

            // Find compatible server for next customer (heterogeneous-aware)
            val serverSelection = findFreeServerForClass(queueIdx, nextClassId)
            if (serverSelection.serverId >= 0) {
                markServerBusy(queueIdx, serverSelection.serverId, serverSelection.serverTypeId)
                customersInService[queueIdx]++
                nextCustomer.assignedServerType = serverSelection.serverTypeId

                lastBusyUpdateTime[queueIdx][nextClassId] = Sim.time()
                currentBusyServers[queueIdx][nextClassId]++

                val serviceTime = if (nextCustomer.serviceTime > 0) {
                    nextCustomer.serviceTime
                } else {
                    generateHeteroServiceTime(queueIdx, nextClassId, serverSelection.serverTypeId)
                }
                val departureEvent = Departure(queueIdx, serverSelection.serverId, nextCustomer)
                departureEvent.schedule(serviceTime)
                inServiceJobs[Pair(queueIdx, serverSelection.serverId)] = InServiceJob(nextCustomer, departureEvent)
            } else {
                // No compatible server available - put customer back in queue
                // This shouldn't happen if the just-freed server was compatible
                waitQueues[queueIdx].add(nextCustomer)
            }
        }

        return classId
    }

    /**
     * Removes a random job from a standard FCFS/LCFS/SIRO queue.
     * Selects uniformly at random from waiting jobs + in-service jobs.
     */
    private fun removeJobFromStandardQueueBySignal(queueIdx: Int): Int {
        val waitingCount = waitQueues[queueIdx].size
        val inServiceCount = customersInService[queueIdx]
        val totalJobs = waitingCount + inServiceCount

        if (totalJobs == 0) return -1

        // Select random victim index
        val victimIdx = (siroRng.nextDouble() * totalJobs).toInt().coerceIn(0, totalJobs - 1)

        return if (victimIdx < waitingCount) {
            // Remove from waiting queue
            removeFromWaitingQueue(queueIdx, victimIdx)
        } else {
            // Remove from in-service (cancel departure event)
            val serverIdx = victimIdx - waitingCount
            removeFromInService(queueIdx, serverIdx)
        }
    }

    /**
     * Removes a job at the given index from the waiting queue.
     */
    private fun removeFromWaitingQueue(queueIdx: Int, victimIdx: Int): Int {
        val waitList = waitQueues[queueIdx].toMutableList()
        if (victimIdx >= waitList.size) return -1

        val victim = waitList[victimIdx]
        waitQueues[queueIdx].remove(victim)

        val classId = victim.classId

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]--

        // Update region counts if applicable
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]--
        }

        return classId
    }

    /**
     * Removes a job in service at the given server index and cancels its departure event.
     */
    private fun removeFromInService(queueIdx: Int, serverIdx: Int): Int {
        // Find the nth busy server
        var busyCount = 0
        var targetServerId = -1
        for (sid in 0 until numServers[queueIdx]) {
            if (serverBusy[queueIdx][sid]) {
                if (busyCount == serverIdx) {
                    targetServerId = sid
                    break
                }
                busyCount++
            }
        }

        if (targetServerId < 0) return -1

        // Get and cancel the departure event
        val inServiceJob = inServiceJobs.remove(Pair(queueIdx, targetServerId)) ?: return -1
        inServiceJob.departureEvent.cancel()

        val classId = inServiceJob.customer.classId

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]--

        // Update region counts if applicable
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]--
        }

        // Update busy time before decrementing
        updateBusyStats(queueIdx, classId)
        currentBusyServers[queueIdx][classId]--

        // Free server (with heterogeneous tracking)
        markServerIdle(queueIdx, targetServerId)
        customersInService[queueIdx]--

        // Start service for next customer in queue (if any) - heterogeneous-aware
        if (waitQueues[queueIdx].isNotEmpty()) {
            val nextCustomer = waitQueues[queueIdx].poll()
            val nextClassId = nextCustomer.classId

            // Cancel any scheduled reneging event for this customer
            cancelRenegingIfScheduled(queueIdx, nextCustomer)

            // Find compatible server for next customer
            val serverSelection = findFreeServerForClass(queueIdx, nextClassId)
            if (serverSelection.serverId >= 0) {
                markServerBusy(queueIdx, serverSelection.serverId, serverSelection.serverTypeId)
                customersInService[queueIdx]++
                nextCustomer.assignedServerType = serverSelection.serverTypeId

                // Reset busy time tracking for the new customer's class to prevent phantom busy time accumulation
                lastBusyUpdateTime[queueIdx][nextClassId] = Sim.time()
                currentBusyServers[queueIdx][nextClassId]++

                val serviceTime = if (nextCustomer.serviceTime > 0) {
                    nextCustomer.serviceTime
                } else {
                    generateHeteroServiceTime(queueIdx, nextClassId, serverSelection.serverTypeId)
                }
                val departureEvent = Departure(queueIdx, serverSelection.serverId, nextCustomer)
                departureEvent.schedule(serviceTime)
                inServiceJobs[Pair(queueIdx, serverSelection.serverId)] = InServiceJob(nextCustomer, departureEvent)
            } else {
                // No compatible server - put customer back
                waitQueues[queueIdx].add(nextCustomer)
            }
        }

        return classId
    }

    /**
     * Removes a random job from a PS (Processor Sharing) queue.
     */
    private fun removeJobFromPSBySignal(queueIdx: Int): Int {
        val jobs = psJobsInService[queueIdx]
        if (jobs.isEmpty()) return -1

        // Update remaining work before removal
        updatePSRemainingWork(queueIdx, Sim.time())

        // Select random victim
        val victimIdx = (siroRng.nextDouble() * jobs.size).toInt().coerceIn(0, jobs.size - 1)
        val victim = jobs[victimIdx]

        // Cancel scheduled departure event
        victim.scheduledDepartureEvent?.cancel()

        // Remove from service
        jobs.removeAt(victimIdx)

        val classId = victim.classId

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]--

        // Update region counts
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]--
        }

        // Update busy stats
        updatePSBusyStats(queueIdx)

        // Reschedule PS departures after state change
        rescheduleAllPSDepartures(queueIdx)

        return classId
    }

    /**
     * Removes a random job from a Delay node (infinite servers).
     */
    private fun removeJobFromDelayBySignal(queueIdx: Int): Int {
        // Filter delay jobs for this queue
        val queueDelayJobs = delayJobs.entries.filter { it.value.queueIdx == queueIdx }
        if (queueDelayJobs.isEmpty()) return -1

        // Select random victim
        val victimIdx = (siroRng.nextDouble() * queueDelayJobs.size).toInt().coerceIn(0, queueDelayJobs.size - 1)
        val victimEntry = queueDelayJobs[victimIdx]
        val victim = victimEntry.value

        // Cancel departure event
        victim.departureEvent.cancel()

        // Remove from tracking
        delayJobs.remove(victimEntry.key)

        val classId = victim.customer.classId

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]--

        // Update region counts
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]--
        }

        // Update busy time before decrementing
        updateBusyStats(queueIdx, classId)
        currentBusyServers[queueIdx][classId]--

        // Decrement customers in service
        customersInService[queueIdx]--

        return classId
    }

    /**
     * Removes a random job from a preemptive LCFS queue.
     */
    private fun removeJobFromPreemptiveBySignal(queueIdx: Int): Int {
        val jobs = preemptiveJobsInService[queueIdx]
        if (jobs.isEmpty()) return -1

        // Select random victim
        val victimIdx = (siroRng.nextDouble() * jobs.size).toInt().coerceIn(0, jobs.size - 1)
        val victim = jobs[victimIdx]

        // Cancel scheduled departure
        victim.scheduledDepartureEvent?.cancel()

        // Remove from service
        jobs.removeAt(victimIdx)

        val classId = victim.classId

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]--

        // Update region counts
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]--
        }

        // Update busy time before decrementing
        updateBusyStats(queueIdx, classId)
        currentBusyServers[queueIdx][classId]--

        // Free the server (with heterogeneous tracking)
        markServerIdle(queueIdx, victim.serverId)
        customersInService[queueIdx]--

        return classId
    }

    /**
     * Removes a random job from a polling queue.
     */
    private fun removeJobFromPollingBySignal(queueIdx: Int): Int {
        // Collect all jobs across all class queues
        val allJobs = mutableListOf<Pair<Int, Customer>>()  // (classId, customer)
        for (k in 0 until numClasses) {
            for (customer in pollingQueues[queueIdx][k]) {
                allJobs.add(Pair(k, customer))
            }
        }

        if (allJobs.isEmpty()) return -1

        // Select random victim
        val victimIdx = (siroRng.nextDouble() * allJobs.size).toInt().coerceIn(0, allJobs.size - 1)
        val (classId, victim) = allJobs[victimIdx]

        // Remove from class queue
        pollingQueues[queueIdx][classId].remove(victim)

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]--

        // Update region counts
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]--
        }

        return classId
    }

    /**
     * Routes a signal to its next destination after processing at the current node.
     */
    private fun routeSignalToNextDestination(queueIdx: Int, signalClassId: Int, systemArrivalTime: Double) {
        val currentNode = serviceNodes[queueIdx]
        val routingResult = selectDestination(currentNode, signalClassId)
        val destNode = routingResult.destNode
        val destClassId = routingResult.destClassId

        when {
            destNode < 0 -> {
                // No destination - signal ends here
            }
            sinkNodes.contains(destNode) -> {
                // Signal reaches sink - record system response time
                if (warmupDone) {
                    val respTime = Sim.time() - systemArrivalTime
                    systemResponseTimeTally[signalClassId].add(respTime)
                    systemCompletedCustomers[signalClassId]++
                }
            }
            forkNodes.contains(destNode) -> {
                // Signal at Fork - treat as regular fork arrival
                val parentJobId = nextJobId++
                handleForkArrival(destNode, parentJobId, destClassId, systemArrivalTime)
            }
            else -> {
                // Route to next queue
                val nextQueueIdx = serviceNodes.indexOf(destNode)
                if (nextQueueIdx >= 0) {
                    if (isNegativeSignal[destClassId] || isCatastropheSignal[destClassId]) {
                        // Continue as negative signal or catastrophe
                        handleNegativeSignalArrival(nextQueueIdx, destClassId, systemArrivalTime)
                    } else {
                        // Signal (including REPLY after unblocking) becomes regular customer
                        // for downstream routing. The REPLY has already done its job (unblocking).
                        // Assign job ID if destination class expects a reply (for synchronous call tracking)
                        val needsJobId = synchCallReplyClass[destClassId] >= 0
                        val customer = Customer(destClassId, classPrio[destClassId], systemArrivalTime, Sim.time(), siroRng.nextDouble(),
                            jobId = if (needsJobId) nextJobId++ else -1L)
                        arriveAtQueue(nextQueueIdx, customer)
                    }
                }
            }
        }
    }

    /**
     * Customer arrives at a preemptive LCFS service node.
     * Handles preemption logic for LCFSPR and LCFSPI scheduling.
     */
    private fun arriveAtPreemptiveLCFSQueue(queueIdx: Int, customer: Customer): Boolean {
        val classId = customer.classId
        val strategy = schedStrategies[queueIdx]

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]++

        // Update region tracking if needed
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]++
            updateRegionArrivalTracking(regionIdx, classId)
        }

        // Track max queue length reached
        val totalAtStation = getTotalCustomersAtStation(queueIdx)
        if (totalAtStation > maxQueueLengthReached) {
            maxQueueLengthReached = totalAtStation
        }

        // For size-based preemptive policies (SRPT, PSJF, FB, LRPT): generate service time at arrival
        if (isSizeBasedPreemptiveScheduling(strategy)) {
            // Generate service time for scheduling decision
            customer.serviceTime = generateServiceTime(queueIdx, classId)

            // Check if we should preempt current job
            val shouldPreempt = shouldPreemptForSizeBasedPolicy(queueIdx, customer, strategy)
            if (shouldPreempt) {
                // Preempt current job
                val victimJob = findJobToPreemptForSizeBasedPolicy(queueIdx, customer, strategy)
                if (victimJob != null) {
                    preemptJob(queueIdx, victimJob)
                    // Start service for arriving customer (inherit victim's server type)
                    startPreemptiveService(queueIdx, customer, victimJob.serverId, victimJob.assignedServerType)
                } else {
                    // No job to preempt, join wait queue
                    waitQueues[queueIdx].add(customer)
                }
            } else {
                // Find free server or join queue (heterogeneous-aware)
                val serverSelection = findFreeServerForClass(queueIdx, classId)
                if (serverSelection.serverId >= 0) {
                    startPreemptiveService(queueIdx, customer, serverSelection.serverId, serverSelection.serverTypeId)
                } else {
                    waitQueues[queueIdx].add(customer)
                }
            }
        } else {
            // LCFS preemptive handling
            // Check if we should preempt current job
            val shouldPreempt = shouldPreemptForLCFS(queueIdx, customer)

            if (shouldPreempt) {
                // Preempt current job
                val victimJob = findJobToPreempt(queueIdx, customer, strategy)
                if (victimJob != null) {
                    preemptJob(queueIdx, victimJob)
                    // Start service for arriving customer (inherit victim's server type)
                    startPreemptiveService(queueIdx, customer, victimJob.serverId, victimJob.assignedServerType)
                } else {
                    // No job to preempt, join wait queue
                    waitQueues[queueIdx].add(customer)
                }
            } else {
                // Find free server or join queue (heterogeneous-aware)
                val serverSelection = findFreeServerForClass(queueIdx, classId)
                if (serverSelection.serverId >= 0) {
                    startPreemptiveService(queueIdx, customer, serverSelection.serverId, serverSelection.serverTypeId)
                } else {
                    waitQueues[queueIdx].add(customer)
                }
            }
        }

        logEvent("ARRIVAL", stationIdx, classId,
            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])

        return true
    }

    /**
     * Determines if an arriving customer should preempt a currently serving job.
     */
    private fun shouldPreemptForLCFS(queueIdx: Int, arrivingCustomer: Customer): Boolean {
        // No free servers AND someone is currently in service (heterogeneous-aware)
        val serverSelection = findFreeServerForClass(queueIdx, arrivingCustomer.classId)
        if (serverSelection.serverId >= 0) return false  // Free compatible server available, no preemption needed

        val jobsInService = preemptiveJobsInService[queueIdx]
        if (jobsInService.isEmpty()) return false  // No one to preempt

        // For LCFS/LCFSPR/LCFSPI: always preempt (new arrival has priority)
        // For priority variants: check priority comparison
        val strategy = schedStrategies[queueIdx]
        if (strategy == SchedStrategy.LCFSPRPRIO || strategy == SchedStrategy.LCFSPIPRIO) {
            // Preempt if there's a job with lower OR EQUAL priority (higher or equal numeric value)
            // LINE/JMT convention: lower priority value = higher priority
            // - Lower priority (higher value): higher priority job always preempts
            // - Same priority: follows LCFSPR/LCFSPI semantics (new arrival preempts)
            return jobsInService.any { it.priority >= arrivingCustomer.priority }
        }

        // Non-priority LCFS variants: always preempt if servers full
        return true
    }

    /**
     * Finds the job to preempt based on scheduling strategy.
     */
    private fun findJobToPreempt(
        queueIdx: Int,
        arrivingCustomer: Customer,
        strategy: SchedStrategy
    ): PreemptiveCustomer? {
        val jobsInService = preemptiveJobsInService[queueIdx]

        return when (strategy) {
            SchedStrategy.LCFSPR, SchedStrategy.LCFSPI -> {
                // For basic LCFS: preempt oldest arrival (first to start service)
                jobsInService.minByOrNull { it.queueArrivalTime }
            }
            SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO -> {
                // For priority variants (LINE/JMT convention: lower value = higher priority):
                // 1. Prefer to preempt a lower priority job (higher numeric value) if one exists
                // 2. If all jobs have same priority as arriving customer, apply LCFSPR/PI semantics
                //    (preempt the job that started service first)
                val lowerPriorityJobs = jobsInService.filter { it.priority > arrivingCustomer.priority }
                if (lowerPriorityJobs.isNotEmpty()) {
                    // Preempt the lowest priority job (highest numeric value); if tied, preempt oldest arrival
                    lowerPriorityJobs.minWithOrNull(compareBy({ -it.priority }, { it.queueArrivalTime }))
                } else {
                    // All jobs have same or higher priority - check for same priority
                    val samePriorityJobs = jobsInService.filter { it.priority == arrivingCustomer.priority }
                    if (samePriorityJobs.isNotEmpty()) {
                        // Same priority: LCFSPR/PI semantics - preempt job that started service first
                        samePriorityJobs.minByOrNull { it.serviceStartTime }
                    } else {
                        // All jobs have higher priority (lower numeric value) - cannot preempt
                        null
                    }
                }
            }
            else -> null
        }
    }

    /**
     * Checks if arriving customer should trigger SRPT preemption.
     * Preempts if arriving job has less remaining work than any in-service job.
     */
    private fun shouldPreemptForSRPT(queueIdx: Int, arrivingCustomer: Customer): Boolean {
        val strategy = schedStrategies[queueIdx]
        if (strategy != SchedStrategy.SRPT && strategy != SchedStrategy.SRPTPRIO) {
            return false
        }

        // Check for free servers first (heterogeneous-aware)
        val serverSelection = findFreeServerForClass(queueIdx, arrivingCustomer.classId)
        if (serverSelection.serverId >= 0) return false

        val jobsInService = preemptiveJobsInService[queueIdx]
        if (jobsInService.isEmpty()) return false

        // For SRPTPRIO: only compare within same or lower priority (same or higher numeric value)
        // LINE/JMT convention: lower priority value = higher priority
        if (strategy == SchedStrategy.SRPTPRIO) {
            val eligibleJobs = jobsInService.filter { it.priority >= arrivingCustomer.priority }
            if (eligibleJobs.isEmpty()) return false

            // Preempt if arriving job has less work than any eligible in-service job
            val jobWithMostWork = eligibleJobs.maxByOrNull { it.remainingServiceWork }
            return jobWithMostWork != null &&
                   arrivingCustomer.serviceTime < jobWithMostWork.remainingServiceWork
        }

        // For SRPT: preempt if arriving job has less work than job with most work
        val jobWithMostWork = jobsInService.maxByOrNull { it.remainingServiceWork }
        return jobWithMostWork != null &&
               arrivingCustomer.serviceTime < jobWithMostWork.remainingServiceWork
    }

    /**
     * Finds the job to preempt for SRPT - job with maximum remaining work.
     */
    private fun findJobToPreemptForSRPT(queueIdx: Int, arrivingCustomer: Customer): PreemptiveCustomer? {
        val strategy = schedStrategies[queueIdx]
        val jobsInService = preemptiveJobsInService[queueIdx]

        if (strategy == SchedStrategy.SRPTPRIO) {
            // Find job with most work within same or lower priority (same or higher numeric value)
            // LINE/JMT convention: lower priority value = higher priority
            return jobsInService
                .filter { it.priority >= arrivingCustomer.priority }
                .maxByOrNull { it.remainingServiceWork }
        }

        // Find job with most remaining work
        return jobsInService.maxByOrNull { it.remainingServiceWork }
    }

    /**
     * Determines if an arriving customer should preempt a currently serving job
     * for size-based preemptive policies (SRPT, PSJF, FB, LRPT).
     */
    private fun shouldPreemptForSizeBasedPolicy(
        queueIdx: Int,
        arrivingCustomer: Customer,
        strategy: SchedStrategy
    ): Boolean {
        // Check for free servers first (heterogeneous-aware)
        val serverSelection = findFreeServerForClass(queueIdx, arrivingCustomer.classId)
        if (serverSelection.serverId >= 0) return false

        val jobsInService = preemptiveJobsInService[queueIdx]
        if (jobsInService.isEmpty()) return false

        val currentTime = Sim.time()

        return when {
            isSRPTScheduling(strategy) -> {
                // SRPT: preempt if arriving job has shorter remaining time than job with most CURRENT remaining work
                val jobWithMostWork = jobsInService.maxByOrNull {
                    it.remainingServiceWork - (currentTime - it.serviceStartTime)
                }
                if (jobWithMostWork != null) {
                    val currentRemaining = jobWithMostWork.remainingServiceWork - (currentTime - jobWithMostWork.serviceStartTime)
                    arrivingCustomer.serviceTime < currentRemaining
                } else false
            }
            isPSJFScheduling(strategy) -> {
                // PSJF: preempt if arriving job has smaller original size
                val jobWithLargestOriginal = jobsInService.maxByOrNull { it.totalServiceRequirement }
                jobWithLargestOriginal != null && arrivingCustomer.serviceTime < jobWithLargestOriginal.totalServiceRequirement
            }
            isFBScheduling(strategy) -> {
                // FB/LAS: preempt if arriving job has less attained service
                // Current attained = stored elapsed + time in current service period
                val jobWithMostProgress = jobsInService.maxByOrNull {
                    it.elapsedServiceTime + (currentTime - it.serviceStartTime)
                }
                if (jobWithMostProgress != null) {
                    val currentAttained = jobWithMostProgress.elapsedServiceTime + (currentTime - jobWithMostProgress.serviceStartTime)
                    currentAttained > 0.0  // New arrivals have 0 attained service
                } else false
            }
            isLRPTScheduling(strategy) -> {
                // LRPT: preempt if arriving job has longer remaining time than job with least CURRENT remaining work
                val jobWithLeastWork = jobsInService.minByOrNull {
                    it.remainingServiceWork - (currentTime - it.serviceStartTime)
                }
                if (jobWithLeastWork != null) {
                    val currentRemaining = jobWithLeastWork.remainingServiceWork - (currentTime - jobWithLeastWork.serviceStartTime)
                    arrivingCustomer.serviceTime > currentRemaining
                } else false
            }
            else -> false
        }
    }

    /**
     * Finds the job to preempt for size-based preemptive policies.
     * Uses CURRENT remaining/elapsed time (accounting for time elapsed since service started).
     */
    private fun findJobToPreemptForSizeBasedPolicy(
        queueIdx: Int,
        arrivingCustomer: Customer,
        strategy: SchedStrategy
    ): PreemptiveCustomer? {
        val jobsInService = preemptiveJobsInService[queueIdx]
        val currentTime = Sim.time()

        return when {
            isSRPTScheduling(strategy) -> {
                // SRPT: preempt job with most CURRENT remaining work
                jobsInService.maxByOrNull { it.remainingServiceWork - (currentTime - it.serviceStartTime) }
            }
            isPSJFScheduling(strategy) -> {
                // PSJF: preempt job with largest original size
                jobsInService.maxByOrNull { it.totalServiceRequirement }
            }
            isFBScheduling(strategy) -> {
                // FB/LAS: preempt job with most CURRENT attained service (most progress)
                jobsInService.maxByOrNull { it.elapsedServiceTime + (currentTime - it.serviceStartTime) }
            }
            isLRPTScheduling(strategy) -> {
                // LRPT: preempt job with least CURRENT remaining work
                jobsInService.minByOrNull { it.remainingServiceWork - (currentTime - it.serviceStartTime) }
            }
            else -> null
        }
    }

    /**
     * Preempts a job currently in service, canceling its departure event
     * and adding it back to the wait queue.
     */
    private fun preemptJob(queueIdx: Int, job: PreemptiveCustomer) {
        val currentTime = Sim.time()
        val strategy = schedStrategies[queueIdx]

        // Cancel scheduled departure event
        job.scheduledDepartureEvent?.cancel()

        // Calculate how much time has elapsed since this job started/resumed service
        val timeInCurrentPeriod = currentTime - job.serviceStartTime

        // Update remaining work based on elapsed time
        val newRemainingWork = maxOf(0.0, job.remainingServiceWork - timeInCurrentPeriod)
        job.remainingServiceWork = newRemainingWork
        job.elapsedServiceTime += timeInCurrentPeriod

        // Remove from in-service list
        preemptiveJobsInService[queueIdx].remove(job)

        // Save complete preemption state for resume
        // Use queue arrival time as stable key (identity hash changes when Customer object is recreated)
        // Note: PSJF, FB/LAS, and LRPT are all preemptive-resume (like SRPT and LCFSPR)
        if (strategy == SchedStrategy.LCFSPR || strategy == SchedStrategy.LCFSPRPRIO ||
            isPSJFScheduling(strategy) || isFBScheduling(strategy) || isLRPTScheduling(strategy) ||
            isSRPTScheduling(strategy)) {
            val customerKey = Triple(queueIdx, job.systemArrivalTime, job.queueArrivalTime)
            // Capture distribution type and phase parameter for PH distribution tracking
            val distType = serviceProcessType[queueIdx][job.classId]
            var phaseParam: Int? = null

            // For Erlang distributions, extract the shape parameter (number of phases)
            if (distType == ProcessType.ERLANG) {
                val gen = serviceGens[queueIdx][job.classId]
                if (gen != null) {
                    try {
                        // Attempt to get the shape parameter from ErlangGen using reflection
                        // ErlangGen has a 'k' field representing the number of phases
                        val kField = gen.javaClass.getDeclaredField("k")
                        kField.isAccessible = true
                        phaseParam = (kField.get(gen) as? Number)?.toInt()
                    } catch (e: Exception) {
                        // If we can't access the field, just use null
                        phaseParam = null
                    }
                }
            }

            // Save residual work, original total (for accurate statistics), elapsed time, and distribution tracking info
            preemptedJobHistory[customerKey] = PreemptionRecord(
                remainingWork = newRemainingWork,
                originalTotal = job.totalServiceRequirement,
                elapsedTime = job.elapsedServiceTime,
                distType = distType,
                phaseParam = phaseParam
            )
        }

        // Add back to wait queue as a regular customer with remaining work preserved
        val preemptedCustomer = Customer(
            job.classId,
            job.priority,
            job.systemArrivalTime,
            job.queueArrivalTime,
            job.randomRank,
            serviceTime = newRemainingWork  // Preserve remaining work for SRPT comparisons
        )
        waitQueues[queueIdx].add(preemptedCustomer)

        // Update statistics (job leaves service but stays in queue)
        updateBusyStats(queueIdx, job.classId)
        currentBusyServers[queueIdx][job.classId]--
        customersInService[queueIdx]--

        // DO NOT update queue length - job stays in system

        logEvent("PREEMPT", serviceStations[queueIdx], job.classId,
            currentQueueLength[queueIdx][job.classId], currentBusyServers[queueIdx][job.classId])
    }

    /**
     * Starts preemptive service for a customer.
     * For new arrivals, generates fresh service time.
     * For resumed jobs, uses saved remaining work (LCFSPR) or resamples (LCFSPI).
     */
    private fun startPreemptiveService(queueIdx: Int, customer: Customer, serverId: Int, serverTypeId: Int = -1) {
        val classId = customer.classId
        val currentTime = Sim.time()
        val strategy = schedStrategies[queueIdx]

        // Mark server busy (with heterogeneous tracking)
        markServerBusy(queueIdx, serverId, serverTypeId)
        customersInService[queueIdx]++

        // Store assigned server type on customer
        customer.assignedServerType = serverTypeId

        // Update busy statistics
        updateBusyStats(queueIdx, classId)
        currentBusyServers[queueIdx][classId]++

        // Determine service time based on whether this is a resumed job
        // Use stable key based on arrival times (matches key used in preemptJob)
        val customerKey = Triple(queueIdx, customer.systemArrivalTime, customer.queueArrivalTime)
        val savedRecord = preemptedJobHistory.remove(customerKey)

        val serviceTime: Double
        val totalRequirement: Double
        val remainingWork: Double
        val elapsedTime: Double

        if (savedRecord != null && (strategy == SchedStrategy.LCFSPR || strategy == SchedStrategy.LCFSPRPRIO ||
                                    isSizeBasedPreemptiveScheduling(strategy))) {
            // Resume with saved state (LCFSPR, SRPT, PSJF, FB, LRPT) - use residual remaining work for scheduling
            // All size-based preemptive policies are preemptive-resume
            // Distribution type and phase parameter are tracked for PH distributions (no correction applied)
            serviceTime = savedRecord.remainingWork
            totalRequirement = savedRecord.originalTotal  // FIXED: use original total for statistics
            remainingWork = savedRecord.remainingWork
            elapsedTime = savedRecord.elapsedTime
        } else {
            // Fresh service time (new arrival or LCFSPI resample) - heterogeneous-aware
            val freshServiceTime = generateHeteroServiceTime(queueIdx, classId, serverTypeId)
            serviceTime = freshServiceTime
            totalRequirement = freshServiceTime
            remainingWork = freshServiceTime
            elapsedTime = 0.0
        }

        // Create preemptive customer record
        val preemptiveCustomer = PreemptiveCustomer(
            classId = classId,
            priority = customer.priority,
            systemArrivalTime = customer.systemArrivalTime,
            queueArrivalTime = customer.queueArrivalTime,
            randomRank = customer.randomRank,
            totalServiceRequirement = totalRequirement,
            remainingServiceWork = remainingWork,
            elapsedServiceTime = elapsedTime,  // Preserve elapsed time from before preemption
            serviceStartTime = currentTime,
            scheduledDepartureEvent = null,
            serverId = serverId,
            assignedServerType = serverTypeId
        )

        // Schedule departure
        val departureEvent = PreemptiveDeparture(queueIdx, serverId, preemptiveCustomer)
        departureEvent.schedule(serviceTime)
        preemptiveCustomer.scheduledDepartureEvent = departureEvent

        // Track in service
        preemptiveJobsInService[queueIdx].add(preemptiveCustomer)
    }

    // =========================================================================
    // POLLING SCHEDULING IMPLEMENTATION
    // =========================================================================

    /**
     * Customer arrives at a polling station.
     * Jobs are added to per-class queues and served according to polling discipline.
     */
    private fun arriveAtPollingQueue(queueIdx: Int, customer: Customer): Boolean {
        val classId = customer.classId

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]++

        // Track max queue length reached
        val totalAtStation = getTotalCustomersAtStation(queueIdx)
        if (totalAtStation > maxQueueLengthReached) {
            maxQueueLengthReached = totalAtStation
        }

        // Update region job counts if in a region
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]++
            updateRegionArrivalTracking(regionIdx, classId)
        }

        // Add to per-class polling queue
        pollingQueues[queueIdx][classId].add(customer)

        // Log arrival event
        logEvent("ARRIVAL", serviceStations[queueIdx], classId,
            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])

        // If server is idle and not in switchover, try to start service
        if (customersInService[queueIdx] == 0 && !pollingInSwitchover[queueIdx]) {
            // Check if the current class has jobs or find next class with jobs
            pollingTryStartService(queueIdx)
        }

        return true
    }

    /**
     * Try to start service at a polling station.
     * Called when server becomes idle or after switchover completes.
     */
    private fun pollingTryStartService(queueIdx: Int) {
        val type = pollingType[queueIdx] ?: PollingType.EXHAUSTIVE

        // For GATED: if this is the start of a visit to a new class, record gate size
        if (type == PollingType.GATED && pollingJobsServedInRound[queueIdx] == 0) {
            val currentClass = pollingCurrentClass[queueIdx]
            pollingGateSize[queueIdx] = pollingQueues[queueIdx][currentClass].size
        }

        // Try to get a job from the current class
        val currentClass = pollingCurrentClass[queueIdx]
        val queue = pollingQueues[queueIdx][currentClass]

        if (pollingCanServeFromCurrentClass(queueIdx, type)) {
            // Serve the next job from current class
            val job = queue.poll()
            if (job != null) {
                pollingStartService(queueIdx, job)
            }
        } else {
            // Need to switch to next class
            pollingInitiateSwitchover(queueIdx)
        }
    }

    /**
     * Check if we can serve another job from the current class based on polling type.
     */
    private fun pollingCanServeFromCurrentClass(queueIdx: Int, type: PollingType): Boolean {
        val currentClass = pollingCurrentClass[queueIdx]
        val queueSize = pollingQueues[queueIdx][currentClass].size

        return when (type) {
            PollingType.EXHAUSTIVE -> queueSize > 0
            PollingType.GATED -> {
                val gateSize = pollingGateSize[queueIdx]
                val served = pollingJobsServedInRound[queueIdx]
                queueSize > 0 && served < gateSize
            }
            PollingType.KLIMITED -> {
                val k = pollingK[queueIdx]
                val served = pollingJobsServedInRound[queueIdx]
                queueSize > 0 && served < k
            }
        }
    }

    /**
     * Start service for a job at a polling station.
     */
    private fun pollingStartService(queueIdx: Int, customer: Customer) {
        val classId = customer.classId

        // Find compatible server for heterogeneous-aware polling (if applicable)
        val serverSelection = findFreeServerForClass(queueIdx, classId)
        val serverId = if (serverSelection.serverId >= 0) serverSelection.serverId else 0
        val serverTypeId = serverSelection.serverTypeId

        // Mark server busy with heterogeneous tracking
        markServerBusy(queueIdx, serverId, serverTypeId)
        customersInService[queueIdx]++
        customer.assignedServerType = serverTypeId

        // Reset busy time tracking for this class to prevent phantom busy time accumulation
        // from when the server was serving a different class
        lastBusyUpdateTime[queueIdx][classId] = Sim.time()
        currentBusyServers[queueIdx][classId]++

        // Generate and schedule service time (heterogeneous-aware)
        val serviceTime = generateHeteroServiceTime(queueIdx, classId, serverTypeId)
        PollingDeparture(queueIdx, customer).schedule(serviceTime)
    }

    /**
     * Initiate switchover to the next class at a polling station.
     */
    private fun pollingInitiateSwitchover(queueIdx: Int) {
        // Find the next class with jobs (round-robin)
        var foundNextClass = false
        var nextClass = pollingCurrentClass[queueIdx]

        for (i in 0 until numClasses) {
            nextClass = (nextClass + 1) % numClasses
            if (pollingQueues[queueIdx][nextClass].isNotEmpty()) {
                foundNextClass = true
                break
            }
        }

        if (!foundNextClass) {
            // No jobs at any class - server stays idle at current class
            return
        }

        // Get switchover time for transitioning to next class
        val switchoverGen = pollingSwitchoverGens[queueIdx][nextClass]
        val switchoverTime = switchoverGen?.nextDouble() ?: 0.0

        if (switchoverTime > 0) {
            // Schedule switchover completion
            pollingInSwitchover[queueIdx] = true
            PollingSwitchover(queueIdx, nextClass).schedule(switchoverTime)
        } else {
            // No switchover time - switch immediately
            pollingCurrentClass[queueIdx] = nextClass
            pollingJobsServedInRound[queueIdx] = 0
            pollingTryStartService(queueIdx)
        }
    }

    /**
     * Departure event for polling station.
     */
    private inner class PollingDeparture(
        private val queueIdx: Int,
        private val customer: Customer
    ) : Event() {
        override fun actions() {
            val classId = customer.classId

            // Record queue response time
            val queueResponseTime = Sim.time() - customer.queueArrivalTime
            responseTimeTally[queueIdx][classId].add(queueResponseTime)
            responseTimeSamples[queueIdx][classId].add(queueResponseTime)
            completedCustomers[queueIdx][classId]++

            // Check event count for stopping/warmup/MSER sampling
            checkEventCountStop()

            // Update queue length statistics
            updateQueueStats(queueIdx, classId)
            currentQueueLength[queueIdx][classId]--

            // Update region job counts if in a region and try releasing blocked customers
            val stationIdx = serviceStations[queueIdx]
            if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
                val regionIdx = fcRegionIndices[stationIdx]
                updateRegionTimeWeightedStats(regionIdx)
                currentJobsInRegion[regionIdx][classId]--
                if (warmupDone) {
                    regionCompletions[regionIdx][classId]++
                }
                // Try to release blocked customers now that space is available
                tryReleaseBlockedCustomers(regionIdx)
            }

            // Update busy time before decrementing
            updateBusyStats(queueIdx, classId)
            currentBusyServers[queueIdx][classId]--

            // Free server (with heterogeneous tracking)
            markServerIdle(queueIdx, 0)
            customersInService[queueIdx]--

            // Increment jobs served in this round
            pollingJobsServedInRound[queueIdx]++

            // Log departure event
            logEvent("DEPARTURE", stationIdx, classId,
                currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])

            // Route customer to next destination
            val currentNode = serviceNodes[queueIdx]
            val routingResult = selectDestination(currentNode, classId)
            val destNode = routingResult.destNode
            val destClassId = routingResult.destClassId

            if (destNode >= 0 && sinkNodes.contains(destNode)) {
                // Customer leaves system
                val systemResponseTime = Sim.time() - customer.systemArrivalTime
                systemResponseTimeTally[classId].add(systemResponseTime)
                systemCompletedCustomers[classId]++
            } else if (destNode >= 0 && forkNodes.contains(destNode)) {
                // Destination is a Fork node
                val parentJobId = nextJobId++
                handleForkArrival(destNode, parentJobId, destClassId, customer.systemArrivalTime)
            } else if (destNode >= 0 && joinNodes.contains(destNode)) {
                // Destination is a Join node
                val forkedJobKey = Pair(queueIdx, customer.queueArrivalTime)
                val forkedJob = forkedCustomerMap.remove(forkedJobKey)
                if (forkedJob != null) {
                    val updatedForkedJob = ForkedJob(
                        forkJobId = forkedJob.forkJobId,
                        parentJobId = forkedJob.parentJobId,
                        classId = destClassId,
                        priority = forkedJob.priority,
                        systemArrivalTime = forkedJob.systemArrivalTime,
                        queueArrivalTime = Sim.time(),
                        randomRank = forkedJob.randomRank
                    )
                    handleJoinArrival(destNode, updatedForkedJob)
                } else {
                    handleUnknownJoinArrival(destNode, destClassId, customer.systemArrivalTime)
                }
            } else if (destNode >= 0) {
                // Route to another queue (possibly with switched class)
                val nextQueueIdx = serviceNodes.indexOf(destNode)
                if (nextQueueIdx >= 0) {
                    val forkedJobKey = Pair(queueIdx, customer.queueArrivalTime)
                    val forkedJob = forkedCustomerMap.remove(forkedJobKey)
                    val nextCustomer = Customer(destClassId, classPrio[destClassId], customer.systemArrivalTime, Sim.time(), siroRng.nextDouble(),
                        absoluteDeadline = customer.absoluteDeadline)
                    if (forkedJob != null) {
                        val nextForkedJob = ForkedJob(
                            forkJobId = forkedJob.forkJobId,
                            parentJobId = forkedJob.parentJobId,
                            classId = destClassId,
                            priority = forkedJob.priority,
                            systemArrivalTime = forkedJob.systemArrivalTime,
                            queueArrivalTime = Sim.time(),
                            randomRank = forkedJob.randomRank
                        )
                        arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob)
                    } else {
                        arriveAtQueue(nextQueueIdx, nextCustomer)
                    }
                }
            }

            // Try to serve next job or initiate switchover
            if (!pollingInSwitchover[queueIdx]) {
                pollingTryStartService(queueIdx)
            }
        }
    }

    /**
     * Switchover event for polling station.
     * Completes transition to the next class.
     */
    private inner class PollingSwitchover(
        private val queueIdx: Int,
        private val nextClass: Int
    ) : Event() {
        override fun actions() {
            pollingInSwitchover[queueIdx] = false
            pollingCurrentClass[queueIdx] = nextClass
            pollingJobsServedInRound[queueIdx] = 0

            // Try to start service in the new class
            pollingTryStartService(queueIdx)
        }
    }

    // =========================================================================
    // END OF POLLING SCHEDULING IMPLEMENTATION
    // =========================================================================

    /**
     * Returns the total number of customers at a station (waiting + in service).
     * Note: currentQueueLength already includes ALL jobs at the station (both waiting and in service).
     */
    private fun getTotalCustomersAtStation(queueIdx: Int): Int {
        var total = 0
        for (k in 0 until numClasses) {
            total += currentQueueLength[queueIdx][k]
        }
        return total
    }

    // ==================== Impatience Helper Functions ====================

    /**
     * Updates time-weighted orbit statistics before changing orbit size.
     */
    private fun updateOrbitTimeStats(queueIdx: Int, classId: Int) {
        val currentTime = Sim.time()
        val elapsed = currentTime - lastOrbitUpdateTime[queueIdx][classId]
        if (elapsed > 0 && warmupDone) {
            totalOrbitTime[queueIdx][classId] += currentOrbitSize[queueIdx][classId] * elapsed
        }
        lastOrbitUpdateTime[queueIdx][classId] = currentTime
    }

    /**
     * Schedules a reneging event for an impatient customer joining a queue.
     * Called when a customer joins a queue with reneging configured.
     */
    private fun scheduleRenegingEvent(queueIdx: Int, customer: Customer) {
        val classId = customer.classId
        if (!hasPatienceConfig[queueIdx][classId]) return

        val patienceGen = patienceGens[queueIdx][classId] ?: return
        val patienceTime = patienceGen.nextDouble()
        if (patienceTime <= 0) return

        val currentTime = Sim.time()
        val deadline = currentTime + patienceTime

        val key = Triple(queueIdx, customer.systemArrivalTime.toBits(), classId)
        val renegingEvent = RenegingEvent(queueIdx, customer, key)
        renegingEvent.schedule(patienceTime)

        val impatient = ImpatientCustomer(customer, patienceTime, deadline, renegingEvent)
        waitingImpatientCustomers[key] = impatient
    }

    /**
     * Cancels a scheduled reneging event when customer starts service.
     */
    private fun cancelRenegingIfScheduled(queueIdx: Int, customer: Customer) {
        val classId = customer.classId
        if (!hasPatienceConfig[queueIdx][classId]) return

        val key = Triple(queueIdx, customer.systemArrivalTime.toBits(), classId)
        val impatient = waitingImpatientCustomers.remove(key)
        impatient?.renegingEvent?.cancel()
    }

    /**
     * Schedules a retrial event for a customer entering the orbit.
     */
    private fun scheduleRetrial(queueIdx: Int, customer: Customer, attempts: Int, maxAttempts: Int) {
        val classId = customer.classId
        val retrialGen = retrialGens[queueIdx][classId] ?: return

        val retrialDelay = retrialGen.nextDouble()
        if (retrialDelay <= 0) return

        val orbitJob = OrbitJob(
            customer = customer,
            orbitEntryTime = Sim.time(),
            retrialAttempts = attempts,
            maxAttempts = maxAttempts,
            sourceQueueIdx = queueIdx
        )

        val retrialEvent = RetrialEvent(queueIdx, orbitJob)
        retrialEvent.schedule(retrialDelay)
        orbitJob.retrialEvent = retrialEvent

        orbitJobs[queueIdx].add(orbitJob)

        // Update orbit statistics
        updateOrbitTimeStats(queueIdx, classId)
        currentOrbitSize[queueIdx][classId]++
    }

    /**
     * Handle customer re-entering queue after successful retrial.
     * Similar to regular arrival but skips balking check (they already decided to join).
     */
    private fun arriveAtQueueFromRetrial(queueIdx: Int, customer: Customer) {
        val classId = customer.classId

        // Update queue length statistics
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]++

        // Track successful retrial
        retriedCustomers[queueIdx][classId]++

        // Try to start service or join queue
        if (isDelayNode[queueIdx]) {
            // Delay node: always start service immediately
            customersInService[queueIdx]++
            updateBusyStats(queueIdx, classId)
            currentBusyServers[queueIdx][classId]++
            lastBusyUpdateTime[queueIdx][classId] = Sim.time()

            val serviceTime = generateServiceTime(queueIdx, classId)
            val departureEvent = DelayDeparture(queueIdx, customer)
            departureEvent.schedule(serviceTime)
        } else {
            val freeServer = findFreeServer(queueIdx)
            if (freeServer >= 0) {
                // Start service immediately
                serverBusy[queueIdx][freeServer] = true
                customersInService[queueIdx]++
                updateBusyStats(queueIdx, classId)
                currentBusyServers[queueIdx][classId]++
                lastBusyUpdateTime[queueIdx][classId] = Sim.time()

                val serviceTime = generateServiceTime(queueIdx, classId)
                val departureEvent = Departure(queueIdx, freeServer, customer)
                departureEvent.schedule(serviceTime)
            } else {
                // Join wait queue
                waitQueues[queueIdx].add(customer)

                // Schedule reneging if configured
                scheduleRenegingEvent(queueIdx, customer)
            }
        }
    }

    /**
     * Determines if an arriving customer should balk based on configured thresholds.
     */
    private fun shouldBalk(queueIdx: Int, classId: Int): Boolean {
        if (!hasBalkingConfig[queueIdx][classId]) return false

        val thresholds = sn.balkingThresholds?.get(sn.stations[serviceStations[queueIdx]])?.get(sn.jobclasses[classId])
            ?: return false

        val currentQueueLen = getTotalCustomersAtStation(queueIdx)

        for (threshold in thresholds) {
            if (threshold.matches(currentQueueLen)) {
                // Within range - balk with given probability
                return routingRng.nextDouble() < threshold.probability
            }
        }
        return false
    }

    /**
     * Checks if BAS (Blocking After Service) is enabled for the given station and class.
     * JMT convention: blocking policy is specified on the DESTINATION queue (the queue with limited capacity).
     * When a job completes service at an upstream queue and the destination is full, the server blocks.
     */
    private fun hasBASBlocking(destQueueIdx: Int, classId: Int): Boolean {
        if (destQueueIdx < 0 || destQueueIdx >= stationDropRule.size) return false
        if (classId < 0 || classId >= stationDropRule[destQueueIdx].size) return false
        return stationDropRule[destQueueIdx][classId] == DropStrategy.BlockingAfterService.getID()
    }

    /**
     * Checks if BBS (Blocking Before Service) is enabled for the given station and class.
     * JMT convention: blocking policy is specified on the DESTINATION queue (the queue with limited capacity).
     */
    private fun hasBBSBlocking(destQueueIdx: Int, classId: Int): Boolean {
        if (destQueueIdx < 0 || destQueueIdx >= stationDropRule.size) return false
        if (classId < 0 || classId >= stationDropRule[destQueueIdx].size) return false
        return stationDropRule[destQueueIdx][classId] == DropStrategy.BlockingBeforeService.getID()
    }

    /**
     * Checks if RSRD (Re-Service on Rejection with Delay) is enabled for the given station and class.
     */
    private fun hasRSRDBlocking(destQueueIdx: Int, classId: Int): Boolean {
        if (destQueueIdx < 0 || destQueueIdx >= stationDropRule.size) return false
        if (classId < 0 || classId >= stationDropRule[destQueueIdx].size) return false
        return stationDropRule[destQueueIdx][classId] == DropStrategy.ReServiceOnRejection.getID()
    }

    /**
     * Checks if a destination queue has capacity to accept a new customer.
     * Returns true if the destination can accept, false if full.
     */
    private fun destinationHasCapacity(destQueueIdx: Int, destClassId: Int): Boolean {
        if (destQueueIdx < 0 || destQueueIdx >= numServiceNodes) return true  // Unknown destination, allow

        // Check total station capacity
        val currentTotal = getTotalCustomersAtStation(destQueueIdx)
        if (currentTotal >= bufferCapacities[destQueueIdx]) {
            return false
        }

        // Check per-class capacity constraint
        if (classCapacities[destQueueIdx][destClassId] < Int.MAX_VALUE) {
            val currentClassCount = currentQueueLength[destQueueIdx][destClassId]
            if (currentClassCount >= classCapacities[destQueueIdx][destClassId]) {
                return false
            }
        }

        return true
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
    private fun tryUnblockBBSServers(destQueueIdx: Int) {
        // Get all sources with blocked servers waiting for this destination
        val sourcesWithBlocked = bbsDestinationToSources[destQueueIdx] ?: return
        if (sourcesWithBlocked.isEmpty()) return

        // Keep trying to unblock servers while destination has capacity
        while (true) {
            // Find oldest blocked server across ALL sources (global FIFO by block start time)
            var oldestBlocked: BBSBlockedServer? = null
            var oldestSourceIdx: Int = -1
            var oldestTimestamp: Double = Double.POSITIVE_INFINITY

            for (sourceIdx in sourcesWithBlocked.toList()) {
                val blockedList = bbsBlockedServers[sourceIdx] ?: continue

                // Find first blocked server for this destination from this source
                for (blocked in blockedList) {
                    if (blocked.destQueueIdx == destQueueIdx) {
                        if (blocked.blockStartTime < oldestTimestamp) {
                            oldestTimestamp = blocked.blockStartTime
                            oldestBlocked = blocked
                            oldestSourceIdx = sourceIdx
                        }
                        break  // Only check first match per source (FIFO per source)
                    }
                }
            }

            // No more blocked servers waiting for this destination
            if (oldestBlocked == null) {
                break
            }

            // Check if destination has capacity for this specific class
            if (!destinationHasCapacity(destQueueIdx, oldestBlocked.destClassId)) {
                break  // No capacity for this class
            }

            // Unblock the oldest blocked server
            val sourceQueueIdx = oldestSourceIdx
            val customer = oldestBlocked.customer
            val destClassId = oldestBlocked.destClassId
            val serverId = oldestBlocked.serverId
            val sourceClassId = oldestBlocked.sourceClassId

            // Remove from blocked servers list
            val blockedList = bbsBlockedServers[sourceQueueIdx]!!
            blockedList.remove(oldestBlocked)

            // Clean up empty lists and indices
            if (blockedList.none { it.destQueueIdx == destQueueIdx }) {
                sourcesWithBlocked.remove(sourceQueueIdx)
            }
            if (blockedList.isEmpty()) {
                bbsBlockedServers.remove(sourceQueueIdx)
            }

            val sourceStationIdx = serviceStations[sourceQueueIdx]
            logEvent("BBS_UNBLOCK", sourceStationIdx, sourceClassId, serverId, destQueueIdx)

            // BBS: blocked job was counted at DESTINATION, so decrement there
            // Update destination queue stats before decrementing blocked count (to properly time-weight)
            updateQueueStats(destQueueIdx, destClassId)
            bbsBlockedAtDest[destQueueIdx][destClassId]--

            // Unblock server at source
            serverBlocked[sourceQueueIdx][serverId] = false

            // Create customer for destination
            val nextCustomer = Customer(
                destClassId,
                classPrio[destClassId],
                customer.systemArrivalTime,
                Sim.time(),  // Queue arrival time at destination is now
                siroRng.nextDouble(),
                jobId = customer.jobId,
                absoluteDeadline = customer.absoluteDeadline
            )

            // Send job to destination queue
            arriveAtQueue(destQueueIdx, nextCustomer, skipCapacityCheck = true)

            // Free server and start next customer at source queue (with heterogeneous tracking)
            markServerIdle(sourceQueueIdx, serverId)
            customersInService[sourceQueueIdx]--

            // Start service for next customer in source queue (if any) - heterogeneous-aware
            if (waitQueues[sourceQueueIdx].isNotEmpty()) {
                val nextWaiting = waitQueues[sourceQueueIdx].poll()
                val nextClassId = nextWaiting.classId

                // Find compatible server for next customer
                val serverSelection = findFreeServerForClass(sourceQueueIdx, nextClassId)
                if (serverSelection.serverId >= 0) {
                    markServerBusy(sourceQueueIdx, serverSelection.serverId, serverSelection.serverTypeId)
                    customersInService[sourceQueueIdx]++
                    nextWaiting.assignedServerType = serverSelection.serverTypeId

                    // Update busy time for the new customer's class
                    lastBusyUpdateTime[sourceQueueIdx][nextClassId] = Sim.time()
                    currentBusyServers[sourceQueueIdx][nextClassId]++

                    val serviceTime = if (nextWaiting.serviceTime > 0) {
                        nextWaiting.serviceTime
                    } else {
                        generateHeteroServiceTime(sourceQueueIdx, nextClassId, serverSelection.serverTypeId)
                    }
                    val departureEvent = Departure(sourceQueueIdx, serverSelection.serverId, nextWaiting)
                    departureEvent.schedule(serviceTime)

                    // Track for signal-based removal (G-networks)
                    if (hasNegativeSignals) {
                        inServiceJobs[Pair(sourceQueueIdx, serverSelection.serverId)] = InServiceJob(nextWaiting, departureEvent)
                    }

                    logEvent("BBS_NEXT_SERVICE", sourceStationIdx, nextClassId,
                        currentQueueLength[sourceQueueIdx][nextClassId], currentBusyServers[sourceQueueIdx][nextClassId])
                } else {
                    // No compatible server - put customer back
                    waitQueues[sourceQueueIdx].add(nextWaiting)
                }
            }
        }

        // Clean up reverse index if no more sources blocked for this destination
        if (sourcesWithBlocked.isEmpty()) {
            bbsDestinationToSources.remove(destQueueIdx)
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
    private fun tryAdmitBASWaitingJobs(destQueueIdx: Int) {
        // Get all sources with waiting jobs for this destination
        val sourcesWithWaiting = basDestinationToSources[destQueueIdx] ?: return
        if (sourcesWithWaiting.isEmpty()) return

        // Keep admitting jobs while destination has capacity
        while (true) {
            // Find oldest waiting job across ALL sources (global FIFO by arrival time)
            var oldestWaiting: BASWaitingJob? = null
            var oldestSourceIdx: Int = -1
            var oldestTimestamp: Double = Double.POSITIVE_INFINITY

            for (sourceIdx in sourcesWithWaiting.toList()) {
                val outgoingBuffer = basOutgoingBuffer[sourceIdx] ?: continue

                // Find first waiting job for this destination from this source
                for (waiting in outgoingBuffer) {
                    if (waiting.destQueueIdx == destQueueIdx) {
                        if (waiting.arrivalTime < oldestTimestamp) {
                            oldestTimestamp = waiting.arrivalTime
                            oldestWaiting = waiting
                            oldestSourceIdx = sourceIdx
                        }
                        break  // Only check first match per source (FIFO per source)
                    }
                }
            }

            // No more waiting jobs for this destination
            if (oldestWaiting == null) {
                break
            }

            // Check if destination has capacity for this specific class
            if (!destinationHasCapacity(destQueueIdx, oldestWaiting.destClassId)) {
                break  // No capacity for this class
            }

            // Remove from source's outgoing buffer
            val sourceQueueIdx = oldestSourceIdx
            val outgoingBuffer = basOutgoingBuffer[sourceQueueIdx]!!
            outgoingBuffer.remove(oldestWaiting)

            val customer = oldestWaiting.customer
            val destClassId = oldestWaiting.destClassId
            val serverId = oldestWaiting.serverId
            val sourceClassId = oldestWaiting.sourceClassId

            val sourceStationIdx = serviceStations[sourceQueueIdx]
            logEvent("BAS_UNBLOCK", sourceStationIdx, sourceClassId, serverId, destQueueIdx)
            basUnblockCount++

            // DEBUG: Print first few BAS unblocks
            // if (basUnblockCount <= 3) {
            //     println("BAS_UNBLOCK #$basUnblockCount at T=${Sim.time()}")
            //     println("  Source Q1: currentQueueLength=${currentQueueLength[sourceQueueIdx][sourceClassId]}")
            //     println("  Dest Q2: currentQueueLength=${currentQueueLength[destQueueIdx][destClassId]}, basBlockedAtDest=${basBlockedAtDest[destQueueIdx][destClassId]}")
            // }

            // BAS: waiting job was counted at DESTINATION, so decrement blocked count
            // Update destination queue stats before decrementing (to properly time-weight)
            updateQueueStats(destQueueIdx, destClassId)
            basBlockedAtDest[destQueueIdx][destClassId]--

            // Unblock server at source
            serverBlocked[sourceQueueIdx][serverId] = false

            // Create customer for destination
            // Use the original buffer arrival time so response time includes waiting in BAS buffer
            // (job is conceptually at destination from qlen perspective, so waiting counts)
            val nextCustomer = Customer(
                destClassId,
                classPrio[destClassId],
                customer.systemArrivalTime,
                oldestWaiting.arrivalTime,  // When job arrived at BAS buffer (counts toward dest response time)
                siroRng.nextDouble(),
                jobId = customer.jobId,
                absoluteDeadline = customer.absoluteDeadline
            )

            // Send job to destination queue
            arriveAtQueue(destQueueIdx, nextCustomer, skipCapacityCheck = true)

            // Free server and start next customer at source queue (with heterogeneous tracking)
            markServerIdle(sourceQueueIdx, serverId)
            customersInService[sourceQueueIdx]--

            // Start service for next customer in source queue (if any) - heterogeneous-aware
            if (waitQueues[sourceQueueIdx].isNotEmpty()) {
                val nextWaiting = waitQueues[sourceQueueIdx].poll()
                val nextClassId = nextWaiting.classId

                // Find compatible server for next customer
                val serverSelection = findFreeServerForClass(sourceQueueIdx, nextClassId)
                if (serverSelection.serverId >= 0) {
                    markServerBusy(sourceQueueIdx, serverSelection.serverId, serverSelection.serverTypeId)
                    customersInService[sourceQueueIdx]++
                    nextWaiting.assignedServerType = serverSelection.serverTypeId

                    // Update busy time for the new customer's class
                    lastBusyUpdateTime[sourceQueueIdx][nextClassId] = Sim.time()
                    currentBusyServers[sourceQueueIdx][nextClassId]++

                    val serviceTime = if (nextWaiting.serviceTime > 0) {
                        nextWaiting.serviceTime
                    } else {
                        generateHeteroServiceTime(sourceQueueIdx, nextClassId, serverSelection.serverTypeId)
                    }
                    val departureEvent = Departure(sourceQueueIdx, serverSelection.serverId, nextWaiting)
                    departureEvent.schedule(serviceTime)

                    // Track for signal-based removal (G-networks)
                    if (hasNegativeSignals) {
                        inServiceJobs[Pair(sourceQueueIdx, serverSelection.serverId)] = InServiceJob(nextWaiting, departureEvent)
                    }

                    logEvent("BAS_NEXT_SERVICE", sourceStationIdx, nextClassId,
                        currentQueueLength[sourceQueueIdx][nextClassId], currentBusyServers[sourceQueueIdx][nextClassId])
                } else {
                    // No compatible server - put customer back
                    waitQueues[sourceQueueIdx].add(nextWaiting)
                }
            }

            // Clean up empty lists and indices
            if (outgoingBuffer.none { it.destQueueIdx == destQueueIdx }) {
                sourcesWithWaiting.remove(sourceQueueIdx)
            }
            if (outgoingBuffer.isEmpty()) {
                basOutgoingBuffer.remove(sourceQueueIdx)
            }
        }

        // Clean up reverse index if no more sources with waiting jobs for this destination
        if (sourcesWithWaiting.isEmpty()) {
            basDestinationToSources.remove(destQueueIdx)
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
    private fun getEffectivePSServerCount(queueIdx: Int): Double {
        if (isLoadDependent[queueIdx]) {
            val totalJobs = getTotalCustomersAtStation(queueIdx)
            if (totalJobs > 0) {
                val scalingArray = lldScaling[queueIdx]
                if (scalingArray != null) {
                    // lldscaling[station, n-1] gives effective server count when n jobs are present
                    val scalingIdx = minOf(totalJobs - 1, scalingArray.size - 1)
                    return scalingArray[scalingIdx]
                }
            }
            // When 0 jobs, use 1 as effective server count (first job gets full service)
            return 1.0
        }
        return numServers[queueIdx].toDouble()
    }

    /**
     * Returns the total number of customers in a finite capacity region (all classes).
     */
    private fun getTotalCustomersInRegion(regionIdx: Int): Int {
        var total = 0
        for (k in 0 until numClasses) {
            total += currentJobsInRegion[regionIdx][k]
        }
        return total
    }

    /**
     * Returns the number of customers of a specific class in a finite capacity region.
     */
    private fun getClassJobsInRegion(regionIdx: Int, classId: Int): Int {
        return if (regionIdx >= 0 && regionIdx < numRegions) {
            currentJobsInRegion[regionIdx][classId]
        } else {
            0
        }
    }

    /**
     * Updates time-weighted statistics for a finite capacity region.
     * Should be called before changing currentJobsInRegion or blockedInRegion.
     * JMT semantics: FCR QLen includes only jobs inside the region (blocked jobs are NOT included).
     * In JMT, BlockingRegion.increaseOccupation() is only called when a job successfully enters.
     */
    private fun updateRegionTimeWeightedStats(regionIdx: Int) {
        if (regionIdx < 0 || regionIdx >= numRegions) return

        val currentTime = Sim.time()
        val elapsed = currentTime - lastRegionUpdateTime[regionIdx]
        if (elapsed > 0 && warmupDone) {
            for (k in 0 until numClasses) {
                // JMT semantics: FCR QLen = weighted occupation (jobs * classWeight), blocked jobs NOT included
                // In JMT, actualOccupation += classWeights[c] when job enters region
                val weightedJobsForClass = currentJobsInRegion[regionIdx][k].toDouble() * fcRegionClassWeights[regionIdx][k]
                totalRegionJobTime[regionIdx][k] += weightedJobsForClass * elapsed
            }
        }
        lastRegionUpdateTime[regionIdx] = currentTime
    }

    /**
     * Updates FCR arrival rate tracking for a region.
     * Should be called when a job successfully enters a region (after updating currentJobsInRegion).
     * Tracks inter-arrival times for computing arrival rate = 1 / mean(inter-arrival time).
     */
    private fun updateRegionArrivalTracking(regionIdx: Int, classId: Int) {
        if (regionIdx < 0 || regionIdx >= numRegions || !warmupDone) return

        val currentTime = Sim.time()
        if (regionArrivalCount[regionIdx][classId] > 0) {
            // Not the first arrival - record inter-arrival time
            val interArrivalTime = currentTime - lastRegionArrivalTime[regionIdx][classId]
            regionInterArrivalTimeSum[regionIdx][classId] += interArrivalTime
        }
        regionArrivalCount[regionIdx][classId]++
        lastRegionArrivalTime[regionIdx][classId] = currentTime
    }

    /**
     * Attempts to release blocked customers from an FCR's waiting queue after a departure.
     * Called when a customer departs from a station in a finite capacity region.
     * Releases customers in FIFO order as long as capacity permits.
     */
    private fun tryReleaseBlockedCustomers(regionIdx: Int) {
        if (regionIdx < 0 || regionIdx >= numRegions) return

        val blockedQueue = fcRegionBlockedQueue[regionIdx]
        while (blockedQueue.isNotEmpty()) {
            val blockedCustomer = blockedQueue.peek()
            val customer = blockedCustomer.customer
            val destQueueIdx = blockedCustomer.destQueueIdx
            val classId = customer.classId

            // Check if we can admit this customer now
            val currentRegionTotal = getTotalCustomersInRegion(regionIdx)
            val globalMax = fcRegionGlobalMax.getOrNull(regionIdx) ?: Int.MAX_VALUE

            // Check global region capacity
            if (currentRegionTotal >= globalMax) {
                // Still full, can't release more
                break
            }

            // Check per-class region capacity
            val regionClassMax = fcRegionClassMax.get(regionIdx, classId).toInt()
            if (regionClassMax < Int.MAX_VALUE) {
                val currentRegionClass = getClassJobsInRegion(regionIdx, classId)
                if (currentRegionClass >= regionClassMax) {
                    // Per-class still full for this customer's class
                    // Skip this customer but continue checking others (different class might be releasable)
                    // To maintain FIFO ordering properly, we stop here - blocked customers must be released in order
                    break
                }
            }

            // Remove from blocked queue and route to destination
            blockedQueue.poll()

            // Update region stats before decrementing blocked count
            updateRegionTimeWeightedStats(regionIdx)
            blockedInRegion[regionIdx][classId]--

            // Update queue stats before decrementing FCR blocked count
            updateQueueStats(destQueueIdx, classId)
            fcrBlockedAtDest[destQueueIdx][classId]--

            // Preserve original queue arrival time (JMT semantics: blocking time is included in response time)
            val releasedCustomer = Customer(
                classId,
                customer.priority,
                customer.systemArrivalTime,
                customer.queueArrivalTime,  // Keep original arrival time for response time calculation
                customer.randomRank
            )

            logEvent("UNBLOCK_REGION", serviceStations[destQueueIdx], classId,
                currentQueueLength[destQueueIdx][classId], currentBusyServers[destQueueIdx][classId])

            // Now route the released customer to its destination queue
            // This will do all the normal processing (including updating region counts)
            processReleasedCustomerArrival(destQueueIdx, releasedCustomer)
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
    private fun processReleasedCustomerArrival(queueIdx: Int, customer: Customer) {
        val classId = customer.classId
        val stationIdx = serviceStations[queueIdx]

        // Dispatch to PS handler for Processor Sharing scheduling
        // Pass fromBlocked=true to skip queue length increment
        if (isPSScheduling(schedStrategies[queueIdx])) {
            arriveAtPSQueueFromBlocked(queueIdx, customer)
            return
        }

        // Dispatch to preemptive LCFS handler
        // Pass fromBlocked=true to skip queue length increment
        if (isPreemptiveScheduling[queueIdx]) {
            arriveAtPreemptiveLCFSQueueFromBlocked(queueIdx, customer)
            return
        }

        // Dispatch to polling handler
        // Pass fromBlocked=true to skip queue length increment
        if (isPollingStation[queueIdx]) {
            arriveAtPollingQueueFromBlocked(queueIdx, customer)
            return
        }

        // For SJF/LJF, we must generate service time upon arrival to sort the queue
        val strategy = schedStrategies[queueIdx]
        if (strategy == SchedStrategy.SJF || strategy == SchedStrategy.LJF) {
            customer.serviceTime = generateServiceTime(queueIdx, classId)
        }

        // JMT semantics: blocked jobs are NOT counted in queue QLen, only in FCR QLen
        // Now that customer is entering the queue, increment queue length
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]++

        // Track max queue length reached
        val totalAtStation = getTotalCustomersAtStation(queueIdx)
        if (totalAtStation > maxQueueLengthReached) {
            maxQueueLengthReached = totalAtStation
        }

        // Update region job counts NOW that customer is actually entering the region
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]++
            updateRegionArrivalTracking(regionIdx, classId)
        }

        if (isDelayNode[queueIdx]) {
            // Delay node (infinite server): always start service immediately
            customersInService[queueIdx]++

            // Update busy time statistics before incrementing
            updateBusyStats(queueIdx, classId)
            currentBusyServers[queueIdx][classId]++
            // Reset busy time tracking for this job - prevents including warmup time in first job of class
            lastBusyUpdateTime[queueIdx][classId] = Sim.time()

            val serviceTime = generateServiceTime(queueIdx, classId)
            val departureEvent = DelayDeparture(queueIdx, customer)
            departureEvent.schedule(serviceTime)
            // Track for signal-based removal (G-networks)
            if (hasNegativeSignals) {
                val jobId = nextDelayJobId++
                delayJobs[jobId] = DelayJob(queueIdx, customer, departureEvent)
            }
        } else {
            // Queue node: check if server is available (heterogeneous-aware)
            val serverSelection = findFreeServerForClass(queueIdx, classId)
            val freeServer = serverSelection.serverId
            val serverTypeId = serverSelection.serverTypeId
            if (freeServer >= 0) {
                // Start service immediately with heterogeneous tracking
                markServerBusy(queueIdx, freeServer, serverTypeId)
                customersInService[queueIdx]++
                customer.assignedServerType = serverTypeId

                // Update busy time statistics before incrementing
                updateBusyStats(queueIdx, classId)
                currentBusyServers[queueIdx][classId]++
                // Reset busy time tracking for this job - prevents including warmup time in first job of class
                lastBusyUpdateTime[queueIdx][classId] = Sim.time()

                val serviceTime = if (customer.serviceTime > 0) {
                    customer.serviceTime
                } else {
                    generateHeteroServiceTime(queueIdx, classId, serverTypeId)
                }
                val departureEvent = Departure(queueIdx, freeServer, customer)
                departureEvent.schedule(serviceTime)
                // Track for signal-based removal (G-networks)
                if (hasNegativeSignals) {
                    inServiceJobs[Pair(queueIdx, freeServer)] = InServiceJob(customer, departureEvent)
                }
            } else {
                // Join queue (priority queue sorts by priority, then FCFS)
                waitQueues[queueIdx].add(customer)
            }
        }

        // Log arrival event
        logEvent("RELEASED_ARRIVAL", stationIdx, classId,
            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
    }

    /**
     * Arrival at PS queue from blocked state.
     * JMT semantics: Queue length and region counts are incremented HERE when customer enters.
     */
    private fun arriveAtPSQueueFromBlocked(queueIdx: Int, customer: Customer) {
        val classId = customer.classId
        val currentTime = Sim.time()

        // Update remaining work for all current jobs based on elapsed time
        updatePSRemainingWork(queueIdx, currentTime)

        // Generate service requirement for new customer
        val serviceRequirement = generateServiceTime(queueIdx, classId)

        // Create PS customer
        val psCustomer = PSCustomer(
            classId = classId,
            priority = customer.priority,
            systemArrivalTime = customer.systemArrivalTime,
            queueArrivalTime = customer.queueArrivalTime,
            totalServiceRequirement = serviceRequirement,
            remainingServiceWork = serviceRequirement,
            scheduledDepartureEvent = null,
            forkedJob = null
        )

        // Update busy time statistics for PS BEFORE adding the new job
        updatePSBusyStats(queueIdx)

        // Add to jobs in service (in PS, all jobs are always "in service")
        psJobsInService[queueIdx].add(psCustomer)

        // JMT semantics: blocked jobs are NOT counted in queue QLen, only in FCR QLen
        // Now that customer is entering the queue, increment queue length
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]++

        // Track total customers in service
        customersInService[queueIdx]++

        // Update region job counts NOW that customer is actually entering the region
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]++
            updateRegionArrivalTracking(regionIdx, classId)
        }

        // Reschedule all departures with new rates
        rescheduleAllPSDepartures(queueIdx)

        logEvent("PS_RELEASED_ARRIVAL", serviceStations[queueIdx], classId,
            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
    }

    /**
     * Arrival at preemptive LCFS queue from blocked state.
     * JMT semantics: Queue length and region counts are incremented HERE when customer enters.
     */
    private fun arriveAtPreemptiveLCFSQueueFromBlocked(queueIdx: Int, customer: Customer) {
        val classId = customer.classId
        val strategy = schedStrategies[queueIdx]

        // JMT semantics: blocked jobs are NOT counted in queue QLen, only in FCR QLen
        // Now that customer is entering the queue, increment queue length
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]++

        // Update region job counts NOW that customer is actually entering the region
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]++
            updateRegionArrivalTracking(regionIdx, classId)
        }

        // Track max queue length reached
        val totalAtStation = getTotalCustomersAtStation(queueIdx)
        if (totalAtStation > maxQueueLengthReached) {
            maxQueueLengthReached = totalAtStation
        }

        // Check if we should preempt current job
        val shouldPreempt = shouldPreemptForLCFS(queueIdx, customer)

        if (shouldPreempt) {
            // Preempt current job
            val victimJob = findJobToPreempt(queueIdx, customer, strategy)
            if (victimJob != null) {
                preemptJob(queueIdx, victimJob)
                // Inherit victim's server type
                startPreemptiveService(queueIdx, customer, victimJob.serverId, victimJob.assignedServerType)
            } else {
                waitQueues[queueIdx].add(customer)
            }
        } else {
            // Find free server (heterogeneous-aware)
            val serverSelection = findFreeServerForClass(queueIdx, classId)
            if (serverSelection.serverId >= 0) {
                startPreemptiveService(queueIdx, customer, serverSelection.serverId, serverSelection.serverTypeId)
            } else {
                waitQueues[queueIdx].add(customer)
            }
        }

        logEvent("LCFS_RELEASED_ARRIVAL", serviceStations[queueIdx], classId,
            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])
    }

    /**
     * Arrival at polling queue from blocked state.
     * JMT semantics: Queue length and region counts are incremented HERE when customer enters.
     */
    private fun arriveAtPollingQueueFromBlocked(queueIdx: Int, customer: Customer) {
        val classId = customer.classId

        // JMT semantics: blocked jobs are NOT counted in queue QLen, only in FCR QLen
        // Now that customer is entering the queue, increment queue length
        updateQueueStats(queueIdx, classId)
        currentQueueLength[queueIdx][classId]++

        // Track max queue length reached
        val totalAtStation = getTotalCustomersAtStation(queueIdx)
        if (totalAtStation > maxQueueLengthReached) {
            maxQueueLengthReached = totalAtStation
        }

        // Update region job counts NOW that customer is actually entering the region
        val stationIdx = serviceStations[queueIdx]
        if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
            val regionIdx = fcRegionIndices[stationIdx]
            updateRegionTimeWeightedStats(regionIdx)
            currentJobsInRegion[regionIdx][classId]++
            updateRegionArrivalTracking(regionIdx, classId)
        }

        // Add to per-class polling queue
        pollingQueues[queueIdx][classId].add(customer)

        logEvent("POLL_RELEASED_ARRIVAL", serviceStations[queueIdx], classId,
            currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])

        // If server is idle and not in switchover, try to start service
        if (customersInService[queueIdx] == 0 && !pollingInSwitchover[queueIdx]) {
            pollingTryStartService(queueIdx)
        }
    }

    /**
     * Departure event when a customer finishes service at a queue.
     */
    private inner class Departure(
        private val queueIdx: Int,
        private val serverId: Int,
        private val customer: Customer
    ) : Event() {
        override fun actions() {
            val classId = customer.classId

            // Record queue response time
            val queueResponseTime = Sim.time() - customer.queueArrivalTime
            responseTimeTally[queueIdx][classId].add(queueResponseTime)
            responseTimeSamples[queueIdx][classId].add(queueResponseTime)
            completedCustomers[queueIdx][classId]++

            // Record queue tardiness (relative to deadline)
            val queueTardiness = maxOf(0.0, Sim.time() - customer.absoluteDeadline)
            tardinessTally[queueIdx][classId].add(queueTardiness)

            // Check event count for stopping/warmup/MSER sampling
            checkEventCountStop()

            // Update queue length statistics
            updateQueueStats(queueIdx, classId)
            currentQueueLength[queueIdx][classId]--

            // Update region job counts if in a region and try releasing blocked customers
            val stationIdx = serviceStations[queueIdx]
            if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
                val regionIdx = fcRegionIndices[stationIdx]
                updateRegionTimeWeightedStats(regionIdx)
                currentJobsInRegion[regionIdx][classId]--
                if (warmupDone) {
                    regionCompletions[regionIdx][classId]++
                }
                // Try to release blocked customers now that space is available
                tryReleaseBlockedCustomers(regionIdx)
            }

            // Try to admit blocked jobs now that this queue has space (must be done after queue length decremented)
            // BAS: Unblock servers at source and admit their jobs (source server was blocked)
            tryAdmitBASWaitingJobs(queueIdx)
            // BBS: Unblock servers at source and admit their jobs (source server was blocked)
            tryUnblockBBSServers(queueIdx)

            // Update busy time before decrementing
            updateBusyStats(queueIdx, classId)
            currentBusyServers[queueIdx][classId]--

            // Clean up in-service tracking (G-networks)
            if (hasNegativeSignals) {
                inServiceJobs.remove(Pair(queueIdx, serverId))
            }

            // Log departure event
            logEvent("DEPARTURE", stationIdx, classId,
                currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])

            // Check if this class expects a reply signal (synchronous call semantics)
            val expectsReply = synchCallReplyClass[classId] >= 0
            val jobId = customer.jobId

            // Route customer to next destination
            val currentNode = serviceNodes[queueIdx]
            val routingResult = selectDestination(currentNode, classId)
            val destNode = routingResult.destNode
            val destClassId = routingResult.destClassId

            if (destNode >= 0 && sinkNodes.contains(destNode)) {
                // Customer leaves system (count in original class for consistency)
                val systemResponseTime = Sim.time() - customer.systemArrivalTime
                systemResponseTimeTally[classId].add(systemResponseTime)

                // Record system tardiness (max(0, completion_time - deadline))
                val systemTardiness = maxOf(0.0, Sim.time() - customer.absoluteDeadline)
                systemTardinessTally[classId].add(systemTardiness)

                systemCompletedCustomers[classId]++
            } else if (destNode >= 0 && forkNodes.contains(destNode)) {
                // Destination is a Fork node - create forked children
                val parentJobId = nextJobId++
                handleForkArrival(destNode, parentJobId, destClassId, customer.systemArrivalTime)
            } else if (destNode >= 0 && joinNodes.contains(destNode)) {
                // Destination is a Join node - check if this is a forked job
                val forkedJobKey = Pair(queueIdx, customer.queueArrivalTime)
                val forkedJob = forkedCustomerMap.remove(forkedJobKey)
                if (forkedJob != null) {
                    // Create a new forked job with current state for join handling
                    val updatedForkedJob = ForkedJob(
                        forkJobId = forkedJob.forkJobId,
                        parentJobId = forkedJob.parentJobId,
                        classId = destClassId,
                        priority = forkedJob.priority,
                        systemArrivalTime = forkedJob.systemArrivalTime,
                        queueArrivalTime = Sim.time(),
                        randomRank = forkedJob.randomRank
                    )
                    handleJoinArrival(destNode, updatedForkedJob)
                } else {
                    // Not a tracked forked job - may be from a different path
                    // Try to find matching fork info by system arrival time
                    handleUnknownJoinArrival(destNode, destClassId, customer.systemArrivalTime)
                }
            } else if (destNode >= 0) {
                // Route to another queue (possibly with switched class)
                val nextQueueIdx = serviceNodes.indexOf(destNode)
                if (nextQueueIdx >= 0) {
                    // Check if class-switched to a negative signal or catastrophe
                    if ((hasNegativeSignals && isNegativeSignal[destClassId]) ||
                        (hasCatastropheSignals && isCatastropheSignal[destClassId])) {
                        handleNegativeSignalArrival(nextQueueIdx, destClassId, customer.systemArrivalTime)
                    } else if (hasReplySignals && isReplySignal[destClassId]) {
                        // Class-switched to a REPLY signal - route as reply with job ID
                        handleReplySignalArrival(nextQueueIdx, destClassId, jobId, customer.systemArrivalTime)
                    } else {
                        // Check blocking policies at BOTH source and destination queues
                        // LINE convention: BAS can be set on source queue (jobs block when leaving if downstream is full)
                        // JMT convention: BAS can also be set on destination queue (the queue with limited capacity)
                        // Support both: check source OR destination for blocking policy
                        val hasBBS = hasBBSBlocking(queueIdx, classId) || hasBBSBlocking(nextQueueIdx, destClassId)
                        val hasBAS = hasBASBlocking(queueIdx, classId) || hasBASBlocking(nextQueueIdx, destClassId)
                        val destHasCap = destinationHasCapacity(nextQueueIdx, destClassId)

                        if (hasBBS && !destHasCap) {
                            // BBS (Blocking Before Service): Block server at source until destination has capacity
                            // Server cannot process new jobs while blocked
                            val blockedServer = BBSBlockedServer(
                                customer = customer,
                                destQueueIdx = nextQueueIdx,
                                destClassId = destClassId,
                                sourceQueueIdx = queueIdx,
                                serverId = serverId,
                                sourceClassId = classId,
                                blockStartTime = Sim.time()
                            )

                            // Add to blocked servers list (FIFO order)
                            val blockedList = bbsBlockedServers.getOrPut(queueIdx) { mutableListOf() }
                            blockedList.add(blockedServer)

                            // Update reverse index for efficient unblocking
                            val sourcesForDest = bbsDestinationToSources.getOrPut(nextQueueIdx) { mutableSetOf() }
                            sourcesForDest.add(queueIdx)

                            // Mark server as blocked (server stays busy but blocked)
                            serverBlocked[queueIdx][serverId] = true

                            // Blocked job counts towards DESTINATION queue length
                            updateQueueStats(nextQueueIdx, destClassId)
                            bbsBlockedAtDest[nextQueueIdx][destClassId]++
                            updateQueueStats(nextQueueIdx, destClassId)

                            logEvent("BBS_BLOCK", stationIdx, classId, serverId, nextQueueIdx)

                            // Do NOT start next customer - server is blocked
                            return
                        } else if (hasBAS && !destHasCap) {
                            // BAS (Blocking After Service): Block server at source until destination has capacity
                            // Job waits at source occupying the server space; server cannot process new jobs
                            val waitingJob = BASWaitingJob(
                                customer = customer,
                                destQueueIdx = nextQueueIdx,
                                destClassId = destClassId,
                                sourceQueueIdx = queueIdx,
                                serverId = serverId,
                                sourceClassId = classId,
                                arrivalTime = Sim.time()
                            )

                            // Add to source's outgoing buffer (FIFO order)
                            val outgoingBuffer = basOutgoingBuffer.getOrPut(queueIdx) { mutableListOf() }
                            outgoingBuffer.add(waitingJob)

                            // Update reverse index for efficient admission
                            val sourcesForDest = basDestinationToSources.getOrPut(nextQueueIdx) { mutableSetOf() }
                            sourcesForDest.add(queueIdx)

                            // Mark server as blocked (server stays busy but blocked)
                            serverBlocked[queueIdx][serverId] = true

                            // Job counts towards DESTINATION queue length (even though physically at source)
                            // Increment blocked count BEFORE updateQueueStats so it's included in time-weighting
                            basBlockedAtDest[nextQueueIdx][destClassId]++
                            updateQueueStats(nextQueueIdx, destClassId)

                            logEvent("BAS_BLOCK", stationIdx, classId, serverId, nextQueueIdx)
                            basBlockCount++

                            // DEBUG: Print first few BAS blocks
                            // if (basBlockCount <= 3) {
                            //     println("BAS_BLOCK #$basBlockCount at T=${Sim.time()}")
                            //     println("  Q1 currentQueueLength=${currentQueueLength[queueIdx][classId]}")
                            //     println("  Q2 currentQueueLength=${currentQueueLength[nextQueueIdx][destClassId]}, basBlockedAtDest=${basBlockedAtDest[nextQueueIdx][destClassId]}")
                            //     println("  Q1 totalQueueTime=${totalQueueTime[queueIdx][classId]}")
                            //     println("  Q2 totalQueueTime=${totalQueueTime[nextQueueIdx][destClassId]}")
                            // }

                            // Do NOT start next customer - server is blocked
                            return
                        } else {
                            // Check for immediate feedback (self-loop staying in service)
                            val isImmediateFeedback = nextQueueIdx == queueIdx &&
                                sn.immfeed != null &&
                                stationIdx >= 0 && stationIdx < sn.immfeed.numRows &&
                                destClassId >= 0 && destClassId < sn.immfeed.numCols &&
                                sn.immfeed.get(stationIdx, destClassId) > 0.0

                            if (isImmediateFeedback) {
                                // Immediate feedback: job stays in service on self-loop
                                // Generate new service time for destination class and schedule departure
                                // Do NOT free server, do NOT call arriveAtQueue

                                // Update busy time tracking for class switch
                                if (classId != destClassId) {
                                    updateBusyStats(queueIdx, classId)
                                    currentBusyServers[queueIdx][classId]--
                                    lastBusyUpdateTime[queueIdx][destClassId] = Sim.time()
                                    currentBusyServers[queueIdx][destClassId]++
                                }

                                // Create new customer with updated class for stats tracking
                                val nextCustomer = Customer(
                                    destClassId, classPrio[destClassId],
                                    customer.systemArrivalTime, Sim.time(),
                                    siroRng.nextDouble(), jobId = jobId,
                                    absoluteDeadline = customer.absoluteDeadline
                                )

                                // Generate service time and schedule departure on same server
                                val serviceTime = generateHeteroServiceTime(queueIdx, destClassId, customer.assignedServerType)
                                val departureEvent = Departure(queueIdx, serverId, nextCustomer)
                                departureEvent.schedule(serviceTime)

                                // Track for signal-based removal (G-networks)
                                if (hasNegativeSignals) {
                                    inServiceJobs[Pair(queueIdx, serverId)] = InServiceJob(nextCustomer, departureEvent)
                                }

                                logEvent("IMMFEED", stationIdx, destClassId,
                                    currentQueueLength[queueIdx][destClassId], currentBusyServers[queueIdx][destClassId])

                                // Do NOT free server or start next customer - job continues in service
                                return
                            }

                            // Normal routing: destination has capacity or no blocking policy
                            // Check if this customer was a forked job and carry forward tracking
                            val forkedJobKey = Pair(queueIdx, customer.queueArrivalTime)
                            val forkedJob = forkedCustomerMap.remove(forkedJobKey)
                            // Create next customer with same job ID for reply tracking
                            val nextCustomer = Customer(destClassId, classPrio[destClassId], customer.systemArrivalTime, Sim.time(), siroRng.nextDouble(), jobId = jobId,
                                absoluteDeadline = customer.absoluteDeadline)
                            if (forkedJob != null) {
                                // Carry forward forked job tracking
                                val nextForkedJob = ForkedJob(
                                    forkJobId = forkedJob.forkJobId,
                                    parentJobId = forkedJob.parentJobId,
                                    classId = destClassId,
                                    priority = forkedJob.priority,
                                    systemArrivalTime = forkedJob.systemArrivalTime,
                                    queueArrivalTime = Sim.time(),
                                    randomRank = forkedJob.randomRank
                                )
                                arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob)
                            } else {
                                arriveAtQueue(nextQueueIdx, nextCustomer)
                            }
                        }
                    }
                }
            }

            // Handle synchronous call blocking: if this class expects a reply, block server
            // Only block if this is the first server to handle this job ID (not downstream servers).
            // If pendingReplyMap already has this job ID, we're at a downstream server and should not block.
            // Also, if we're class-switching TO a REPLY signal, don't block - this is the server producing the reply.
            val classSwitchedToReply = hasReplySignals && destClassId >= 0 && isReplySignal[destClassId]
            if (expectsReply && jobId >= 0 && !pendingReplyMap.containsKey(jobId) && !classSwitchedToReply) {
                // Record pending reply expectation
                val pendingReply = PendingReply(
                    jobId = jobId,
                    originalClassId = classId,
                    queueIdx = queueIdx,
                    serverId = serverId,
                    blockStartTime = Sim.time()
                )
                pendingReplyMap[jobId] = pendingReply

                // Mark server as blocked (server stays busy but blocked)
                // Note: serverBusy is already true from when the job was being served.
                // We don't decrement customersInService because the server remains occupied (blocked).
                serverBlocked[queueIdx][serverId] = true
                // serverBusy[queueIdx][serverId] is already true

                // Update queue stats before incrementing blocked count (to properly time-weight)
                updateQueueStats(queueIdx, classId)
                // Track blocked server count for this class (for utilization and queue length with blocking)
                currentBlockedServers[queueIdx][classId]++

                logEvent("SYNCH_BLOCK", stationIdx, classId, serverId, 0)

                // Do NOT start next customer - server is blocked waiting for reply
                return
            }

            // Free server (normal case - no reply expected, with heterogeneous tracking)
            markServerIdle(queueIdx, serverId)
            customersInService[queueIdx]--

            // Start service for next customer in queue (if any)
            // Priority queue will return highest priority customer first
            if (waitQueues[queueIdx].isNotEmpty()) {
                val nextCustomer = waitQueues[queueIdx].poll()
                val nextClassId = nextCustomer.classId

                // Find compatible server for next customer (heterogeneous-aware)
                val serverSelection = findFreeServerForClass(queueIdx, nextClassId)
                if (serverSelection.serverId >= 0) {
                    markServerBusy(queueIdx, serverSelection.serverId, serverSelection.serverTypeId)
                    customersInService[queueIdx]++
                    nextCustomer.assignedServerType = serverSelection.serverTypeId

                    // Reset busy time tracking for the new customer's class to prevent phantom busy time accumulation
                    // from the interval before this job started service. Do NOT call updateBusyStats here since the
                    // server was already busy (serving a different class) - we only start tracking THIS class's busy time.
                    lastBusyUpdateTime[queueIdx][nextClassId] = Sim.time()
                    currentBusyServers[queueIdx][nextClassId]++

                    val serviceTime = if (nextCustomer.serviceTime > 0) {
                        nextCustomer.serviceTime
                    } else {
                        generateHeteroServiceTime(queueIdx, nextClassId, serverSelection.serverTypeId)
                    }
                    val departureEvent = Departure(queueIdx, serverSelection.serverId, nextCustomer)
                    departureEvent.schedule(serviceTime)
                    // Track for signal-based removal (G-networks)
                    if (hasNegativeSignals) {
                        inServiceJobs[Pair(queueIdx, serverSelection.serverId)] = InServiceJob(nextCustomer, departureEvent)
                    }
                } else {
                    // No compatible server available - put customer back in queue
                    waitQueues[queueIdx].add(nextCustomer)
                }
            } else {
                // Queue is empty - initiate delayoff if enabled
                if (isSetupDelayoffEnabled(queueIdx, classId)) {
                    startServerDelayoff(queueIdx, serverId, classId)
                }
            }
        }
    }

    /**
     * Departure event for preemptive LCFS scheduling.
     * Handles job completion and starts service for next waiting customer.
     */
    private inner class PreemptiveDeparture(
        private val queueIdx: Int,
        private val serverId: Int,
        private val customer: PreemptiveCustomer
    ) : Event() {
        override fun actions() {
            val classId = customer.classId
            val currentTime = Sim.time()

            // Remove from in-service tracking - if not present, job was preempted
            val wasInService = preemptiveJobsInService[queueIdx].remove(customer)
            if (!wasInService) {
                // Job was already preempted, this is a stale departure event
                return
            }

            // Record queue response time
            val queueResponseTime = currentTime - customer.queueArrivalTime
            if (warmupDone) {
                responseTimeTally[queueIdx][classId].add(queueResponseTime)
                responseTimeSamples[queueIdx][classId].add(queueResponseTime)
                completedCustomers[queueIdx][classId]++
            }

            // Check event count for stopping/warmup/MSER sampling (always count, regardless of warmup)
            checkEventCountStop()

            // Update queue length statistics
            updateQueueStats(queueIdx, classId)
            currentQueueLength[queueIdx][classId]--

            // Update region job counts if in a region and try releasing blocked customers
            val stationIdx = serviceStations[queueIdx]
            if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
                val regionIdx = fcRegionIndices[stationIdx]
                updateRegionTimeWeightedStats(regionIdx)
                currentJobsInRegion[regionIdx][classId]--
                if (warmupDone) {
                    regionCompletions[regionIdx][classId]++
                }
                // Try to release blocked customers now that space is available
                tryReleaseBlockedCustomers(regionIdx)
            }

            // Update busy time before decrementing
            updateBusyStats(queueIdx, classId)
            currentBusyServers[queueIdx][classId]--

            // Free server (with heterogeneous tracking)
            markServerIdle(queueIdx, serverId)
            customersInService[queueIdx]--

            // Log departure event
            logEvent("DEPARTURE", stationIdx, classId,
                currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])

            // Route customer to next destination
            val currentNode = serviceNodes[queueIdx]
            val routingResult = selectDestination(currentNode, classId)
            val destNode = routingResult.destNode
            val destClassId = routingResult.destClassId

            if (destNode >= 0 && sinkNodes.contains(destNode)) {
                // Customer leaves system (count in original class for consistency)
                val systemResponseTime = currentTime - customer.systemArrivalTime
                if (warmupDone) {
                    systemResponseTimeTally[classId].add(systemResponseTime)
                    systemCompletedCustomers[classId]++
                }
            } else if (destNode >= 0 && forkNodes.contains(destNode)) {
                // Destination is a Fork node - create forked children
                val parentJobId = nextJobId++
                handleForkArrival(destNode, parentJobId, destClassId, customer.systemArrivalTime)
            } else if (destNode >= 0 && joinNodes.contains(destNode)) {
                // Destination is a Join node
                val forkedJobKey = Pair(queueIdx, customer.queueArrivalTime)
                val forkedJob = forkedCustomerMap.remove(forkedJobKey)
                if (forkedJob != null) {
                    val updatedForkedJob = ForkedJob(
                        forkJobId = forkedJob.forkJobId,
                        parentJobId = forkedJob.parentJobId,
                        classId = destClassId,
                        priority = forkedJob.priority,
                        systemArrivalTime = forkedJob.systemArrivalTime,
                        queueArrivalTime = currentTime,
                        randomRank = forkedJob.randomRank
                    )
                    handleJoinArrival(destNode, updatedForkedJob)
                } else {
                    handleUnknownJoinArrival(destNode, destClassId, customer.systemArrivalTime)
                }
            } else if (destNode >= 0) {
                // Route to another queue (possibly with switched class)
                val nextQueueIdx = serviceNodes.indexOf(destNode)
                if (nextQueueIdx >= 0) {
                    val forkedJobKey = Pair(queueIdx, customer.queueArrivalTime)
                    val forkedJob = forkedCustomerMap.remove(forkedJobKey)
                    val nextCustomer = Customer(
                        destClassId, classPrio[destClassId], customer.systemArrivalTime,
                        currentTime, customer.randomRank
                    )
                    if (forkedJob != null) {
                        val nextForkedJob = ForkedJob(
                            forkJobId = forkedJob.forkJobId,
                            parentJobId = forkedJob.parentJobId,
                            classId = destClassId,
                            priority = forkedJob.priority,
                            systemArrivalTime = forkedJob.systemArrivalTime,
                            queueArrivalTime = currentTime,
                            randomRank = forkedJob.randomRank
                        )
                        arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob)
                    } else {
                        arriveAtQueue(nextQueueIdx, nextCustomer)
                    }
                }
            }

            // Start service for next customer in queue (if any) - heterogeneous-aware
            if (waitQueues[queueIdx].isNotEmpty()) {
                val nextCustomer = waitQueues[queueIdx].poll()
                val nextClassId = nextCustomer.classId
                val serverSelection = findFreeServerForClass(queueIdx, nextClassId)
                if (serverSelection.serverId >= 0) {
                    startPreemptiveService(queueIdx, nextCustomer, serverSelection.serverId, serverSelection.serverTypeId)
                } else {
                    // No compatible server available - put customer back in queue
                    waitQueues[queueIdx].add(nextCustomer)
                }
            }
        }
    }

    /**
     * Setup completion event when a server finishes its cold start/initialization phase.
     * The server transitions from SETUP to ACTIVE state and begins serving the waiting job.
     */
    private inner class SetupCompletion(
        private val queueIdx: Int,
        private val serverId: Int,
        private val classId: Int
    ) : Event() {
        override fun actions() {
            // Update setup statistics before changing state
            updateSetupStats(queueIdx, classId)
            currentServersInSetup[queueIdx][classId]--

            // Transition server to ACTIVE state
            serverState[queueIdx][serverId] = ServerState.ACTIVE

            // Log setup completion
            val stationIdx = serviceStations[queueIdx]
            logEvent("SETUP_COMPLETE", stationIdx, classId, serverId, 0)

            // Check if there are jobs waiting in queue for this server
            if (waitQueues[queueIdx].isNotEmpty()) {
                // Find a waiting job (prefer same class that triggered setup)
                val nextCustomer = waitQueues[queueIdx].poll()

                // Start service for this job
                startService(queueIdx, serverId, nextCustomer)
            }
            // If no jobs waiting, server remains ACTIVE but idle
        }
    }

    /**
     * Delayoff completion event when a server finishes its teardown phase.
     * The server transitions from DELAYOFF to OFF state.
     */
    private inner class DelayoffCompletion(
        private val queueIdx: Int,
        private val serverId: Int
    ) : Event() {
        override fun actions() {
            // Verify server is still in DELAYOFF state (could have been interrupted)
            if (serverState[queueIdx][serverId] != ServerState.DELAYOFF) {
                return  // Server was resumed, delayoff was cancelled
            }

            // Get the class that initiated delayoff
            val classId = serverLastClass[queueIdx][serverId]
            if (classId >= 0) {
                // Update delayoff statistics before changing state
                updateDelayoffStats(queueIdx, classId)
                currentServersInDelayoff[queueIdx][classId]--
            }

            // Transition server to OFF state
            serverState[queueIdx][serverId] = ServerState.OFF

            // Clear pending event reference
            pendingDelayoffEvents[queueIdx][serverId] = null

            // Log delayoff completion
            val stationIdx = serviceStations[queueIdx]
            logEvent("DELAYOFF_COMPLETE", stationIdx, classId, serverId, 0)
        }
    }

    /**
     * Departure event when a customer finishes service at a Delay node (infinite server).
     * No queue management needed since all customers are served immediately.
     */
    private inner class DelayDeparture(
        private val queueIdx: Int,
        private val customer: Customer
    ) : Event() {
        override fun actions() {
            val classId = customer.classId

            // Record response time (equals service time for Delay nodes)
            val responseTime = Sim.time() - customer.queueArrivalTime
            responseTimeTally[queueIdx][classId].add(responseTime)
            responseTimeSamples[queueIdx][classId].add(responseTime)
            completedCustomers[queueIdx][classId]++

            // Check event count for stopping/warmup/MSER sampling
            checkEventCountStop()

            // Update queue length statistics
            updateQueueStats(queueIdx, classId)
            currentQueueLength[queueIdx][classId]--

            // Update region job counts if in a region and try releasing blocked customers
            val stationIdx = serviceStations[queueIdx]
            if (stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
                val regionIdx = fcRegionIndices[stationIdx]
                updateRegionTimeWeightedStats(regionIdx)
                currentJobsInRegion[regionIdx][classId]--
                if (warmupDone) {
                    regionCompletions[regionIdx][classId]++
                }
                // Try to release blocked customers now that space is available
                tryReleaseBlockedCustomers(regionIdx)
            }

            // Update busy time before decrementing
            updateBusyStats(queueIdx, classId)
            currentBusyServers[queueIdx][classId]--

            // Decrement customers in service
            customersInService[queueIdx]--

            // Clean up delay job tracking (G-networks)
            if (hasNegativeSignals) {
                // Find and remove this customer from delayJobs by matching customer reference
                val jobIdToRemove = delayJobs.entries.find {
                    it.value.queueIdx == queueIdx && it.value.customer === customer
                }?.key
                if (jobIdToRemove != null) {
                    delayJobs.remove(jobIdToRemove)
                }
            }

            // Log delay departure event
            logEvent("DELAY_DEPARTURE", stationIdx, classId,
                currentQueueLength[queueIdx][classId], currentBusyServers[queueIdx][classId])

            // Route customer to next destination
            val currentNode = serviceNodes[queueIdx]
            val routingResult = selectDestination(currentNode, classId)
            val destNode = routingResult.destNode
            val destClassId = routingResult.destClassId

            if (destNode >= 0 && sinkNodes.contains(destNode)) {
                // Customer leaves system (count in original class for consistency)
                val systemResponseTime = Sim.time() - customer.systemArrivalTime
                systemResponseTimeTally[classId].add(systemResponseTime)

                // Record system tardiness (max(0, completion_time - deadline))
                val systemTardiness = maxOf(0.0, Sim.time() - customer.absoluteDeadline)
                systemTardinessTally[classId].add(systemTardiness)

                systemCompletedCustomers[classId]++
            } else if (destNode >= 0 && forkNodes.contains(destNode)) {
                // Destination is a Fork node - create forked children
                val parentJobId = nextJobId++
                handleForkArrival(destNode, parentJobId, destClassId, customer.systemArrivalTime)
            } else if (destNode >= 0 && joinNodes.contains(destNode)) {
                // Destination is a Join node - check if this is a forked job
                val forkedJobKey = Pair(queueIdx, customer.queueArrivalTime)
                val forkedJob = forkedCustomerMap.remove(forkedJobKey)
                if (forkedJob != null) {
                    val updatedForkedJob = ForkedJob(
                        forkJobId = forkedJob.forkJobId,
                        parentJobId = forkedJob.parentJobId,
                        classId = destClassId,
                        priority = forkedJob.priority,
                        systemArrivalTime = forkedJob.systemArrivalTime,
                        queueArrivalTime = Sim.time(),
                        randomRank = forkedJob.randomRank
                    )
                    handleJoinArrival(destNode, updatedForkedJob)
                } else {
                    handleUnknownJoinArrival(destNode, destClassId, customer.systemArrivalTime)
                }
            } else if (destNode >= 0) {
                // Route to another queue/delay (possibly with switched class)
                val nextQueueIdx = serviceNodes.indexOf(destNode)
                if (nextQueueIdx >= 0) {
                    // Check if class-switched to a negative signal or catastrophe
                    if ((hasNegativeSignals && isNegativeSignal[destClassId]) ||
                        (hasCatastropheSignals && isCatastropheSignal[destClassId])) {
                        handleNegativeSignalArrival(nextQueueIdx, destClassId, customer.systemArrivalTime)
                    } else if (hasReplySignals && isReplySignal[destClassId]) {
                        // Class-switched to a REPLY signal - route as reply with job ID
                        handleReplySignalArrival(nextQueueIdx, destClassId, customer.jobId, customer.systemArrivalTime)
                    } else {
                        val forkedJobKey = Pair(queueIdx, customer.queueArrivalTime)
                        val forkedJob = forkedCustomerMap.remove(forkedJobKey)
                        // Create next customer with same job ID for reply tracking
                        val nextCustomer = Customer(destClassId, classPrio[destClassId], customer.systemArrivalTime, Sim.time(), siroRng.nextDouble(), jobId = customer.jobId,
                            absoluteDeadline = customer.absoluteDeadline)
                        if (forkedJob != null) {
                            val nextForkedJob = ForkedJob(
                                forkJobId = forkedJob.forkJobId,
                                parentJobId = forkedJob.parentJobId,
                                classId = destClassId,
                                priority = forkedJob.priority,
                                systemArrivalTime = forkedJob.systemArrivalTime,
                                queueArrivalTime = Sim.time(),
                                randomRank = forkedJob.randomRank
                            )
                            arriveAtQueueForked(nextQueueIdx, nextCustomer, nextForkedJob)
                        } else {
                            arriveAtQueue(nextQueueIdx, nextCustomer)
                        }
                    }
                }
            }
        }
    }

    /**
     * End of warmup event - resets statistics and schedules end of simulation.
     */
    private inner class EndOfWarmup(private val remainingSimTime: Double) : Event() {
        override fun actions() {
            // Reset all statistics after warmup
            resetStatistics()
            warmupDone = true

            // Schedule end of simulation
            EndOfSimulation().schedule(remainingSimTime)
        }
    }

    /**
     * Resets all statistics probes after warmup period.
     */
    private fun resetStatistics() {
        warmupEndTime = Sim.time()
        for (qIdx in 0 until numServiceNodes) {
            for (k in 0 until numClasses) {
                responseTimeTally[qIdx][k].init()
                completedCustomers[qIdx][k] = 0
                totalQueueTime[qIdx][k] = 0.0
                lastQueueUpdateTime[qIdx][k] = Sim.time()
                // Reset busy time tracking
                totalBusyTime[qIdx][k] = 0.0
                lastBusyUpdateTime[qIdx][k] = Sim.time()
                // Reset blocking time tracking
                totalBlockingTime[qIdx][k] = 0.0
                // Note: currentQueueLength, currentBusyServers, currentBlockedServers are NOT reset - they reflect actual state
            }
        }
        for (k in 0 until numClasses) {
            systemResponseTimeTally[k].init()
            systemCompletedCustomers[k] = 0
        }
        // Reset Join node statistics
        for (joinListIdx in joinNodes.indices) {
            for (k in 0 until numClasses) {
                joinResponseTimeTally[joinListIdx][k].init()
                joinCompletions[joinListIdx][k] = 0
                totalJoinQueueTime[joinListIdx][k] = 0.0
                lastJoinUpdateTime[joinListIdx][k] = Sim.time()
                arrivedAtJoin[joinListIdx][k] = 0
                // Note: currentJoinQueueLength is NOT reset - it reflects actual state
            }
        }
    }

    /**
     * Updates busy time statistics for a queue and class.
     */
    private fun updateBusyStats(queueIdx: Int, classId: Int) {
        val currentTime = Sim.time()
        val elapsed = currentTime - lastBusyUpdateTime[queueIdx][classId]
        if (elapsed > 0) {
            totalBusyTime[queueIdx][classId] += currentBusyServers[queueIdx][classId] * elapsed
            lastBusyUpdateTime[queueIdx][classId] = currentTime
        }
    }

    /**
     * Updates setup time statistics for a service node and class.
     * Tracks time-weighted average of servers in setup phase.
     * Setup time is NOT counted as busy time (it's initialization overhead).
     */
    private fun updateSetupStats(queueIdx: Int, classId: Int) {
        if (!warmupDone) return  // Don't accumulate statistics during warmup

        val currentTime = Sim.time()
        val elapsed = currentTime - lastSetupUpdateTime[queueIdx][classId]
        if (elapsed > 0) {
            totalSetupTime[queueIdx][classId] += currentServersInSetup[queueIdx][classId] * elapsed
            lastSetupUpdateTime[queueIdx][classId] = currentTime
        }
    }

    /**
     * Updates delayoff time statistics for a service node and class.
     * Tracks time-weighted average of servers in delayoff (teardown) phase.
     * Delayoff time is NOT counted as busy time (it's teardown overhead).
     */
    private fun updateDelayoffStats(queueIdx: Int, classId: Int) {
        if (!warmupDone) return  // Don't accumulate statistics during warmup

        val currentTime = Sim.time()
        val elapsed = currentTime - lastDelayoffUpdateTime[queueIdx][classId]
        if (elapsed > 0) {
            totalDelayoffTime[queueIdx][classId] += currentServersInDelayoff[queueIdx][classId] * elapsed
            lastDelayoffUpdateTime[queueIdx][classId] = currentTime
        }
    }

    /**
     * Updates busy time statistics for a PS queue.
     * Uses the actual allocated rates to determine utilization per class.
     */
    private fun updatePSBusyStats(queueIdx: Int) {
        val currentTime = Sim.time()
        val elapsed = currentTime - psLastBusyUpdateTime[queueIdx]
        if (elapsed <= 0) {
            psLastBusyUpdateTime[queueIdx] = currentTime
            return
        }

        val jobs = psJobsInService[queueIdx]
        if (jobs.isNotEmpty()) {
            val strategy = schedStrategies[queueIdx]
            // For load-dependent PS queues, use the effective server count based on current population
            val c = getEffectivePSServerCount(queueIdx)

            // Calculate effective rates (speed factors) for each job
            // Note: these rates sum up to at most c (total server capacity)
            val rates = calculatePSRates(queueIdx, strategy, jobs, c)

            // Aggregate utilization by class
            for ((idx, job) in jobs.withIndex()) {
                val rate = rates[idx]
                if (rate > 0) {
                    totalBusyTime[queueIdx][job.classId] += rate * elapsed
                }
            }
        }

        psLastBusyUpdateTime[queueIdx] = currentTime
    }

    /**
     * Returns the actual simulation time (excluding warmup).
     */
    private fun getActualSimTime(): Double {
        return Sim.time() - warmupEndTime
    }

    /**
     * End of simulation event.
     */
    private inner class EndOfSimulation : Event() {
        override fun actions() {
            // Final update of queue and busy statistics
            for (qIdx in 0 until numServiceNodes) {
                // For PS queues, use the PS-specific busy stats update
                if (isPSScheduling(schedStrategies[qIdx])) {
                    updatePSBusyStats(qIdx)
                    for (k in 0 until numClasses) {
                        updateQueueStats(qIdx, k)
                    }
                } else {
                    for (k in 0 until numClasses) {
                        updateQueueStats(qIdx, k)
                        updateBusyStats(qIdx, k)
                    }
                }
            }

            // Apply MSER-5 truncation to determine warmup period
            if (mserEnabled) {
                applyMSER5Truncation()
            }

            // Close trace writer
            closeTracing()

            // Close Logger file writers
            closeLoggers()

            if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                // Print final sample count before newline
                System.out.printf("\b\b\b\b\b\b\b %6d", options.samples)
                println()
            }

            Sim.stop()
        }
    }

    /**
     * Progress reporting event.
     */
    private inner class ProgressEvent : Event() {
        var printedStart = false

        override fun actions() {
            val t = Sim.time()
            val samples = t.toInt()
            val totalSamples = options.samples

            if (options.verbose == VerboseLevel.STD || options.verbose == VerboseLevel.DEBUG) {
                if (!printedStart) {
                    System.out.printf("DES samples: %6d ", samples)
                    System.out.flush()
                    printedStart = true
                } else {
                    System.out.printf("\b\b\b\b\b\b\b %6d", samples)
                    System.out.flush()
                }
            }

            val updateInterval = totalSamples.toDouble() / 50.0
            if (updateInterval > 0 && t + updateInterval <= totalSamples) {
                schedule(updateInterval)
            }
        }
    }

    // ==================== Impatience Event Classes ====================

    /**
     * Reneging event - customer abandons queue due to expired patience.
     * Scheduled when a customer joins a queue with reneging configured.
     * Cancelled if the customer starts service before this event fires.
     */
    private inner class RenegingEvent(
        private val queueIdx: Int,
        private val customer: Customer,
        private val impatientKey: Triple<Int, Long, Int>
    ) : Event() {
        override fun actions() {
            // Check if customer is still waiting (may have started service already)
            val impatient = waitingImpatientCustomers.remove(impatientKey) ?: return

            val classId = customer.classId
            val currentTime = Sim.time()

            // Remove customer from wait queue
            val removed = waitQueues[queueIdx].removeIf {
                it.classId == classId &&
                it.systemArrivalTime == customer.systemArrivalTime
            }

            if (removed) {
                // Update queue length statistics before changing count
                updateQueueStats(queueIdx, classId)
                currentQueueLength[queueIdx][classId]--

                // Record reneging statistics
                if (warmupDone) {
                    renegedCustomers[queueIdx][classId]++
                    val waitTime = currentTime - customer.queueArrivalTime
                    totalRenegingWaitTime[queueIdx][classId] += waitTime
                }

                // Update FCR region stats if applicable
                val stationIdx = serviceStations[queueIdx]
                if (fcRegionIndices.isNotEmpty() && stationIdx < fcRegionIndices.size && fcRegionIndices[stationIdx] >= 0) {
                    val regionIdx = fcRegionIndices[stationIdx]
                    updateRegionTimeWeightedStats(regionIdx)
                    currentJobsInRegion[regionIdx][classId]--
                }

                // Log event if verbose
                if (options.verbose == VerboseLevel.DEBUG) {
                    println("RENEGE: Queue $queueIdx class $classId at time $currentTime, queue length now ${currentQueueLength[queueIdx][classId]}")
                }

                // Customer leaves system (no forwarding for reneged customers)
                // Note: Future enhancement could route reneged customers to a "balk sink"
            }
        }
    }

    /**
     * Retrial event - job in orbit attempts to re-enter queue.
     * Scheduled when a customer is rejected due to capacity and retrial is configured.
     */
    private inner class RetrialEvent(
        private val queueIdx: Int,
        private val orbitJob: OrbitJob
    ) : Event() {
        override fun actions() {
            val classId = orbitJob.customer.classId
            val currentTime = Sim.time()

            // Remove from orbit tracking
            orbitJobs[queueIdx].removeIf { it === orbitJob }

            // Update orbit statistics before changing count
            updateOrbitTimeStats(queueIdx, classId)
            currentOrbitSize[queueIdx][classId]--

            // Check if queue can accept the job now
            val currentTotal = getTotalCustomersAtStation(queueIdx)
            val capacity = bufferCapacities[queueIdx]

            if (currentTotal < capacity) {
                // Queue has capacity - customer enters
                val newCustomer = Customer(
                    classId = orbitJob.customer.classId,
                    priority = orbitJob.customer.priority,
                    systemArrivalTime = orbitJob.customer.systemArrivalTime,
                    queueArrivalTime = currentTime,  // New queue arrival time
                    randomRank = siroRng.nextDouble(),
                    jobId = orbitJob.customer.jobId,
                    absoluteDeadline = orbitJob.customer.absoluteDeadline
                )

                if (warmupDone) {
                    retriedCustomers[queueIdx][classId]++
                }

                // Log successful retrial
                if (options.verbose == VerboseLevel.DEBUG) {
                    println("RETRIAL_SUCCESS: Queue $queueIdx class $classId at time $currentTime")
                }

                // Re-enter queue (this will handle reneging scheduling if configured)
                arriveAtQueueFromRetrial(queueIdx, newCustomer)
            } else {
                // Queue still full - check retry limit
                val newAttempts = orbitJob.retrialAttempts + 1
                if (orbitJob.maxAttempts >= 0 && newAttempts >= orbitJob.maxAttempts) {
                    // Max retries exceeded - drop customer
                    if (warmupDone) {
                        maxRetriesExceeded[queueIdx][classId]++
                    }

                    if (options.verbose == VerboseLevel.DEBUG) {
                        println("RETRIAL_DROPPED: Queue $queueIdx class $classId at time $currentTime after $newAttempts attempts")
                    }
                } else {
                    // Schedule another retry
                    scheduleRetrial(queueIdx, orbitJob.customer, newAttempts, orbitJob.maxAttempts)

                    if (options.verbose == VerboseLevel.DEBUG) {
                        println("RETRIAL_RESCHEDULE: Queue $queueIdx class $classId attempt $newAttempts at time $currentTime")
                    }
                }
            }
        }
    }

    /**
     * Select destination node and class based on routing strategy.
     * For self-looping classes, always returns the reference station node.
     * Handles Logger, Router, ClassSwitch, and Cache nodes as pass-through.
     * Supports implicit class switching via routing matrix and explicit via ClassSwitch nodes.
     */
    private fun selectDestination(fromNode: Int, classId: Int): RoutingResult {
        // If pass-through nodes exist (Logger, Router, ClassSwitch, Cache), use pass-through routing
        if (loggerNodes.isNotEmpty() || routerNodes.isNotEmpty() || classSwitchNodes.isNotEmpty() || cacheNodes.isNotEmpty()) {
            val jobId = nextJobId++
            return routeThroughPassthroughNodes(fromNode, classId, jobId)
        }

        // Direct routing with class switching support
        return selectDestinationWithClassSwitch(fromNode, classId)
    }

    /**
     * Handle a job arriving at a Fork node.
     * Splits the parent job into multiple forked children and routes each to its destination.
     *
     * @param forkNodeIdx The Fork node index
     * @param parentJobId Unique ID for the parent job
     * @param classId Class of the parent job
     * @param systemArrivalTime When the parent job entered the system
     */
    private fun handleForkArrival(forkNodeIdx: Int, parentJobId: Long, classId: Int, systemArrivalTime: Double) {
        val forkListIdx = forkNodes.indexOf(forkNodeIdx)
        if (forkListIdx < 0) return

        val fanOut = forkFanOut[forkListIdx]
        val currentTime = Sim.time()

        // Get all output destinations from the Fork node
        val R = numClasses
        val I = numNodes
        val outputDestinations = mutableListOf<RoutingResult>()

        for (toNode in 0 until I) {
            for (toClass in 0 until R) {
                val prob = sn.rtnodes[forkNodeIdx * R + classId, toNode * R + toClass]
                if (prob > 0) {
                    outputDestinations.add(RoutingResult(toNode, toClass))
                }
            }
        }

        if (outputDestinations.isEmpty()) {
            // No valid outputs from fork - this shouldn't happen in a valid model
            return
        }

        // Calculate total number of forked tasks
        val totalTasks = outputDestinations.size * fanOut

        // Create fork job info for tracking synchronization at Join
        val forkInfo = ForkJobInfo(
            parentJobId = parentJobId,
            parentClassId = classId,
            parentSystemArrivalTime = systemArrivalTime,
            forkNodeIdx = forkNodeIdx,
            totalTasks = totalTasks
        )
        forkJobInfoMap[parentJobId] = forkInfo

        // Note: Queue length is now tracked when forked jobs arrive at Join buffer,
        // not when parent job enters fork section. This matches JMT semantics.

        // Create and route forked children (tasks)
        for (dest in outputDestinations) {
            for (taskIdx in 0 until fanOut) {
                val forkedJobId = nextForkedJobId++
                forkedJobParentMap[forkedJobId] = parentJobId

                val forkedJob = ForkedJob(
                    forkJobId = forkedJobId,
                    parentJobId = parentJobId,
                    classId = dest.destClassId,
                    priority = classPrio[dest.destClassId],
                    systemArrivalTime = systemArrivalTime,
                    queueArrivalTime = currentTime,
                    randomRank = siroRng.nextDouble()
                )

                // Route the forked job to its destination
                routeForkedJob(forkedJob, dest.destNode)
            }
        }
    }

    /**
     * Route a forked job to its destination node.
     * Handles Queue, Delay, Join, Sink, and pass-through nodes.
     */
    private fun routeForkedJob(forkedJob: ForkedJob, destNode: Int) {
        // Handle pass-through nodes (Logger, Router, ClassSwitch)
        var currentNode = destNode
        var currentClass = forkedJob.classId
        var maxIterations = 100

        while (maxIterations > 0) {
            maxIterations--

            if (sinkNodes.contains(currentNode)) {
                // Forked job reaches sink - this is unusual for a fork-join pattern
                // but handle it by removing the forked job (it won't reach join)
                return
            }

            if (joinNodes.contains(currentNode)) {
                // Forked job arrives at Join node
                handleJoinArrival(currentNode, forkedJob)
                return
            }

            if (forkNodes.contains(currentNode)) {
                // Nested fork - create sub-fork
                handleForkArrival(currentNode, forkedJob.forkJobId, currentClass, forkedJob.systemArrivalTime)
                return
            }

            val queueIdx = serviceNodes.indexOf(currentNode)
            if (queueIdx >= 0) {
                // Arrived at a service node (Queue or Delay)
                val customer = Customer(
                    currentClass,
                    forkedJob.priority,
                    forkedJob.systemArrivalTime,
                    Sim.time(),
                    forkedJob.randomRank,
                    -1.0
                )
                // Store the mapping from customer to forked job for tracking
                arriveAtQueueForked(queueIdx, customer, forkedJob)
                return
            }

            // Handle pass-through nodes
            if (loggerNodes.contains(currentNode)) {
                logJobPassage(currentNode, currentClass, forkedJob.forkJobId)
                val result = selectDestinationWithClassSwitch(currentNode, currentClass)
                if (result.destNode < 0) return
                currentNode = result.destNode
                currentClass = result.destClassId
                continue
            }

            if (routerNodes.contains(currentNode)) {
                val result = selectDestinationWithClassSwitch(currentNode, currentClass)
                if (result.destNode < 0) return
                currentNode = result.destNode
                currentClass = result.destClassId
                continue
            }

            if (classSwitchNodes.contains(currentNode)) {
                // ClassSwitch: the rtnodes matrix already encodes class switching
                // Entry format: rtnodes[CS*K+inputClass, dest*K+outputClass]
                // So we look up with currentClass (input) and get output class from result
                val result = selectDestinationWithClassSwitch(currentNode, currentClass)
                if (result.destNode < 0) return
                currentNode = result.destNode
                currentClass = result.destClassId
                continue
            }

            // Unknown node type - stop routing
            return
        }
    }

    /**
     * Handle arrival of a forked job at a service node.
     * Tracks the forked job for later routing to Join.
     */
    private fun arriveAtQueueForked(queueIdx: Int, customer: Customer, forkedJob: ForkedJob) {
        // For PS queues, pass forked job directly (PS customers store forked job internally)
        if (isPSScheduling(schedStrategies[queueIdx])) {
            arriveAtPSQueue(queueIdx, customer, forkedJob)
            return
        }

        // For non-PS queues, record forked job mapping for departure handling
        val key = Pair(queueIdx, customer.queueArrivalTime)
        forkedCustomerMap[key] = forkedJob

        // Use normal queue arrival handling
        arriveAtQueue(queueIdx, customer)
    }

    // Map to track forked jobs by queue and arrival time
    private val forkedCustomerMap = mutableMapOf<Pair<Int, Double>, ForkedJob>()

    /**
     * Handle a forked job arriving at a Join node.
     * Checks synchronization condition and routes parent job when complete.
     *
     * @param joinNodeIdx The Join node index
     * @param forkedJob The arriving forked job
     */
    private fun handleJoinArrival(joinNodeIdx: Int, forkedJob: ForkedJob) {
        val parentJobId = forkedJob.parentJobId
        val forkInfo = forkJobInfoMap[parentJobId] ?: return

        val joinListIdx = joinNodes.indexOf(joinNodeIdx)
        if (joinListIdx < 0) return

        // Use forked job's class for statistics and routing from Join
        // This is required for class-switching fork-join models where forked jobs
        // have different classes than the parent job
        val classId = forkedJob.classId
        val currentTime = Sim.time()

        // Record first arrival time at Join for this parent job
        if (forkInfo.firstJoinArrivalTime < 0) {
            forkInfo.firstJoinArrivalTime = currentTime
        }
        // Record this forked job's arrival time and class for per-job response time calculation
        forkInfo.forkedJobJoinArrivalTimes.add(currentTime)
        forkInfo.forkedJobJoinClasses.add(classId)

        // JMT semantics: Increment queue length when forked job arrives at Join buffer
        // (counts forked jobs waiting for synchronization, not parent jobs)
        if (warmupDone) {
            val elapsed = currentTime - lastJoinUpdateTime[joinListIdx][classId]
            if (elapsed > 0) {
                totalJoinQueueTime[joinListIdx][classId] += currentJoinQueueLength[joinListIdx][classId] * elapsed
            }
            lastJoinUpdateTime[joinListIdx][classId] = currentTime
            currentJoinQueueLength[joinListIdx][classId]++
            // Track arrivals for arrival rate calculation
            arrivedAtJoin[joinListIdx][classId]++
        }

        // Increment completed tasks count
        forkInfo.completedTasks++

        // Get join strategy and required count for this class
        val strategy = joinStrategies[joinListIdx][classId]
        val required = joinRequired[joinListIdx][classId]

        // Check if synchronization condition is met
        val syncComplete = when (strategy) {
            JoinStrategy.STD -> {
                // Wait for all forked tasks
                forkInfo.completedTasks >= forkInfo.totalTasks
            }
            JoinStrategy.PARTIAL, JoinStrategy.Quorum -> {
                // Wait for specified number of tasks
                if (required > 0) {
                    forkInfo.completedTasks >= required
                } else {
                    // -1 means wait for all
                    forkInfo.completedTasks >= forkInfo.totalTasks
                }
            }
            JoinStrategy.Guard -> {
                // Guard strategy - per-class requirements (treat as STD for now)
                forkInfo.completedTasks >= forkInfo.totalTasks
            }
            else -> {
                // Default: wait for all tasks
                forkInfo.completedTasks >= forkInfo.totalTasks
            }
        }

        if (syncComplete) {
            // Record Join statistics before removing fork info (only after warmup)
            if (warmupDone) {
                // JMT semantics: Record response time for each forked job at Join per class
                // Response time = time from arrival at Join to sync completion
                for (i in forkInfo.forkedJobJoinArrivalTimes.indices) {
                    val arrivalTime = forkInfo.forkedJobJoinArrivalTimes[i]
                    val forkedClassId = forkInfo.forkedJobJoinClasses[i]
                    val joinResponseTime = currentTime - arrivalTime
                    joinResponseTimeTally[joinListIdx][forkedClassId].add(joinResponseTime)
                }

                // Count forked jobs per class for proper queue length decrement
                val classCountMap = mutableMapOf<Int, Int>()
                for (forkedClassId in forkInfo.forkedJobJoinClasses) {
                    classCountMap[forkedClassId] = (classCountMap[forkedClassId] ?: 0) + 1
                }

                // JMT semantics: Update time-weighted stats and decrement queue length per class
                for ((forkedClassId, count) in classCountMap) {
                    val elapsed = currentTime - lastJoinUpdateTime[joinListIdx][forkedClassId]
                    if (elapsed > 0) {
                        totalJoinQueueTime[joinListIdx][forkedClassId] += currentJoinQueueLength[joinListIdx][forkedClassId] * elapsed
                    }
                    lastJoinUpdateTime[joinListIdx][forkedClassId] = currentTime

                    // Decrement by the number of forked jobs of this class that were waiting at Join
                    currentJoinQueueLength[joinListIdx][forkedClassId] -= count
                    if (currentJoinQueueLength[joinListIdx][forkedClassId] < 0) {
                        currentJoinQueueLength[joinListIdx][forkedClassId] = 0
                    }
                }

                // Increment completion count for parent class (one per parent job sync)
                // JMT shows throughput at Join for the original parent class, not forked classes
                joinCompletions[joinListIdx][forkInfo.parentClassId]++
            }

            // Remove fork info as synchronization is complete
            forkJobInfoMap.remove(parentJobId)

            // JMT semantics: When a job is forked, forked copies may change class during
            // processing, but at the Join the synchronized job RETURNS to the original
            // parent class. Routing from Join uses the parent class, not the forked classes.
            val routingResult = selectDestination(joinNodeIdx, forkInfo.parentClassId)
            var destNode = routingResult.destNode
            var destClassId = routingResult.destClassId

            // JMT compatibility: If routing from Join goes through a ClassSwitch but the class
            // doesn't actually change, this is not a valid route (the routing matrix only has
            // class-switching routes defined for the forked classes, not the parent class).
            // In this case, the job should be dropped (similar to JMT behavior).
            if (destNode >= 0 && !sinkNodes.contains(destNode) && destClassId == forkInfo.parentClassId) {
                // Check if the only route is through a ClassSwitch that doesn't transform the class
                // This happens when the model defines routes from Join for forked classes (Class2, Class3)
                // but not for the parent class (Class1), yet LINE creates default routes for all classes.
                // If the destination class equals parent class, no class switching occurred.
                val isQueueOrDelay = serviceNodes.contains(destNode)
                if (isQueueOrDelay) {
                    // Destination is a service node and class didn't change - check if service is available
                    val stationIdx = sn.nodeToStation[destNode].toInt()
                    val rates = sn.rates?.get(stationIdx, destClassId)
                    if (rates != null && rates.isNaN()) {
                        // Service is disabled for this class - treat as no valid route
                        // This matches JMT behavior where jobs are dropped if no valid route exists
                        destNode = -1
                    }
                }
            }

            if (destNode >= 0 && sinkNodes.contains(destNode)) {
                // Parent job leaves system - use parent's original class for statistics
                val systemResponseTime = Sim.time() - forkInfo.parentSystemArrivalTime
                systemResponseTimeTally[forkInfo.parentClassId].add(systemResponseTime)
                systemCompletedCustomers[forkInfo.parentClassId]++
            } else if (destNode >= 0) {
                // Route to next node
                routeParentJobFromJoin(destNode, destClassId, forkInfo)
            }
        }
        // If sync not complete, forked job is "absorbed" at Join (no further action)
    }

    /**
     * Route the parent job from Join node to next destination.
     */
    private fun routeParentJobFromJoin(destNode: Int, destClassId: Int, forkInfo: ForkJobInfo) {
        val currentTime = Sim.time()

        // Handle pass-through and various node types
        var currentNode = destNode
        var currentClass = destClassId
        var maxIterations = 100

        while (maxIterations > 0) {
            maxIterations--

            if (sinkNodes.contains(currentNode)) {
                // Parent job leaves system
                val systemResponseTime = currentTime - forkInfo.parentSystemArrivalTime
                systemResponseTimeTally[forkInfo.parentClassId].add(systemResponseTime)
                systemCompletedCustomers[forkInfo.parentClassId]++
                return
            }

            if (forkNodes.contains(currentNode)) {
                // Another fork - create new fork operation
                val newParentJobId = nextJobId++
                handleForkArrival(currentNode, newParentJobId, currentClass, forkInfo.parentSystemArrivalTime)
                return
            }

            if (joinNodes.contains(currentNode)) {
                // Nested join - shouldn't happen in normal patterns, but handle gracefully
                val routingResult = selectDestinationWithClassSwitch(currentNode, currentClass)
                if (routingResult.destNode < 0) return
                currentNode = routingResult.destNode
                currentClass = routingResult.destClassId
                continue
            }

            if (placeNodes.contains(currentNode)) {
                // Route to Place node
                handlePlaceArrival(currentNode, currentClass, forkInfo.parentSystemArrivalTime)
                return
            }

            if (transitionNodes.contains(currentNode)) {
                // Route to Transition node
                checkAndFireTransitions()
                return
            }

            val queueIdx = serviceNodes.indexOf(currentNode)
            if (queueIdx >= 0) {
                // Route to service node
                val customer = Customer(
                    currentClass,
                    classPrio[currentClass],
                    forkInfo.parentSystemArrivalTime,
                    currentTime,
                    siroRng.nextDouble()
                )
                arriveAtQueue(queueIdx, customer)
                return
            }

            // Handle pass-through nodes
            if (loggerNodes.contains(currentNode)) {
                logJobPassage(currentNode, currentClass, forkInfo.parentJobId)
                val result = selectDestinationWithClassSwitch(currentNode, currentClass)
                if (result.destNode < 0) return
                currentNode = result.destNode
                currentClass = result.destClassId
                continue
            }

            if (routerNodes.contains(currentNode)) {
                val result = selectDestinationWithClassSwitch(currentNode, currentClass)
                if (result.destNode < 0) return
                currentNode = result.destNode
                currentClass = result.destClassId
                continue
            }

            if (classSwitchNodes.contains(currentNode)) {
                // ClassSwitch: the rtnodes matrix already encodes class switching
                // Entry format: rtnodes[CS*K+inputClass, dest*K+outputClass]
                // So we look up with currentClass (input) and get output class from result
                val result = selectDestinationWithClassSwitch(currentNode, currentClass)
                if (result.destNode < 0) return
                currentNode = result.destNode
                currentClass = result.destClassId
                continue
            }

            // Unknown node type
            return
        }
    }

    /**
     * Check if a destination node is a Fork node and handle accordingly.
     * Returns true if the destination was a Fork and handling was done.
     */
    private fun handleDestinationIfFork(destNode: Int, classId: Int, systemArrivalTime: Double): Boolean {
        if (forkNodes.contains(destNode)) {
            val parentJobId = nextJobId++
            handleForkArrival(destNode, parentJobId, classId, systemArrivalTime)
            return true
        }
        return false
    }

    /**
     * Check if a destination node is a Join node.
     * This should not be called directly for regular customers - only forked jobs arrive at Join.
     * Returns true if the destination is a Join node.
     */
    private fun isJoinNode(destNode: Int): Boolean {
        return joinNodes.contains(destNode)
    }

    /**
     * Handle a job arriving at Join when forked job tracking was lost.
     * Attempts to find matching parent by system arrival time.
     */
    private fun handleUnknownJoinArrival(joinNodeIdx: Int, classId: Int, systemArrivalTime: Double) {
        val joinListIdx = joinNodes.indexOf(joinNodeIdx)
        if (joinListIdx < 0) return

        // Find the fork that this join is paired with
        val forkListIdx = joinToForkMap[joinListIdx]
        if (forkListIdx < 0) return

        val forkNodeIdx = forkNodes[forkListIdx]

        // Look for a matching fork info by system arrival time
        val matchingEntry = forkJobInfoMap.entries.find { (_, info) ->
            info.forkNodeIdx == forkNodeIdx &&
            kotlin.math.abs(info.parentSystemArrivalTime - systemArrivalTime) < 1e-9
        }

        if (matchingEntry != null) {
            val (parentJobId, forkInfo) = matchingEntry
            // Create a synthetic forked job for join handling
            val syntheticForkedJob = ForkedJob(
                forkJobId = nextForkedJobId++,
                parentJobId = parentJobId,
                classId = classId,
                priority = classPrio[classId],
                systemArrivalTime = systemArrivalTime,
                queueArrivalTime = Sim.time(),
                randomRank = 0.0
            )
            handleJoinArrival(joinNodeIdx, syntheticForkedJob)
        }
        // If no matching fork info, the job is lost (shouldn't happen in valid models)
    }

    /**
     * Finds a free server for the given queue.
     * Considers both busy status and server state (for setup/delayoff).
     * Priority: ACTIVE idle > DELAYOFF > OFF (if setup enabled)
     * @param queueIdx The service node index
     * @return The server ID if found, -1 otherwise
     */
    private fun findFreeServer(queueIdx: Int): Int {
        // If setup/delayoff not enabled, use original logic
        if (!hasSetupDelayoff[queueIdx]) {
            for (i in serverBusy[queueIdx].indices) {
                if (!serverBusy[queueIdx][i]) return i
            }
            return -1
        }

        // With setup/delayoff enabled, check server states
        // First pass: Find ACTIVE idle servers (highest priority)
        for (i in serverBusy[queueIdx].indices) {
            if (!serverBusy[queueIdx][i] && serverState[queueIdx][i] == ServerState.ACTIVE) {
                return i
            }
        }

        // Second pass: Find DELAYOFF servers (can resume immediately)
        for (i in serverBusy[queueIdx].indices) {
            if (serverState[queueIdx][i] == ServerState.DELAYOFF) {
                return i
            }
        }

        // Third pass: Find OFF servers (requires setup)
        for (i in serverBusy[queueIdx].indices) {
            if (serverState[queueIdx][i] == ServerState.OFF) {
                return i
            }
        }

        // No available servers (all are SETUP or busy ACTIVE)
        return -1
    }

    // ==================== Heterogeneous Server Selection ====================

    /**
     * Data class for server selection result.
     */
    private data class ServerSelection(
        val serverId: Int,      // Global server ID (-1 if none available)
        val serverTypeId: Int   // Server type ID (-1 if homogeneous or none available)
    )

    /**
     * Finds a free server compatible with the given job class.
     * For heterogeneous queues, applies the configured scheduling policy.
     * For homogeneous queues, delegates to findFreeServer().
     */
    private fun findFreeServerForClass(queueIdx: Int, classId: Int): ServerSelection {
        // Check if this queue has heterogeneous servers
        if (numServerTypes[queueIdx] == 0) {
            // Homogeneous queue - use original logic
            val serverId = findFreeServer(queueIdx)
            return ServerSelection(serverId, -1)
        }

        // Get compatible server types with available capacity
        val compatibleTypes = getCompatibleServerTypes(queueIdx, classId)
        if (compatibleTypes.isEmpty()) {
            return ServerSelection(-1, -1)
        }

        // Apply scheduling policy to select server type
        val policy = heteroSchedPolicy[queueIdx] ?: HeteroSchedPolicy.ORDER
        val selectedType = selectServerType(queueIdx, classId, compatibleTypes, policy)

        if (selectedType < 0) {
            return ServerSelection(-1, -1)
        }

        // Find free server within selected type
        val serverId = findFreeServerOfType(queueIdx, selectedType)
        return ServerSelection(serverId, selectedType)
    }

    /**
     * Gets list of server types compatible with given class that have available capacity.
     */
    private fun getCompatibleServerTypes(queueIdx: Int, classId: Int): List<Int> {
        val nTypes = numServerTypes[queueIdx]
        val result = mutableListOf<Int>()

        for (typeId in 0 until nTypes) {
            if (serverCompat[queueIdx][typeId][classId]) {
                val totalServers = serversPerType[queueIdx][typeId]
                val busyServers = busyCountPerType[queueIdx][typeId]
                if (busyServers < totalServers) {
                    result.add(typeId)
                }
            }
        }

        return result
    }

    /**
     * Finds a free server within a specific server type.
     * Handles setup/delayoff state if enabled.
     */
    private fun findFreeServerOfType(queueIdx: Int, typeId: Int): Int {
        // Find the range of server IDs for this type
        var startId = 0
        for (t in 0 until typeId) {
            startId += serversPerType[queueIdx][t]
        }
        val endId = startId + serversPerType[queueIdx][typeId]

        if (!hasSetupDelayoff[queueIdx]) {
            // Simple case: find first free server of this type
            for (serverId in startId until endId) {
                if (!serverBusy[queueIdx][serverId]) {
                    return serverId
                }
            }
        } else {
            // With setup/delayoff: prioritize ACTIVE > DELAYOFF > OFF
            // First pass: ACTIVE idle servers
            for (serverId in startId until endId) {
                if (!serverBusy[queueIdx][serverId] && serverState[queueIdx][serverId] == ServerState.ACTIVE) {
                    return serverId
                }
            }
            // Second pass: DELAYOFF servers
            for (serverId in startId until endId) {
                if (serverState[queueIdx][serverId] == ServerState.DELAYOFF) {
                    return serverId
                }
            }
            // Third pass: OFF servers
            for (serverId in startId until endId) {
                if (serverState[queueIdx][serverId] == ServerState.OFF) {
                    return serverId
                }
            }
        }
        return -1
    }

    /**
     * Selects server type according to the heterogeneous scheduling policy.
     */
    private fun selectServerType(
        queueIdx: Int,
        classId: Int,
        compatibleTypes: List<Int>,
        policy: HeteroSchedPolicy
    ): Int {
        if (compatibleTypes.isEmpty()) return -1
        if (compatibleTypes.size == 1) return compatibleTypes[0]

        return when (policy) {
            HeteroSchedPolicy.ORDER -> selectServerTypeORDER(compatibleTypes)
            HeteroSchedPolicy.ALIS -> selectServerTypeALIS(queueIdx, compatibleTypes)
            HeteroSchedPolicy.ALFS -> selectServerTypeALFS(queueIdx, compatibleTypes)
            HeteroSchedPolicy.FAIRNESS -> selectServerTypeFAIRNESS(queueIdx, compatibleTypes)
            HeteroSchedPolicy.FSF -> selectServerTypeFSF(queueIdx, classId, compatibleTypes)
            HeteroSchedPolicy.RAIS -> selectServerTypeRAIS(compatibleTypes)
        }
    }

    /** ORDER: First compatible server type in definition order */
    private fun selectServerTypeORDER(compatibleTypes: List<Int>): Int {
        return compatibleTypes.minOrNull() ?: -1
    }

    /** ALIS: Assign Longest Idle Server - round-robin among compatible types */
    private fun selectServerTypeALIS(queueIdx: Int, compatibleTypes: List<Int>): Int {
        val order = serverTypeOrder[queueIdx]
        for (typeId in order) {
            if (compatibleTypes.contains(typeId)) {
                // Move used type to end for round-robin
                order.remove(typeId)
                order.add(typeId)
                return typeId
            }
        }
        return compatibleTypes.firstOrNull() ?: -1
    }

    /** ALFS: Prefer types with fewer compatible classes (least flexible first) */
    private fun selectServerTypeALFS(queueIdx: Int, compatibleTypes: List<Int>): Int {
        // Use pre-sorted order (least flexible first)
        val sortedTypes = alfsOrder[queueIdx]
        for (typeId in sortedTypes) {
            if (compatibleTypes.contains(typeId)) {
                return typeId
            }
        }
        return compatibleTypes.firstOrNull() ?: -1
    }

    /** FAIRNESS: Simple round-robin among compatible types */
    private fun selectServerTypeFAIRNESS(queueIdx: Int, compatibleTypes: List<Int>): Int {
        // Same as ALIS for now
        return selectServerTypeALIS(queueIdx, compatibleTypes)
    }

    /** FSF: Fastest Server First - pick type with highest service rate for this class */
    private fun selectServerTypeFSF(queueIdx: Int, classId: Int, compatibleTypes: List<Int>): Int {
        var bestType = -1
        var bestRate = -1.0

        for (typeId in compatibleTypes) {
            val rate = heteroMus[queueIdx][typeId][classId]
            if (rate > bestRate && rate < Double.MAX_VALUE) {
                bestRate = rate
                bestType = typeId
            }
        }

        return if (bestType >= 0) bestType else compatibleTypes.firstOrNull() ?: -1
    }

    /** RAIS: Random Available Idle Server - random weighted by free server count */
    private fun selectServerTypeRAIS(compatibleTypes: List<Int>): Int {
        if (compatibleTypes.isEmpty()) return -1
        val idx = (routingRng.nextDouble() * compatibleTypes.size).toInt()
            .coerceIn(0, compatibleTypes.size - 1)
        return compatibleTypes[idx]
    }

    /**
     * Marks a server as busy and updates heterogeneous tracking.
     */
    private fun markServerBusy(queueIdx: Int, serverId: Int, typeId: Int = -1) {
        serverBusy[queueIdx][serverId] = true
        if (numServerTypes[queueIdx] > 0) {
            val effectiveTypeId = if (typeId >= 0) typeId else serverToType[queueIdx][serverId]
            if (effectiveTypeId >= 0 && effectiveTypeId < busyCountPerType[queueIdx].size) {
                busyCountPerType[queueIdx][effectiveTypeId]++
            }
        }
    }

    /**
     * Marks a server as idle and updates heterogeneous tracking.
     */
    private fun markServerIdle(queueIdx: Int, serverId: Int) {
        serverBusy[queueIdx][serverId] = false
        if (numServerTypes[queueIdx] > 0 && serverId < serverToType[queueIdx].size) {
            val typeId = serverToType[queueIdx][serverId]
            if (typeId >= 0 && typeId < busyCountPerType[queueIdx].size) {
                busyCountPerType[queueIdx][typeId]--
            }
        }
    }

    // ==================== End Heterogeneous Server Selection ====================

    private fun generateServiceTime(queueIdx: Int, classId: Int): Double {
        val procType = serviceProcessType[queueIdx][classId]

        // Check for disabled service - job should not request service from disabled class-station pair
        if (procType == ProcessType.DISABLED) {
            val stationName = sn.stations[serviceStations[queueIdx]].name
            val className = sn.jobclasses[classId].name
            throw RuntimeException("DES: Job of class '$className' requested service at station '$stationName' " +
                "but service is DISABLED for this class-station pair. Check model routing.")
        }

        // Generate base service time
        var baseServiceTime = 0.0

        // Check for Replayer/trace distribution first
        if (procType == ProcessType.REPLAYER) {
            val traceSampler = serviceTraceSamplers[queueIdx][classId]
            if (traceSampler != null) {
                baseServiceTime = traceSampler.nextSample()
            }
        } else {
            val gen = serviceGens[queueIdx][classId]
            if (gen != null) {
                baseServiceTime = gen.nextDouble()
            } else if (procType == ProcessType.MMAP) {
                // MMAP: sample service time using mmap_sample (ignores marks for services)
                val proc = serviceProc[queueIdx][classId]
                val rng = serviceRng[queueIdx][classId]
                if (proc != null && rng != null && proc.size() >= 2) {
                    val mmapResult = mmap_sample(proc, 1, rng)
                    baseServiceTime = mmapResult.samples[0]
                }
            } else if (procType == ProcessType.PH || procType == ProcessType.APH ||
                procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN ||
                procType == ProcessType.COX2 || procType == ProcessType.MAP || procType == ProcessType.MMPP2 ||
                procType == ProcessType.ME || procType == ProcessType.RAP) {
                // Fallback to map_sample for PH/MAP/MMPP/ME/RAP distributions
                val proc = serviceProc[queueIdx][classId]
                val rng = serviceRng[queueIdx][classId]
                if (proc != null && rng != null && proc.size() >= 2) {
                    val D0 = proc.get(0)
                    val D1 = proc.get(1)
                    val samples = map_sample(D0, D1, 1, rng)
                    baseServiceTime = samples[0]
                }
            }
        }

        // Apply load-dependent scaling if this station has load dependence
        // Note: For PS queues, load-dependent scaling is handled in the PS rate calculation
        // (via getEffectivePSServerCount), so we skip scaling here for PS queues.
        if (isLoadDependent[queueIdx] && !isPSScheduling(schedStrategies[queueIdx])) {
            val totalJobs = getTotalCustomersAtStation(queueIdx)
            if (totalJobs > 0) {
                val scalingArray = lldScaling[queueIdx]
                if (scalingArray != null) {
                    // lldscaling[station, n-1] gives scaling factor when n jobs are present
                    // Cap at array length for populations exceeding the defined range
                    val scalingIdx = minOf(totalJobs - 1, scalingArray.size - 1)
                    val scalingFactor = scalingArray[scalingIdx]
                    // Scaling factor multiplies the rate, so divide the service time
                    if (scalingFactor > 0) {
                        baseServiceTime /= scalingFactor
                    }
                }
            }
        }

        // Apply limited joint class-dependent (LJCD) scaling if present
        // LJCD takes precedence over LJD as it provides per-class scaling
        // LJCD tables store throughput values X(n). Service time = 1/X matches MVA's STeff = 1/X.
        // For PS scheduling, the PS algorithm provides the queuing factor (1+Q).
        if (hasLjcd && ljcdScaling != null && ljcdCutoffs != null) {
            val stationScaling = ljcdScaling!![queueIdx]
            val cutoffs = ljcdCutoffs!![queueIdx]
            if (stationScaling != null && cutoffs != null) {
                val classScaling = stationScaling[classId]
                if (classScaling != null) {
                    // Get population vector INCLUDING the arriving job
                    // This is called BEFORE queue length is incremented, but service time
                    // should be based on population AFTER arrival (job is present during service)
                    val nvec = IntArray(numClasses) { c ->
                        currentQueueLength[queueIdx][c] + (if (c == classId) 1 else 0)
                    }
                    // Clamp to cutoffs and compute linearized index
                    for (c in nvec.indices) {
                        nvec[c] = minOf(nvec[c], cutoffs[c])
                    }
                    val idx = ljdLinearize(nvec, cutoffs)
                    if (idx < classScaling.size) {
                        val scalingFactor = classScaling[idx]
                        // Use 1/X as service time (matches MVA's STeff = 1/X)
                        // PS sharing will provide the queuing factor
                        if (scalingFactor > 1e-10) {
                            baseServiceTime = baseServiceTime / scalingFactor
                        }
                        // else: X near 0 means no service, keep base time
                    }
                }
            }
        } else if (hasLjd && ljdScaling != null && ljdCutoffs != null) {
            // Apply limited joint-dependent (LJD) scaling if present
            // LJD uses a single scaling table for all classes
            // Unlike LJCD (which stores throughput values), LJD scaling is a direct
            // time multiplier that multiplies service time (same as in MVA: STeff = ST * ljdterm)
            val scaling = ljdScaling!![queueIdx]
            val cutoffs = ljdCutoffs!![queueIdx]
            if (scaling != null && cutoffs != null) {
                // Get population vector INCLUDING the arriving job
                // This is called BEFORE queue length is incremented, but service time
                // should be based on population AFTER arrival (job is present during service)
                val nvec = IntArray(numClasses) { c ->
                    currentQueueLength[queueIdx][c] + (if (c == classId) 1 else 0)
                }
                // Clamp to cutoffs and compute linearized index
                for (c in nvec.indices) {
                    nvec[c] = minOf(nvec[c], cutoffs[c])
                }
                val idx = ljdLinearize(nvec, cutoffs)
                if (idx < scaling.size) {
                    val scalingFactor = scaling[idx]
                    // LJD scaling is a time multiplier - multiply service time by it
                    if (scalingFactor > 0) {
                        baseServiceTime *= scalingFactor
                    }
                }
            }
        }

        // Track sample for control variates if enabled and after warmup
        if (useControlVariates && warmupDone && baseServiceTime > 0) {
            serviceSampleSum[queueIdx][classId] += baseServiceTime
            serviceSampleCount[queueIdx][classId]++
        }

        return baseServiceTime
    }

    /**
     * Generates service time for heterogeneous servers.
     * Uses server-type-specific distribution if available.
     * Falls back to homogeneous generation if not heterogeneous or typeId is -1.
     */
    private fun generateHeteroServiceTime(queueIdx: Int, classId: Int, typeId: Int): Double {
        // If not heterogeneous or typeId is -1, use homogeneous generation
        if (numServerTypes[queueIdx] == 0 || typeId < 0) {
            return generateServiceTime(queueIdx, classId)
        }

        val procType = heteroServiceProcType[queueIdx][typeId][classId]

        // Check for disabled service
        if (procType == ProcessType.DISABLED) {
            // Fall back to homogeneous (which will throw if also disabled)
            return generateServiceTime(queueIdx, classId)
        }

        var baseServiceTime = 0.0

        // Try to use the heterogeneous generator
        val gen = heteroServiceGens[queueIdx][typeId][classId]
        if (gen != null) {
            baseServiceTime = gen.nextDouble()
        } else if (procType == ProcessType.MMAP) {
            // Use mmap_sample for MMAP distributions
            val proc = heteroServiceProc[queueIdx][typeId][classId]
            val rng = heteroServiceRng[queueIdx][typeId][classId]
            if (proc != null && rng != null && proc.size() >= 2) {
                val mmapResult = mmap_sample(proc, 1, rng)
                baseServiceTime = mmapResult.samples[0]
            }
        } else if (procType == ProcessType.PH || procType == ProcessType.APH ||
            procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN ||
            procType == ProcessType.COX2 || procType == ProcessType.MAP) {
            // Use map_sample for PH/MAP distributions
            val proc = heteroServiceProc[queueIdx][typeId][classId]
            val rng = heteroServiceRng[queueIdx][typeId][classId]
            if (proc != null && rng != null && proc.size() >= 2) {
                val D0 = proc.get(0)
                val D1 = proc.get(1)
                val samples = map_sample(D0, D1, 1, rng)
                baseServiceTime = samples[0]
            }
        } else {
            // Fallback to rate-based exponential
            val rate = heteroMus[queueIdx][typeId][classId]
            if (rate > 0 && rate < Double.MAX_VALUE) {
                baseServiceTime = -Math.log(routingRng.nextDouble()) / rate
            }
        }

        // Apply load-dependent scaling if this station has load dependence
        if (isLoadDependent[queueIdx] && !isPSScheduling(schedStrategies[queueIdx])) {
            val totalJobs = getTotalCustomersAtStation(queueIdx)
            if (totalJobs > 0) {
                val scalingArray = lldScaling[queueIdx]
                if (scalingArray != null) {
                    val scalingIdx = minOf(totalJobs - 1, scalingArray.size - 1)
                    val scalingFactor = scalingArray[scalingIdx]
                    if (scalingFactor > 0) {
                        baseServiceTime /= scalingFactor
                    }
                }
            }
        }

        // Track sample for control variates if enabled and after warmup
        if (useControlVariates && warmupDone && baseServiceTime > 0) {
            serviceSampleSum[queueIdx][classId] += baseServiceTime
            serviceSampleCount[queueIdx][classId]++
        }

        return baseServiceTime
    }

    private fun generateInterarrivalTime(sourceIdx: Int, classId: Int): Double {
        var interarrivalTime = 0.0

        // Check for Replayer/trace distribution first
        val procType = arrivalProcessType[sourceIdx][classId]
        if (procType == ProcessType.REPLAYER) {
            val traceSampler = arrivalTraceSamplers[sourceIdx][classId]
            if (traceSampler != null) {
                interarrivalTime = traceSampler.nextSample()
            }
        } else {
            val gen = arrivalGens[sourceIdx][classId]
            if (gen != null) {
                interarrivalTime = gen.nextDouble()
            } else if (procType == ProcessType.BMAP) {
                // BMAP: sample inter-arrival time and batch size
                val proc = arrivalProc[sourceIdx][classId]
                val rng = arrivalRng[sourceIdx][classId]
                if (proc != null && rng != null && proc.size() >= 3) {
                    val samples = bmap_sample(proc, 1, rng)
                    interarrivalTime = samples[0].interarrivalTime
                    arrivalBatchSize[sourceIdx][classId] = samples[0].batchSize
                }
            } else if (procType == ProcessType.MMAP) {
                // MMAP: sample inter-arrival time using mmap_sample (ignores marks for per-class arrivals)
                val proc = arrivalProc[sourceIdx][classId]
                val rng = arrivalRng[sourceIdx][classId]
                if (proc != null && rng != null && proc.size() >= 2) {
                    val mmapResult = mmap_sample(proc, 1, rng)
                    interarrivalTime = mmapResult.samples[0]
                }
            } else if (procType == ProcessType.PH || procType == ProcessType.APH ||
                procType == ProcessType.HYPEREXP || procType == ProcessType.COXIAN ||
                procType == ProcessType.COX2 || procType == ProcessType.MAP || procType == ProcessType.MMPP2 ||
                procType == ProcessType.ME || procType == ProcessType.RAP) {
                // Fallback to map_sample for PH/MAP/MMPP/ME/RAP distributions
                val proc = arrivalProc[sourceIdx][classId]
                val rng = arrivalRng[sourceIdx][classId]
                if (proc != null && rng != null && proc.size() >= 2) {
                    val D0 = proc.get(0)
                    val D1 = proc.get(1)
                    val samples = map_sample(D0, D1, 1, rng)
                    interarrivalTime = samples[0]
                }
            }
        }

        // Track sample for control variates if enabled and after warmup
        if (useControlVariates && warmupDone && interarrivalTime > 0) {
            arrivalSampleSum[sourceIdx][classId] += interarrivalTime
            arrivalSampleCount[sourceIdx][classId]++
        }

        return interarrivalTime
    }

    /**
     * Generates a setup (cold start) time for a server at the given queue and class.
     * Setup time is sampled from the distribution configured for this queue/class.
     * @param queueIdx The service node index
     * @param classId The job class index
     * @return The sampled setup time
     */
    private fun generateSetupTime(queueIdx: Int, classId: Int): Double {
        val gen = setupGens[queueIdx][classId]
        return gen?.nextDouble() ?: 0.0
    }

    /**
     * Generates a delayoff (teardown) time for a server at the given queue and class.
     * Delayoff time is sampled from the distribution configured for this queue/class.
     * @param queueIdx The service node index
     * @param classId The job class index
     * @return The sampled delayoff time
     */
    private fun generateDelayoffTime(queueIdx: Int, classId: Int): Double {
        val gen = delayoffGens[queueIdx][classId]
        return gen?.nextDouble() ?: 0.0
    }

    /**
     * Checks if setup/delayoff is enabled for the given queue and class.
     * @param queueIdx The service node index
     * @param classId The job class index
     * @return true if setup/delayoff is enabled, false otherwise
     */
    private fun isSetupDelayoffEnabled(queueIdx: Int, classId: Int): Boolean {
        return hasSetupDelayoff[queueIdx] &&
                setupGens[queueIdx][classId] != null &&
                delayoffGens[queueIdx][classId] != null
    }

    /**
     * Initiates server setup (cold start) phase.
     * The server transitions to SETUP state and schedules a SetupCompletion event.
     * @param queueIdx The service node index
     * @param serverId The server ID
     * @param classId The job class that triggered the setup
     */
    private fun startServerSetup(queueIdx: Int, serverId: Int, classId: Int) {
        // Update server state to SETUP
        serverState[queueIdx][serverId] = ServerState.SETUP

        // Track which class triggered this setup
        serverLastClass[queueIdx][serverId] = classId

        // Generate setup time from distribution
        val setupTime = generateSetupTime(queueIdx, classId)

        // Update setup statistics
        updateSetupStats(queueIdx, classId)
        currentServersInSetup[queueIdx][classId]++
        lastSetupUpdateTime[queueIdx][classId] = Sim.time()

        // Schedule setup completion event
        val setupEvent = SetupCompletion(queueIdx, serverId, classId)
        setupEvent.schedule(setupTime)

        // Log setup start
        val stationIdx = serviceStations[queueIdx]
        logEvent("SETUP_START", stationIdx, classId, serverId, 0)
    }

    /**
     * Initiates server delayoff (teardown) phase.
     * The server transitions to DELAYOFF state and schedules a DelayoffCompletion event.
     * @param queueIdx The service node index
     * @param serverId The server ID
     * @param classId The job class associated with this delayoff
     */
    private fun startServerDelayoff(queueIdx: Int, serverId: Int, classId: Int) {
        // Update server state to DELAYOFF
        serverState[queueIdx][serverId] = ServerState.DELAYOFF

        // Track which class initiated this delayoff
        serverLastClass[queueIdx][serverId] = classId

        // Generate delayoff time from distribution
        val delayoffTime = generateDelayoffTime(queueIdx, classId)

        // Update delayoff statistics
        updateDelayoffStats(queueIdx, classId)
        currentServersInDelayoff[queueIdx][classId]++
        lastDelayoffUpdateTime[queueIdx][classId] = Sim.time()

        // Schedule delayoff completion event
        val delayoffEvent = DelayoffCompletion(queueIdx, serverId)
        delayoffEvent.schedule(delayoffTime)

        // Store event reference for potential cancellation
        pendingDelayoffEvents[queueIdx][serverId] = delayoffEvent

        // Log delayoff start
        val stationIdx = serviceStations[queueIdx]
        logEvent("DELAYOFF_START", stationIdx, classId, serverId, 0)
    }

    /**
     * Cancels an ongoing delayoff phase and returns the server to ACTIVE state.
     * This happens when a job arrives while the server is in the delayoff phase.
     * @param queueIdx The service node index
     * @param serverId The server ID
     */
    private fun cancelDelayoff(queueIdx: Int, serverId: Int) {
        // Cancel the pending delayoff event
        val event = pendingDelayoffEvents[queueIdx][serverId]
        if (event != null) {
            event.cancel()
            pendingDelayoffEvents[queueIdx][serverId] = null
        }

        // Get the class that was in delayoff
        val classId = serverLastClass[queueIdx][serverId]
        if (classId >= 0) {
            // Update delayoff statistics
            updateDelayoffStats(queueIdx, classId)
            currentServersInDelayoff[queueIdx][classId]--
        }

        // Update server state to ACTIVE (resume without setup)
        serverState[queueIdx][serverId] = ServerState.ACTIVE

        // Log delayoff cancellation
        val stationIdx = serviceStations[queueIdx]
        logEvent("DELAYOFF_CANCEL", stationIdx, classId, serverId, 0)
    }

    /**
     * Starts service for a customer on a specific server.
     * Updates statistics and schedules a departure event.
     * @param queueIdx The service node index
     * @param serverId The server ID
     * @param customer The customer to serve
     */
    private fun startService(queueIdx: Int, serverId: Int, customer: Customer, serverTypeId: Int = -1) {
        val classId = customer.classId

        // Cancel any scheduled reneging event (customer is starting service, no longer waiting)
        cancelRenegingIfScheduled(queueIdx, customer)

        // Mark server as busy (with heterogeneous tracking)
        markServerBusy(queueIdx, serverId, serverTypeId)
        customersInService[queueIdx]++

        // Store assigned server type on customer
        customer.assignedServerType = serverTypeId

        // Update busy time statistics
        updateBusyStats(queueIdx, classId)
        currentBusyServers[queueIdx][classId]++
        lastBusyUpdateTime[queueIdx][classId] = Sim.time()

        // Generate or use existing service time (heterogeneous-aware)
        val serviceTime = if (customer.serviceTime > 0) {
            customer.serviceTime
        } else {
            generateHeteroServiceTime(queueIdx, classId, serverTypeId)
        }

        // Schedule departure event
        val departureEvent = Departure(queueIdx, serverId, customer)
        departureEvent.schedule(serviceTime)

        // Track for signal-based removal (G-networks)
        if (hasNegativeSignals) {
            inServiceJobs[Pair(queueIdx, serverId)] = InServiceJob(customer, departureEvent)
        }
    }

    private fun updateQueueStats(queueIdx: Int, classId: Int) {
        val currentTime = Sim.time()
        val elapsed = currentTime - lastQueueUpdateTime[queueIdx][classId]
        if (elapsed > 0) {
            // currentQueueLength includes jobs waiting and in service
            // currentBlockedServers includes servers blocked waiting for SYNCH REPLY (at source)
            // basBlockedAtDest includes jobs BAS-blocked waiting to ENTER this queue
            // bbsBlockedAtDest includes jobs BBS-blocked waiting to ENTER this queue
            // fcrBlockedAtDest includes jobs FCR-blocked waiting to ENTER this queue
            val effectiveQueueLength = currentQueueLength[queueIdx][classId] +
                                       currentBlockedServers[queueIdx][classId] +  // SYNCH reply blocking (at source)
                                       basBlockedAtDest[queueIdx][classId] +        // BAS blocking (jobs trying to enter THIS queue)
                                       bbsBlockedAtDest[queueIdx][classId] +        // BBS blocking (jobs trying to enter THIS queue)
                                       fcrBlockedAtDest[queueIdx][classId]          // FCR blocking (jobs trying to enter THIS queue)
            totalQueueTime[queueIdx][classId] += effectiveQueueLength * elapsed
            lastQueueUpdateTime[queueIdx][classId] = currentTime
        }
    }

    // ==================== Place and Transition Handling ====================

    /**
     * Handle a token arriving at a Place node.
     * Adds token to place and triggers check for enabled transitions.
     */
    private fun handlePlaceArrival(placeNodeIdx: Int, classId: Int, systemArrivalTime: Double) {
        val placeListIdx = placeNodes.indexOf(placeNodeIdx)
        if (placeListIdx < 0) return

        // Update time-weighted token count before adding new token
        val currentTime = Sim.time()
        val elapsed = currentTime - lastPlaceUpdateTime[placeListIdx]
        if (elapsed > 0) {
            for (k in 0 until numClasses) {
                totalPlaceTokenTime[placeListIdx][k] += placeTokens[placeListIdx][k] * elapsed
            }
            lastPlaceUpdateTime[placeListIdx] = currentTime
        }

        // Add token to place
        placeTokens[placeListIdx][classId]++

        // Check and fire enabled transitions
        checkAndFireTransitions()
    }

    /**
     * Check all transitions for enabling conditions and fire if enabled.
     * Uses race semantics: all enabled modes sample their firing times concurrently,
     * and the mode with the minimum firing time wins the race.
     *
     * Supports multi-server modes: loops until no more firings are possible,
     * allowing multiple concurrent firings when servers and tokens are available.
     */
    private fun checkAndFireTransitions() {
        // Keep firing until no more firings are possible across all transitions
        var firedAny: Boolean
        do {
            firedAny = false
            for (transListIdx in transitionNodes.indices) {
                val modes = transitionModes[transListIdx]
                if (modes.isEmpty()) continue

                // Find enabled modes (respects server availability and token counts)
                val enabledModes = modes.filter { mode -> isTransitionModeEnabled(transListIdx, mode) }
                if (enabledModes.isEmpty()) continue

                // For priority-based selection, first filter by highest priority
                val maxPriority = enabledModes.maxOf { it.priority }
                val highPriorityModes = enabledModes.filter { it.priority == maxPriority }

                // Race semantics: pick winner among enabled modes
                val (winningMode, winningTime) = if (highPriorityModes.size == 1) {
                    val mode = highPriorityModes[0]
                    val time = if (mode.timingStrategy == TimingStrategy.IMMEDIATE) 0.0
                               else sampleTransitionFiringTime(transListIdx, mode.modeIdx)
                    Pair(mode, time)
                } else {
                    selectModeByRace(transListIdx, highPriorityModes)
                }

                if (winningMode != null) {
                    fireTransition(transListIdx, winningMode, winningTime)
                    firedAny = true
                }
            }
        } while (firedAny)
    }

    /**
     * Select a mode using race semantics: sample firing times for all modes
     * and return the mode with the minimum firing time along with that time.
     * This correctly handles concurrent mode competition in SPNs.
     */
    private fun selectModeByRace(transListIdx: Int, modes: List<TransitionModeInfo>): Pair<TransitionModeInfo?, Double?> {
        if (modes.isEmpty()) return Pair(null, null)
        if (modes.size == 1) {
            val mode = modes[0]
            val time = if (mode.timingStrategy == TimingStrategy.IMMEDIATE) 0.0
                       else sampleTransitionFiringTime(transListIdx, mode.modeIdx)
            return Pair(mode, time)
        }

        var minTime = Double.MAX_VALUE
        var winningMode: TransitionModeInfo? = null

        for (mode in modes) {
            val firingTime = if (mode.timingStrategy == TimingStrategy.IMMEDIATE) {
                // Immediate transitions fire at time 0, use weight for tie-breaking
                -mode.weight  // More negative = higher priority (larger weight wins)
            } else {
                // Timed transition: sample the firing time
                sampleTransitionFiringTime(transListIdx, mode.modeIdx)
            }

            if (firingTime < minTime) {
                minTime = firingTime
                winningMode = mode
            }
        }

        // For immediate transitions with negative times, use 0 for actual firing
        val actualTime = if (minTime < 0) 0.0 else minTime
        return Pair(winningMode, actualTime)
    }

    /**
     * Check if a transition mode is enabled (all enabling conditions met, no inhibiting conditions violated).
     */
    private fun isTransitionModeEnabled(transListIdx: Int, mode: TransitionModeInfo): Boolean {
        // Check server availability
        if (transitionInService[transListIdx][mode.modeIdx] >= mode.numServers) {
            return false
        }

        // Check enabling conditions: each input place must have required tokens
        for (placeListIdx in placeNodes.indices) {
            for (classIdx in 0 until numClasses) {
                val required = mode.enablingConditions[placeListIdx][classIdx]
                val available = placeTokens[placeListIdx][classIdx]
                if (required > 0 && available < required) {
                    return false
                }
            }
        }

        // Check inhibiting conditions: each input place must have fewer tokens than threshold
        for (placeListIdx in placeNodes.indices) {
            for (classIdx in 0 until numClasses) {
                val maxTokens = mode.inhibitingConditions[placeListIdx][classIdx]
                val currentTokens = placeTokens[placeListIdx][classIdx]
                if (maxTokens < Int.MAX_VALUE && currentTokens >= maxTokens) {
                    return false
                }
            }
        }

        return true
    }

    /**
     * Select a transition mode from enabled modes by priority and weight.
     */
    private fun selectTransitionMode(enabledModes: List<TransitionModeInfo>): TransitionModeInfo? {
        if (enabledModes.isEmpty()) return null
        if (enabledModes.size == 1) return enabledModes[0]

        // Find highest priority among enabled modes
        val maxPriority = enabledModes.maxOf { it.priority }
        val highPriorityModes = enabledModes.filter { it.priority == maxPriority }

        if (highPriorityModes.size == 1) return highPriorityModes[0]

        // Select by weight (probabilistic) among equally-prioritized modes
        val totalWeight = highPriorityModes.sumOf { it.weight }
        if (totalWeight <= 0) return highPriorityModes[0]

        val rand = siroRng.nextDouble() * totalWeight
        var cumWeight = 0.0
        for (mode in highPriorityModes) {
            cumWeight += mode.weight
            if (rand <= cumWeight) return mode
        }
        return highPriorityModes.last()
    }

    /**
     * Fire a transition mode: consume input tokens, schedule firing completion.
     * @param preSampledFiringTime Optional pre-sampled firing time from race selection.
     *        If provided, this time is used directly; otherwise, a new time is sampled.
     */
    private fun fireTransition(transListIdx: Int, mode: TransitionModeInfo, preSampledFiringTime: Double?) {
        val currentTime = Sim.time()

        // Update time-weighted token counts and consume input tokens from places
        for (placeListIdx in placeNodes.indices) {
            // Update statistics before consuming tokens
            val elapsed = currentTime - lastPlaceUpdateTime[placeListIdx]
            if (elapsed > 0) {
                for (k in 0 until numClasses) {
                    totalPlaceTokenTime[placeListIdx][k] += placeTokens[placeListIdx][k] * elapsed
                }
                lastPlaceUpdateTime[placeListIdx] = currentTime
            }

            // Consume tokens and track completions
            for (classIdx in 0 until numClasses) {
                val required = mode.enablingConditions[placeListIdx][classIdx]
                if (required > 0) {
                    placeTokens[placeListIdx][classIdx] -= required
                    placeCompletions[placeListIdx][classIdx] += required
                }
            }
        }

        // Update tokens in transit statistics before adding new tokens
        val transitElapsed = currentTime - lastTransitUpdateTime[transListIdx]
        if (transitElapsed > 0) {
            for (k in 0 until numClasses) {
                totalTransitTokenTime[transListIdx][k] += tokensInTransit[transListIdx][k] * transitElapsed
            }
            lastTransitUpdateTime[transListIdx] = currentTime
        }

        // Add consumed tokens to "in transit" and track per-place transit
        for (placeListIdx in placeNodes.indices) {
            for (classIdx in 0 until numClasses) {
                val required = mode.enablingConditions[placeListIdx][classIdx]
                if (required > 0) {
                    tokensInTransit[transListIdx][classIdx] += required
                    // Update place transit statistics before adding tokens
                    val placeTransitElapsed = currentTime - lastPlaceTransitUpdateTime[placeListIdx]
                    if (placeTransitElapsed > 0) {
                        for (k in 0 until numClasses) {
                            placeTransitTokenTime[placeListIdx][k] += placeTokensInTransit[placeListIdx][k] * placeTransitElapsed
                        }
                        lastPlaceTransitUpdateTime[placeListIdx] = currentTime
                    }
                    // Track that these tokens from this place are now in transit
                    placeTokensInTransit[placeListIdx][classIdx] += required
                }
            }
        }

        // Increment in-service count
        transitionInService[transListIdx][mode.modeIdx]++

        // Schedule firing completion
        if (mode.timingStrategy == TimingStrategy.IMMEDIATE) {
            // Immediate transition: fire at current time (schedule with 0 delay)
            TransitionFiring(transListIdx, mode.modeIdx).schedule(0.0)
        } else {
            // Timed transition: use pre-sampled time if available, otherwise sample new
            val firingTime = preSampledFiringTime ?: sampleTransitionFiringTime(transListIdx, mode.modeIdx)
            TransitionFiring(transListIdx, mode.modeIdx).schedule(firingTime)
        }
    }

    /**
     * Sample firing time for a timed transition.
     */
    private fun sampleTransitionFiringTime(transListIdx: Int, modeIdx: Int): Double {
        val gen = transitionFiringGens[transListIdx][modeIdx]
        return if (gen != null) {
            gen.nextDouble()
        } else {
            // Fallback: exponential with rate 1
            -FastMath.log(siroRng.nextDouble())
        }
    }

    /**
     * Transition firing completion event.
     */
    private inner class TransitionFiring(
        private val transListIdx: Int,
        private val modeIdx: Int
    ) : Event() {
        override fun actions() {
            val currentTime = Sim.time()
            val transNodeIdx = transitionNodes[transListIdx]
            val mode = transitionModes[transListIdx][modeIdx]

            // Update tokens in transit statistics before removing tokens
            val transitElapsed = currentTime - lastTransitUpdateTime[transListIdx]
            if (transitElapsed > 0) {
                for (k in 0 until numClasses) {
                    totalTransitTokenTime[transListIdx][k] += tokensInTransit[transListIdx][k] * transitElapsed
                }
                lastTransitUpdateTime[transListIdx] = currentTime
            }

            // Remove tokens from transit (tokens consumed from input places)
            for (placeListIdx in placeNodes.indices) {
                for (classIdx in 0 until numClasses) {
                    val required = mode.enablingConditions[placeListIdx][classIdx]
                    if (required > 0) {
                        tokensInTransit[transListIdx][classIdx] -= required
                        // Update place transit statistics before removing tokens
                        val placeTransitElapsed = currentTime - lastPlaceTransitUpdateTime[placeListIdx]
                        if (placeTransitElapsed > 0) {
                            for (k in 0 until numClasses) {
                                placeTransitTokenTime[placeListIdx][k] += placeTokensInTransit[placeListIdx][k] * placeTransitElapsed
                            }
                            lastPlaceTransitUpdateTime[placeListIdx] = currentTime
                        }
                        // Remove these tokens from place transit tracking
                        placeTokensInTransit[placeListIdx][classIdx] -= required
                    }
                }
            }

            // Decrement in-service count
            transitionInService[transListIdx][modeIdx]--

            // Check event count for stopping/warmup/MSER sampling
            checkEventCountStop()

            // Produce output tokens
            for (nodeIdx in 0 until numNodes) {
                for (classIdx in 0 until numClasses) {
                    val produced = mode.firingOutcomes[nodeIdx][classIdx]
                    if (produced > 0) {
                        // Route produced tokens to destination node
                        for (t in 0 until produced) {
                            routeTokenFromTransition(transNodeIdx, nodeIdx, classIdx)
                        }
                    }
                }
            }

            // Check for more enabled transitions
            checkAndFireTransitions()
        }
    }

    /**
     * Route a token from a transition to a destination node.
     */
    private fun routeTokenFromTransition(transNodeIdx: Int, destNodeIdx: Int, classId: Int) {
        val currentTime = Sim.time()

        when {
            placeNodes.contains(destNodeIdx) -> {
                // Destination is a Place - add token with proper statistics tracking
                val placeListIdx = placeNodes.indexOf(destNodeIdx)
                if (placeListIdx >= 0) {
                    // Update time-weighted token count before adding new token
                    val elapsed = currentTime - lastPlaceUpdateTime[placeListIdx]
                    if (elapsed > 0) {
                        for (k in 0 until numClasses) {
                            totalPlaceTokenTime[placeListIdx][k] += placeTokens[placeListIdx][k] * elapsed
                        }
                        lastPlaceUpdateTime[placeListIdx] = currentTime
                    }
                    // Add token to place
                    placeTokens[placeListIdx][classId]++
                }
                // Note: checkAndFireTransitions() is called at the end of TransitionFiring
            }
            sinkNodes.contains(destNodeIdx) -> {
                // Token leaves system
                systemCompletedCustomers[classId]++
            }
            serviceNodes.contains(destNodeIdx) -> {
                // Route to service node (queue)
                val queueIdx = serviceNodes.indexOf(destNodeIdx)
                val customer = Customer(classId, classPrio[classId], currentTime, currentTime, siroRng.nextDouble(),
                    absoluteDeadline = currentTime + classDeadline[classId])
                arriveAtQueue(queueIdx, customer)
            }
            forkNodes.contains(destNodeIdx) -> {
                // Route to fork
                val parentJobId = nextJobId++
                handleForkArrival(destNodeIdx, parentJobId, classId, currentTime)
            }
            transitionNodes.contains(destNodeIdx) -> {
                // Token goes directly to another transition (unusual but handle it)
                // This would require the transition to check enabling again
                checkAndFireTransitions()
            }
        }
    }

    /**
     * Check if a node is a Place node.
     */
    private fun isPlaceNode(nodeIdx: Int): Boolean = placeNodes.contains(nodeIdx)

    /**
     * Check if a node is a Transition node.
     */
    private fun isTransitionNode(nodeIdx: Int): Boolean = transitionNodes.contains(nodeIdx)

    // Result getters
    fun getAvgQueueLength(queueIdx: Int, classId: Int): Double {
        // Use MSER-5 truncated observations if available
        if (mserEnabled && ::queueLengthObservations.isInitialized &&
            queueIdx < queueLengthObservations.size &&
            classId < queueLengthObservations[queueIdx].size) {
            val obs = queueLengthObservations[queueIdx][classId]
            val truncationIdx = mserTruncationBatch * effectiveMserBatchSize
            if (truncationIdx < obs.size) {
                var sum = 0.0
                var count = 0
                for (i in truncationIdx until obs.size) {
                    sum += obs[i]
                    count++
                }
                if (count > 0) return sum / count
            }
        }
        // Fallback to time-weighted average
        // totalQueueTime now includes blocked servers (via updateQueueStats with currentBlockedServers)
        val simTime = getActualSimTime()
        return if (simTime > 0) {
            totalQueueTime[queueIdx][classId] / simTime
        } else {
            0.0
        }
    }

    /**
     * Gets the number of servers compatible with a given class at a queue.
     * For heterogeneous queues, returns the sum of servers across all compatible types.
     * For homogeneous queues, returns the total number of servers.
     */
    private fun getCompatibleServerCount(queueIdx: Int, classId: Int): Int {
        if (numServerTypes[queueIdx] > 0) {
            var count = 0
            for (typeId in 0 until numServerTypes[queueIdx]) {
                if (serverCompat[queueIdx][typeId][classId]) {
                    count += serversPerType[queueIdx][typeId]
                }
            }
            return if (count > 0) count else numServers[queueIdx]
        }
        return numServers[queueIdx]
    }

    /**
     * Get effective utilization (service time only, excluding blocking time).
     * This represents the fraction of time servers are actually performing service work.
     * For heterogeneous queues, utilization is calculated relative to compatible servers only.
     */
    fun getEffectiveUtilization(queueIdx: Int, classId: Int): Double {
        if (isDelayNode[queueIdx]) {
            // For Delay nodes, effective utilization = service time contribution to queue length
            val simTime = getActualSimTime()
            return if (simTime > 0) {
                totalQueueTime[queueIdx][classId] / simTime
            } else {
                0.0
            }
        }
        // Per-server Utilization = busy time / (simulation time * compatibleServers) (λ × S / c convention, consistent with MVA)
        val simTime = getActualSimTime()
        val compatibleServers = getCompatibleServerCount(queueIdx, classId)
        return if (simTime > 0 && compatibleServers > 0) {
            totalBusyTime[queueIdx][classId] / (simTime * compatibleServers)
        } else {
            0.0
        }
    }

    fun getUtilization(queueIdx: Int, classId: Int): Double {
        // For Delay nodes (infinite servers), utilization is defined as mean queue length
        // (i.e., traffic intensity λ/μ, which equals the average number of busy servers)
        // Include blocking time for REPLY signal semantics
        if (isDelayNode[queueIdx]) {
            return getAvgQueueLength(queueIdx, classId)
        }
        // Per-server Utilization = (busy time + blocking time) / (simulation time * compatibleServers)
        // For heterogeneous queues, utilization is calculated relative to compatible servers only.
        // This gives the per-server utilization (λ × S / c convention, consistent with MVA and LQNS)
        val simTime = getActualSimTime()
        val compatibleServers = getCompatibleServerCount(queueIdx, classId)
        return if (simTime > 0 && compatibleServers > 0) {
            (totalBusyTime[queueIdx][classId] + totalBlockingTime[queueIdx][classId]) / (simTime * compatibleServers)
        } else {
            0.0
        }
    }

    fun getAvgResponseTime(queueIdx: Int, classId: Int): Double {
        return if (responseTimeTally[queueIdx][classId].numberObs() > 0) {
            responseTimeTally[queueIdx][classId].average()
        } else {
            // Use M/M/1 formula as fallback
            val mu = mus[queueIdx][classId]
            val lambda = getThroughput(queueIdx, classId)
            if (mu > lambda && mu < Double.MAX_VALUE) 1.0 / (mu - lambda) else 0.0
        }
    }

    fun getThroughput(queueIdx: Int, classId: Int): Double {
        // Use MSER-5 truncated throughput if available
        if (mserEnabled && ::throughputObservations.isInitialized && ::observationTimes.isInitialized &&
            queueIdx < throughputObservations.size &&
            classId < throughputObservations[queueIdx].size) {
            val tputObs = throughputObservations[queueIdx][classId]
            val truncationIdx = mserTruncationBatch * effectiveMserBatchSize
            if (truncationIdx < tputObs.size && truncationIdx < observationTimes.size) {
                val startCompletions = tputObs[truncationIdx]
                val endCompletions = completedCustomers[queueIdx][classId]
                val startTime = observationTimes[truncationIdx]
                val endTime = Sim.time()
                val elapsed = endTime - startTime
                if (elapsed > 0) {
                    return (endCompletions - startCompletions).toDouble() / elapsed
                }
            }
        }
        // Fallback
        val simTime = getActualSimTime()
        return if (simTime > 0) {
            completedCustomers[queueIdx][classId].toDouble() / simTime
        } else {
            0.0
        }
    }

    /**
     * Returns the arrival rate at a service node for a given class.
     * This includes all arrivals (including those that were dropped).
     */
    fun getArrivalRate(queueIdx: Int, classId: Int): Double {
        val simTime = getActualSimTime()
        return if (simTime > 0) {
            arrivedCustomers[queueIdx][classId].toDouble() / simTime
        } else {
            0.0
        }
    }

    fun getSystemThroughput(classId: Int): Double {
        // Use MSER-5 truncated time
        val simTime = getActualSimTime()
        return if (simTime > 0) {
            systemCompletedCustomers[classId].toDouble() / simTime
        } else {
            0.0
        }
    }

    fun getSystemResponseTime(classId: Int): Double {
        return if (systemResponseTimeTally[classId].numberObs() > 0) {
            systemResponseTimeTally[classId].average()
        } else {
            0.0
        }
    }

    // ==================== OBM Confidence Interval Functions ====================

    /**
     * Data class to hold confidence interval results for all metrics.
     */
    private data class CIResults(
        val QNCI: Matrix,
        val UNCI: Matrix,
        val RNCI: Matrix,
        val TNCI: Matrix,
        val ANCI: Matrix,
        val WNCI: Matrix
    )

    /**
     * Compute overlap adjustment factor for OBM variance estimation.
     * Based on Alexopoulos & Seila (1998) for 50% overlap.
     * For other overlap fractions, uses linear interpolation.
     * @param overlap Overlap fraction (0.0 to 1.0)
     * @return Variance adjustment factor
     */
    private fun computeOverlapAdjustmentFactor(overlap: Double): Double {
        return when {
            overlap <= 0.0 -> 1.0  // Non-overlapping
            overlap >= 0.5 -> 4.0 / 3.0  // 50% overlap (standard OBM)
            else -> 1.0 + (overlap / 0.5) * (4.0 / 3.0 - 1.0)  // Linear interpolation
        }
    }

    /**
     * Compute OBM (Overlapping Batch Means) statistics with configurable overlap.
     * @param observations List of observations to analyze
     * @param batchSize Size of each batch
     * @return Triple(grandMean, stdError, effectiveDf) or null if insufficient data
     */
    private fun computeOBMStatisticsInternal(observations: List<Double>, batchSize: Int): Triple<Double, Double, Int>? {
        val n = observations.size
        if (n < batchSize * 2) return null

        // Use configurable overlap fraction
        val stepSize = maxOf(1, (batchSize * (1.0 - effectiveObmOverlap)).toInt())
        val numBatches = (n - batchSize) / stepSize + 1
        if (numBatches < 2) return null

        // Compute overlapping batch means
        val batchMeans = DoubleArray(numBatches)
        for (i in 0 until numBatches) {
            val startIdx = i * stepSize
            var sum = 0.0
            for (j in 0 until batchSize) {
                if (startIdx + j < n) {
                    sum += observations[startIdx + j]
                }
            }
            batchMeans[i] = sum / batchSize
        }

        // Compute grand mean
        var grandMean = 0.0
        for (i in 0 until numBatches) {
            grandMean += batchMeans[i]
        }
        grandMean /= numBatches

        // Compute sum of squared differences
        var sumSquaredDiff = 0.0
        for (i in 0 until numBatches) {
            val diff = batchMeans[i] - grandMean
            sumSquaredDiff += diff * diff
        }

        // OBM variance with configurable overlap adjustment factor
        val overlapAdjustment = computeOverlapAdjustmentFactor(effectiveObmOverlap)
        val batchMeanVariance = overlapAdjustment * sumSquaredDiff / (numBatches - 1)
        val stdError = kotlin.math.sqrt(batchMeanVariance / numBatches)

        // Effective degrees of freedom for OBM
        val effectiveDf = maxOf(1, ((numBatches - 1) / overlapAdjustment).toInt())

        return Triple(grandMean, stdError, effectiveDf)
    }

    /**
     * Compute BM (Non-overlapping Batch Means) statistics.
     * @param observations List of observations to analyze
     * @param batchSize Size of each batch
     * @return Triple(grandMean, stdError, df) or null if insufficient data
     */
    private fun computeBMStatisticsInternal(observations: List<Double>, batchSize: Int): Triple<Double, Double, Int>? {
        val n = observations.size
        val numBatches = n / batchSize
        if (numBatches < 2) return null

        // Compute non-overlapping batch means
        val batchMeans = DoubleArray(numBatches)
        for (i in 0 until numBatches) {
            val startIdx = i * batchSize
            var sum = 0.0
            for (j in 0 until batchSize) {
                sum += observations[startIdx + j]
            }
            batchMeans[i] = sum / batchSize
        }

        // Compute grand mean
        var grandMean = 0.0
        for (bm in batchMeans) grandMean += bm
        grandMean /= numBatches

        // Sample variance of batch means (no overlap adjustment needed)
        var sumSquaredDiff = 0.0
        for (i in 0 until numBatches) {
            val diff = batchMeans[i] - grandMean
            sumSquaredDiff += diff * diff
        }
        val batchMeanVariance = sumSquaredDiff / (numBatches - 1)
        val stdError = kotlin.math.sqrt(batchMeanVariance / numBatches)

        val df = numBatches - 1

        return Triple(grandMean, stdError, df)
    }

    /**
     * Compute batch means statistics using the configured CI method.
     * @param observations List of observations to analyze
     * @param batchSize Size of each batch
     * @return Triple(grandMean, stdError, df) or null if insufficient data or CI disabled
     */
    private fun computeCIStatistics(observations: List<Double>, batchSize: Int): Triple<Double, Double, Int>? {
        return when (effectiveCiMethod) {
            "obm" -> computeOBMStatisticsInternal(observations, batchSize)
            "bm" -> computeBMStatisticsInternal(observations, batchSize)
            else -> null  // "none" - no CI computation
        }
    }

    /**
     * Get t-distribution critical value for given confidence level and degrees of freedom.
     * @param confintLevel Confidence level (e.g., 0.95 for 95% CI)
     * @param df Degrees of freedom
     * @return Critical value from t-distribution
     */
    private fun getTCriticalValueInternal(confintLevel: Double, df: Int): Double {
        // t-distribution critical values for common confidence levels
        // Table lookup for small df (up to 30)
        val tTable95 = doubleArrayOf(
            12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
            2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
            2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042
        )
        val tTable90 = doubleArrayOf(
            6.314, 2.920, 2.353, 2.132, 2.015, 1.943, 1.895, 1.860, 1.833, 1.812,
            1.796, 1.782, 1.771, 1.761, 1.753, 1.746, 1.740, 1.734, 1.729, 1.725,
            1.721, 1.717, 1.714, 1.711, 1.708, 1.706, 1.703, 1.701, 1.699, 1.697
        )
        val tTable99 = doubleArrayOf(
            63.657, 9.925, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.250, 3.169,
            3.106, 3.055, 3.012, 2.977, 2.947, 2.921, 2.898, 2.878, 2.861, 2.845,
            2.831, 2.819, 2.807, 2.797, 2.787, 2.779, 2.771, 2.763, 2.756, 2.750
        )

        // Table index is df - 1 (df=1 is at index 0)
        val tableIdx = minOf(df, 30) - 1
        if (tableIdx < 0) return 1.96 // Fallback to z-value

        // Select table based on confidence level (checking higher confidence first)
        return when {
            confintLevel >= 0.99 -> if (tableIdx < tTable99.size) tTable99[tableIdx] else 2.576
            confintLevel >= 0.95 -> if (tableIdx < tTable95.size) tTable95[tableIdx] else 1.96
            confintLevel >= 0.90 -> if (tableIdx < tTable90.size) tTable90[tableIdx] else 1.645
            else -> 1.96 // Default to 95% CI z-value
        }
    }

    /**
     * Compute confidence intervals for all metrics using configured CI method.
     * Uses MSER truncation point and applies batch means to post-warmup observations.
     * @param confintLevel Confidence level (e.g., 0.95)
     * @return CIResults containing CI half-widths for all metrics
     */
    private fun computeOBMConfidenceIntervals(confintLevel: Double): CIResults {
        val QNCI = Matrix(numStations, numClasses)
        QNCI.fill(0.0)
        val UNCI = Matrix(numStations, numClasses)
        UNCI.fill(0.0)
        val RNCI = Matrix(numStations, numClasses)
        RNCI.fill(0.0)
        val TNCI = Matrix(numStations, numClasses)
        TNCI.fill(0.0)
        val ANCI = Matrix(numStations, numClasses)
        ANCI.fill(0.0)
        val WNCI = Matrix(numStations, numClasses)
        WNCI.fill(0.0)

        // Skip CI computation if disabled
        if (effectiveCiMethod == "none") {
            return CIResults(QNCI, UNCI, RNCI, TNCI, ANCI, WNCI)
        }

        // Get truncation index from MSER
        val truncationIdx = mserTruncationBatch * effectiveMserBatchSize

        // Compute CI for each service station and class
        for ((svcIdx, serviceStation) in serviceStations.withIndex()) {
            for (k in 0 until numClasses) {
                // Queue length CI from observations
                if (::queueLengthObservations.isInitialized &&
                    svcIdx < queueLengthObservations.size &&
                    k < queueLengthObservations[svcIdx].size) {
                    val allObs = queueLengthObservations[svcIdx][k]
                    if (truncationIdx < allObs.size) {
                        val postWarmupObs = allObs.subList(truncationIdx, allObs.size)
                        if (postWarmupObs.size >= effectiveCiMinObs) {
                            val batchSize = maxOf(effectiveCiMinBatch, kotlin.math.sqrt(postWarmupObs.size.toDouble()).toInt())
                            val stats = computeCIStatistics(postWarmupObs, batchSize)
                            if (stats != null) {
                                val tCrit = getTCriticalValueInternal(confintLevel, stats.third)
                                QNCI[serviceStation, k] = tCrit * stats.second
                            }
                        }
                    }
                }

                // Utilization CI - derive from throughput and service rate
                if (::throughputObservations.isInitialized &&
                    svcIdx < throughputObservations.size &&
                    k < throughputObservations[svcIdx].size) {
                    val compObs = throughputObservations[svcIdx][k]
                    if (truncationIdx < compObs.size && truncationIdx < observationTimes.size) {
                        // Convert cumulative completions to per-interval rates
                        val postWarmupComp = compObs.subList(truncationIdx, compObs.size)
                        val postWarmupTimes = observationTimes.subList(truncationIdx, observationTimes.size)
                        if (postWarmupComp.size >= 2) {
                            val rates = ArrayList<Double>()
                            for (i in 1 until postWarmupComp.size) {
                                val deltaComp = postWarmupComp[i] - postWarmupComp[i - 1]
                                val deltaTime = postWarmupTimes[i] - postWarmupTimes[i - 1]
                                if (deltaTime > 0) {
                                    rates.add(deltaComp.toDouble() / deltaTime)
                                }
                            }
                            if (rates.size >= effectiveCiMinObs) {
                                val batchSize = maxOf(effectiveCiMinBatch, kotlin.math.sqrt(rates.size.toDouble()).toInt())
                                val stats = computeCIStatistics(rates, batchSize)
                                if (stats != null) {
                                    val tCrit = getTCriticalValueInternal(confintLevel, stats.third)
                                    TNCI[serviceStation, k] = tCrit * stats.second
                                    // Utilization CI: convert throughput CI using service rate
                                    val mu = mus[svcIdx][k]
                                    if (mu > 0 && mu < Double.MAX_VALUE) {
                                        val numServers = sn.nservers[serviceStations[svcIdx]].toInt()
                                        UNCI[serviceStation, k] = (tCrit * stats.second) / (mu * numServers)
                                    }
                                }
                            }
                        }
                    }
                }

                // Response time CI from tally variance
                if (svcIdx < responseTimeTally.size && k < responseTimeTally[svcIdx].size) {
                    val tally = responseTimeTally[svcIdx][k]
                    if (tally.numberObs() > 30) {
                        val stdErr = tally.standardDeviation() / kotlin.math.sqrt(tally.numberObs().toDouble())
                        val tCrit = getTCriticalValueInternal(confintLevel, (tally.numberObs() - 1).toInt())
                        RNCI[serviceStation, k] = tCrit * stdErr
                    }
                }
            }
        }

        // WNCI = RNCI (residence time same as response time for single-visit)
        for (i in 0 until numStations) {
            for (k in 0 until numClasses) {
                WNCI[i, k] = RNCI[i, k]
            }
        }

        return CIResults(QNCI, UNCI, RNCI, TNCI, ANCI, WNCI)
    }
}

// ==================== Top-Level OBM Test Helper Functions ====================

/**
 * Compute OBM (Overlapping Batch Means) statistics with 50% overlap.
 * Exposed as top-level function for testing.
 * @param observations List of observations to analyze
 * @param batchSize Size of each batch
 * @return Triple(grandMean, stdError, effectiveDf) or null if insufficient data
 */
fun computeOBMStatistics(observations: List<Double>, batchSize: Int): Triple<Double, Double, Int>? {
    val n = observations.size
    if (n < batchSize * 2) return null

    val halfBatch = batchSize / 2
    val numBatches = (n - batchSize) / halfBatch + 1
    if (numBatches < 2) return null

    // Compute overlapping batch means
    val batchMeans = DoubleArray(numBatches)
    for (i in 0 until numBatches) {
        val startIdx = i * halfBatch
        var sum = 0.0
        for (j in 0 until batchSize) {
            sum += observations[startIdx + j]
        }
        batchMeans[i] = sum / batchSize
    }

    // Compute grand mean
    var grandMean = 0.0
    for (i in 0 until numBatches) {
        grandMean += batchMeans[i]
    }
    grandMean /= numBatches

    // Compute sum of squared differences
    var sumSquaredDiff = 0.0
    for (i in 0 until numBatches) {
        val diff = batchMeans[i] - grandMean
        sumSquaredDiff += diff * diff
    }

    // OBM variance with 4/3 overlap adjustment factor
    val overlapAdjustment = 4.0 / 3.0
    val batchMeanVariance = overlapAdjustment * sumSquaredDiff / (numBatches - 1)
    val stdError = kotlin.math.sqrt(batchMeanVariance / numBatches)

    // Effective degrees of freedom for OBM
    val effectiveDf = maxOf(1, ((numBatches - 1) / overlapAdjustment).toInt())

    return Triple(grandMean, stdError, effectiveDf)
}

/**
 * Get t-distribution critical value for given confidence level and degrees of freedom.
 * Exposed as top-level function for testing.
 * @param confintLevel Confidence level (e.g., 0.95 for 95% CI)
 * @param df Degrees of freedom
 * @return Critical value from t-distribution
 */
fun getTCriticalValue(confintLevel: Double, df: Int): Double {
    // Lookup based on confidence level directly for better precision
    val tTable95 = doubleArrayOf(
        12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
        2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
        2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042
    )
    val tTable90 = doubleArrayOf(
        6.314, 2.920, 2.353, 2.132, 2.015, 1.943, 1.895, 1.860, 1.833, 1.812,
        1.796, 1.782, 1.771, 1.761, 1.753, 1.746, 1.740, 1.734, 1.729, 1.725,
        1.721, 1.717, 1.714, 1.711, 1.708, 1.706, 1.703, 1.701, 1.699, 1.697
    )
    val tTable99 = doubleArrayOf(
        63.657, 9.925, 5.841, 4.604, 4.032, 3.707, 3.499, 3.355, 3.250, 3.169,
        3.106, 3.055, 3.012, 2.977, 2.947, 2.921, 2.898, 2.878, 2.861, 2.845,
        2.831, 2.819, 2.807, 2.797, 2.787, 2.779, 2.771, 2.763, 2.756, 2.750
    )

    // Table index is df - 1 (df=1 is at index 0)
    val tableIdx = minOf(df, 30) - 1
    if (tableIdx < 0) return 1.96

    // Select table based on confidence level (checking higher confidence first)
    return when {
        confintLevel >= 0.99 -> if (tableIdx < tTable99.size) tTable99[tableIdx] else 2.576
        confintLevel >= 0.95 -> if (tableIdx < tTable95.size) tTable95[tableIdx] else 1.96
        confintLevel >= 0.90 -> if (tableIdx < tTable90.size) tTable90[tableIdx] else 1.645
        else -> 1.96
    }
}

/**
 * Compute standard (non-overlapping) batch means statistics.
 * Exposed as top-level function for testing comparison with OBM.
 * @param observations List of observations to analyze
 * @param batchSize Size of each batch
 * @return Triple(grandMean, stdError, effectiveDf) or null if insufficient data
 */
fun computeStandardBatchMeansStatistics(observations: List<Double>, batchSize: Int): Triple<Double, Double, Int>? {
    val n = observations.size
    val numBatches = n / batchSize
    if (numBatches < 2) return null

    // Compute non-overlapping batch means
    val batchMeans = DoubleArray(numBatches)
    for (i in 0 until numBatches) {
        val startIdx = i * batchSize
        var sum = 0.0
        for (j in 0 until batchSize) {
            sum += observations[startIdx + j]
        }
        batchMeans[i] = sum / batchSize
    }

    // Compute grand mean
    var grandMean = 0.0
    for (i in 0 until numBatches) {
        grandMean += batchMeans[i]
    }
    grandMean /= numBatches

    // Compute sum of squared differences
    var sumSquaredDiff = 0.0
    for (i in 0 until numBatches) {
        val diff = batchMeans[i] - grandMean
        sumSquaredDiff += diff * diff
    }

    // Standard batch means variance (no overlap adjustment)
    val batchMeanVariance = sumSquaredDiff / (numBatches - 1)
    val stdError = kotlin.math.sqrt(batchMeanVariance / numBatches)

    return Triple(grandMean, stdError, numBatches - 1)
}
