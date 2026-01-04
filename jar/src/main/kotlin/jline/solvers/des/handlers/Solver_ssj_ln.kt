/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * @file Solver_ssj_ln.kt
 * @brief SSJ-based discrete event simulation engine for LayeredNetwork (LQN) models.
 *
 * @details This file contains the DES simulation engine for Layered Queueing Networks
 * built on the SSJ (Stochastic Simulation in Java) library. It provides steady-state
 * analysis for LQN models with tasks, entries, activities, and synchronous/asynchronous calls.
 *
 * @section ssj_ln_architecture Architecture Overview
 *
 * The simulation engine is organized around the LNSSJSimulator class which manages:
 * - Task lifecycle (creation, service, blocking, completion)
 * - Entry call handling (synchronous blocking, asynchronous fire-and-forget)
 * - Activity execution with precedence constraints (SEQ, AND-fork/join, OR-fork/join)
 * - Host/processor multiplicity and scheduling
 * - Statistics collection per host, task, entry, activity
 *
 * @author QORE Lab, Imperial College London
 * @since 1.0
 */

package jline.solvers.des.handlers

import jline.lang.constant.ActivityPrecedenceType
import jline.lang.constant.CallType
import jline.lang.constant.SchedStrategy
import jline.lang.layered.LayeredNetworkElement
import jline.lang.layered.LayeredNetworkStruct
import jline.lang.processes.Distribution
import jline.solvers.SolverOptions
import jline.solvers.des.LNDESResult
import jline.util.matrix.Matrix
import umontreal.ssj.simevents.Event
import umontreal.ssj.simevents.Sim
import umontreal.ssj.stat.Tally
import java.util.LinkedList
import java.util.PriorityQueue
import java.util.Random

/**
 * Steady-state LayeredNetwork simulation using SSJ library.
 *
 * @param lsn LayeredNetwork structure
 * @param options Solver configuration options
 * @return LNDESResult with steady-state performance metrics
 */
fun solver_ssj_ln(lsn: LayeredNetworkStruct, options: SolverOptions): LNDESResult {
    val simulator = LNSSJSimulator(lsn, options)
    val maxEvents = options.samples.toDouble()
    simulator.simulate(maxEvents)
    return simulator.getLNDESResult()
}

/**
 * Core discrete event simulation engine for LayeredNetwork models using SSJ library.
 */
internal class LNSSJSimulator(
    private val lsn: LayeredNetworkStruct,
    private val options: SolverOptions
) {
    // LQN dimensions
    private val nhosts = lsn.nhosts
    private val ntasks = lsn.ntasks
    private val nentries = lsn.nentries
    private val nacts = lsn.nacts
    private val ncalls = lsn.ncalls

    // Index shifts for navigating lsn indices
    private val hshift = lsn.hshift
    private val tshift = lsn.tshift
    private val eshift = lsn.eshift
    private val ashift = lsn.ashift
    private val cshift = lsn.cshift

    // Random number generator
    private val seed = options.seed.toLong()
    private val random = Random(seed)

    // Host (processor) state
    private val hostBusy = IntArray(nhosts)         // Number of busy servers per host
    private val hostServers = IntArray(nhosts)      // Number of servers (multiplicity) per host
    private lateinit var hostSchedStrategy: Array<SchedStrategy>  // Scheduling strategy per host
    private lateinit var hostQueues: Array<PriorityQueue<HostRequest>>  // Priority queues for hosts

    // Task state
    private val taskBusy = IntArray(ntasks)         // Active task instances
    private val taskMultiplicity = IntArray(ntasks) // Max concurrent instances
    private lateinit var taskSchedStrategy: Array<SchedStrategy>  // Scheduling strategy per task
    private lateinit var taskQueues: Array<PriorityQueue<LNRequest>>  // Priority queues for tasks

    // PS (Processor Sharing) state for hosts
    private lateinit var hostPSJobs: Array<MutableList<PSJob>>  // Active PS jobs per host
    private lateinit var hostPSEvent: Array<PSUpdateEvent?>     // Scheduled PS update events per host

    // PS time slice interval
    private val psTimeSlice = 0.001  // 1ms time slices for PS scheduling

    // Request tracking
    private var nextRequestId: Long = 0L
    private val activeRequests = mutableMapOf<Long, LNRequest>()

    // AND-fork/join tracking
    private val pendingAndJoins = mutableMapOf<Long, AndJoinState>()

    // Sync call tracking: maps caller request ID to queue of pending sync calls (target entry indices)
    // Sync calls are executed sequentially - one at a time
    private val pendingSyncCallQueue = mutableMapOf<Long, LinkedList<Int>>()
    // Maps caller request ID to the activity they were executing when call was made
    private val syncCallActivity = mutableMapOf<Long, Int>()
    // Maps caller request ID to indicate a sync call is currently in flight
    private val syncCallInFlight = mutableMapOf<Long, Boolean>()

    // Statistics per host (indexed 0 to nhosts-1)
    private val hostBusyTime = DoubleArray(nhosts)
    private val hostCompletions = IntArray(nhosts)
    private val lastHostUpdateTime = DoubleArray(nhosts)

    // Statistics per task (indexed 0 to ntasks-1)
    private val taskUtilTime = DoubleArray(ntasks)
    private val taskCompletions = DoubleArray(ntasks)
    private val lastTaskUpdateTime = DoubleArray(ntasks)

    // Statistics per entry (indexed 0 to nentries-1)
    private val entryThroughput = DoubleArray(nentries)
    private val entryResponseTimeTally = Array(nentries) { Tally() }

    // Statistics per activity (indexed 0 to nacts-1)
    private val activityThroughput = DoubleArray(nacts)
    private val activityServiceTimeTally = Array(nacts) { Tally() }
    private val activityBusyTime = DoubleArray(nacts)
    private val lastActivityUpdateTime = DoubleArray(nacts)
    private val activityBusy = IntArray(nacts)  // Count of active instances per activity

    // Track activity start times per request: maps (requestId, activityIdx) -> start time
    private val activityStartTimes = mutableMapOf<Pair<Long, Int>, Double>()

    // Queue length statistics for tasks (time-weighted queue lengths)
    // taskQueueLength[t] = current number of jobs at task t (busy + waiting)
    // taskQueueLengthTime[t] = cumulative time-weighted queue length for task t
    // lastTaskQueueUpdateTime[t] = last time the queue length was updated
    private val taskQueueLength = IntArray(ntasks)
    private val taskQueueLengthTime = DoubleArray(ntasks)
    private val lastTaskQueueUpdateTime = DoubleArray(ntasks)

    // Queue length statistics for hosts (time-weighted queue lengths)
    private val hostQueueLength = IntArray(nhosts)
    private val hostQueueLengthTime = DoubleArray(nhosts)
    private val lastHostQueueUpdateTime = DoubleArray(nhosts)

    // Event counting for stopping criterion
    private var totalEventCount: Long = 0
    private var maxEvents: Long = 0
    private var warmupEventThreshold: Long = 0
    private var warmupDone = false
    private var warmupEndTime = 0.0

    // Mapping from entry index to its bound activity
    private val entryBoundActivity = mutableMapOf<Int, Int>()

    // Mapping from entry index to its reply activity (if any)
    private val entryReplyActivity = mutableMapOf<Int, Int>()

    // Mapping from activity index to parent entry index
    private val activityToEntry = mutableMapOf<Int, Int>()

    // Activity graph successors: aidx -> list of successor activity indices
    private val activitySuccessors = mutableMapOf<Int, MutableList<Int>>()

    // OR-fork probabilities: aidx -> list of (successor_aidx, probability)
    private val orForkProbs = mutableMapOf<Int, MutableList<Pair<Int, Double>>>()

    // Comparators for scheduling strategies
    private val fcfsComparator = Comparator<LNRequest> { r1, r2 ->
        r1.taskArrivalTime.compareTo(r2.taskArrivalTime)
    }

    private val lcfsComparator = Comparator<LNRequest> { r1, r2 ->
        r2.taskArrivalTime.compareTo(r1.taskArrivalTime)  // Later arrival first
    }

    private val hostFcfsComparator = Comparator<HostRequest> { h1, h2 ->
        h1.queueArrivalTime.compareTo(h2.queueArrivalTime)
    }

    private val hostLcfsComparator = Comparator<HostRequest> { h1, h2 ->
        h2.queueArrivalTime.compareTo(h1.queueArrivalTime)  // Later arrival first
    }

    /**
     * Data class for PS (Processor Sharing) jobs at a host.
     */
    data class PSJob(
        val hostRequest: HostRequest,
        var remainingService: Double,   // Remaining service time
        var lastUpdateTime: Double      // When service was last updated
    )

    /**
     * Data class representing a request flowing through the LQN.
     */
    data class LNRequest(
        val requestId: Long,
        val entryIdx: Int,              // Target entry being called
        val callType: CallType,         // SYNC or ASYNC
        val arrivalTime: Double,        // When request entered system
        var taskArrivalTime: Double,    // When request arrived at current task
        val callerRequestId: Long,      // Parent request waiting for reply (-1 if none)
        val callerTaskIdx: Int,         // Task of the caller (-1 if REF task origin)
        var currentActivityIdx: Int,    // Current activity being executed
        var parentForkId: Long          // For AND-fork tracking (-1 if not in fork)
    )

    /**
     * Data class for host service requests.
     */
    data class HostRequest(
        val request: LNRequest,
        val activityIdx: Int,
        val serviceTime: Double,
        val queueArrivalTime: Double
    )

    /**
     * State tracking for AND-join synchronization.
     */
    data class AndJoinState(
        val parentRequest: LNRequest,
        val joinActivityIdx: Int,
        var pendingBranches: Int,
        val branchCompletionTimes: MutableList<Double>
    )

    init {
        initialize()
    }

    /**
     * Initialize simulation state from LQN structure.
     */
    private fun initialize() {
        // Initialize host scheduling strategies and multiplicity
        hostSchedStrategy = Array(nhosts) { SchedStrategy.FCFS }
        hostPSJobs = Array(nhosts) { mutableListOf<PSJob>() }
        hostPSEvent = arrayOfNulls(nhosts)

        for (h in 0 until nhosts) {
            val hidx = h + 1  // Host indices start at 1
            hostServers[h] = lsn.mult.get(0, hidx).toInt().coerceAtLeast(1)
            hostSchedStrategy[h] = lsn.sched[hidx] ?: SchedStrategy.FCFS
        }

        // Initialize host queues with appropriate comparators
        hostQueues = Array(nhosts) { h ->
            PriorityQueue(11, getHostComparator(hostSchedStrategy[h]))
        }

        // Initialize task scheduling strategies and multiplicity
        taskSchedStrategy = Array(ntasks) { SchedStrategy.FCFS }

        for (t in 0 until ntasks) {
            val tidx = tshift + t + 1
            taskMultiplicity[t] = lsn.mult.get(0, tidx).toInt().coerceAtLeast(1)
            val sched = lsn.sched[tidx] ?: SchedStrategy.FCFS
            // REF tasks use think time, not scheduling
            taskSchedStrategy[t] = if (sched == SchedStrategy.REF) SchedStrategy.FCFS else sched
        }

        // Initialize task queues with appropriate comparators
        taskQueues = Array(ntasks) { t ->
            PriorityQueue(11, getTaskComparator(taskSchedStrategy[t]))
        }

        // Build entry to bound activity mapping
        buildEntryActivityMappings()

        // Build activity graph
        buildActivityGraph()
    }

    /**
     * Get comparator for host queue based on scheduling strategy.
     */
    private fun getHostComparator(strategy: SchedStrategy): Comparator<HostRequest> {
        return when (strategy) {
            SchedStrategy.LCFS -> hostLcfsComparator
            SchedStrategy.PS, SchedStrategy.INF -> hostFcfsComparator  // PS handled separately
            else -> hostFcfsComparator  // FCFS is default
        }
    }

    /**
     * Get comparator for task queue based on scheduling strategy.
     */
    private fun getTaskComparator(strategy: SchedStrategy): Comparator<LNRequest> {
        return when (strategy) {
            SchedStrategy.LCFS -> lcfsComparator
            SchedStrategy.PS, SchedStrategy.INF -> fcfsComparator  // PS/INF handled separately
            else -> fcfsComparator  // FCFS is default
        }
    }

    /**
     * Check if a scheduling strategy uses Processor Sharing.
     */
    private fun isPSScheduling(strategy: SchedStrategy): Boolean {
        return strategy == SchedStrategy.PS || strategy == SchedStrategy.DPS || strategy == SchedStrategy.GPS
    }

    /**
     * Check if a scheduling strategy is INF (infinite server / delay).
     */
    private fun isINFScheduling(strategy: SchedStrategy): Boolean {
        return strategy == SchedStrategy.INF
    }

    /**
     * Build mappings from entries to their bound and reply activities.
     */
    private fun buildEntryActivityMappings() {
        // For each entry, find its bound activity by checking replygraph
        // Note: replygraph is indexed by raw activity index (1..nacts) and entry index (1..nentries),
        // not by global shifted indices
        for (e in 1..nentries) {
            val eidx = eshift + e  // Global index for lookups in parent, names, etc.
            val tidx = lsn.parent.get(0, eidx).toInt()
            val activities = lsn.actsof[tidx] ?: continue

            for (aidx in activities) {
                // Convert global activity index to raw activity index for replygraph
                val rawAidx = aidx - ashift
                // Check if this activity is bound to this entry
                // The replygraph matrix indicates reply relationships
                // replygraph is indexed (activity: 1..nacts, entry: 1..nentries)
                if (lsn.replygraph != null && rawAidx >= 0 && rawAidx < lsn.replygraph.numRows &&
                    e >= 0 && e < lsn.replygraph.numCols) {
                    if (lsn.replygraph.get(rawAidx, e) > 0) {
                        entryReplyActivity[eidx] = aidx
                    }
                }
            }

            // The first activity in the task that has this entry as parent is the bound activity
            // We need to find the activity that starts when this entry is called
            // This is typically indicated by the graph structure
            if (activities.isNotEmpty()) {
                // Find the activity with no predecessors within this entry's scope
                // For simplicity, assume the first activity listed is the bound one
                val firstActivity = findBoundActivityForEntry(eidx, activities)
                if (firstActivity > 0) {
                    entryBoundActivity[eidx] = firstActivity
                }
            }
        }
    }

    /**
     * Find the bound (starting) activity for an entry.
     */
    private fun findBoundActivityForEntry(eidx: Int, activities: List<Int>): Int {
        // Check each activity to see if it's the entry point
        // The bound activity typically has PRE_SEQ or no predecessors
        for (aidx in activities) {
            val preType = lsn.actpretype.get(0, aidx).toInt()
            // If this is a starting activity (PRE_SEQ with no predecessors)
            // For now, use simple heuristic: first activity is the bound one
            if (preType == ActivityPrecedenceType.ID_PRE_SEQ || preType == 0) {
                // Check if any other activity points to this one
                var hasPredecessor = false
                for (otherAidx in activities) {
                    if (otherAidx != aidx && lsn.graph.get(otherAidx, aidx) > 0) {
                        hasPredecessor = true
                        break
                    }
                }
                if (!hasPredecessor) {
                    return aidx
                }
            }
        }
        // Fallback: return first activity
        return if (activities.isNotEmpty()) activities[0] else -1
    }

    /**
     * Find the parent entry for an activity.
     * This determines which entry the activity's utilization should be attributed to.
     */
    private fun findParentEntry(aidx: Int, tidx: Int, taskEntries: List<Int>): Int {
        // If task has only one entry, all activities belong to it
        if (taskEntries.size == 1) {
            return taskEntries[0]
        }

        // For multiple entries, find which entry's activity chain includes this activity
        // Check if this activity is a bound activity for any entry
        for (eidx in taskEntries) {
            val boundAct = entryBoundActivity[eidx]
            if (boundAct == aidx) {
                return eidx
            }
        }

        // Check if this activity is reachable from any entry's bound activity
        // by following the activity graph (within the same task)
        for (eidx in taskEntries) {
            val boundAct = entryBoundActivity[eidx] ?: continue
            if (isActivityReachableFrom(boundAct, aidx, tidx)) {
                return eidx
            }
        }

        // Fallback: return first entry
        return if (taskEntries.isNotEmpty()) taskEntries[0] else -1
    }

    /**
     * Check if targetAidx is reachable from startAidx in the activity graph.
     */
    private fun isActivityReachableFrom(startAidx: Int, targetAidx: Int, tidx: Int): Boolean {
        if (startAidx == targetAidx) return true

        val visited = mutableSetOf<Int>()
        val queue = LinkedList<Int>()
        queue.add(startAidx)

        while (queue.isNotEmpty()) {
            val current = queue.poll()
            if (current == targetAidx) return true
            if (visited.contains(current)) continue
            visited.add(current)

            val successors = activitySuccessors[current] ?: continue
            for (succ in successors) {
                // Only follow activities in the same task
                val succParent = lsn.parent.get(0, succ).toInt()
                if (succParent == tidx && !visited.contains(succ)) {
                    queue.add(succ)
                }
            }
        }
        return false
    }

    /**
     * Build the activity execution graph (successors for each activity).
     */
    private fun buildActivityGraph() {
        if (lsn.graph == null) return

        for (a in 1..nacts) {
            val aidx = ashift + a
            val successors = mutableListOf<Int>()

            // Check graph matrix for successors
            for (a2 in 1..nacts) {
                val aidx2 = ashift + a2
                if (lsn.graph.get(aidx, aidx2) > 0) {
                    successors.add(aidx2)
                }
            }

            if (successors.isNotEmpty()) {
                activitySuccessors[aidx] = successors
            }

            // Check for OR-fork probabilities
            val postType = lsn.actposttype.get(0, aidx).toInt()
            if (postType == ActivityPrecedenceType.ID_POST_OR && successors.size > 1) {
                // Extract probabilities from graph weights
                val probs = mutableListOf<Pair<Int, Double>>()
                var totalProb = 0.0
                for (succIdx in successors) {
                    val prob = lsn.graph.get(aidx, succIdx)
                    probs.add(Pair(succIdx, prob))
                    totalProb += prob
                }
                // Normalize if needed
                if (totalProb > 0 && totalProb != 1.0) {
                    for (i in probs.indices) {
                        probs[i] = Pair(probs[i].first, probs[i].second / totalProb)
                    }
                }
                orForkProbs[aidx] = probs
            }
        }
    }

    /**
     * Main simulation loop.
     */
    fun simulate(maxEventsParam: Double) {
        maxEvents = maxEventsParam.toLong()
        warmupEventThreshold = (maxEvents * 0.2).toLong()  // 20% warmup

        Sim.init()

        // Start REF tasks (workload generators)
        startRefTasks()

        // Run simulation
        Sim.start()
    }

    /**
     * Start reference (REF) tasks that generate workload.
     */
    private fun startRefTasks() {
        for (t in 0 until ntasks) {
            val tidx = tshift + t + 1
            val sched = lsn.sched[tidx]

            if (sched == SchedStrategy.REF) {
                // Get the first entry of this REF task
                val entries = lsn.entriesof[tidx]
                if (entries != null && entries.isNotEmpty()) {
                    val entryIdx = entries[0]

                    // Get think time distribution
                    @Suppress("DEPRECATION")
                    val thinkDist = lsn.think[tidx]

                    // Start initial requests based on task multiplicity
                    val mult = taskMultiplicity[t]
                    for (i in 0 until mult) {
                        val thinkTime = sampleDistribution(thinkDist)
                        RefTaskThinkComplete(tidx, entryIdx).schedule(thinkTime)
                    }
                }
            }
        }
    }

    /**
     * Sample from a distribution, returning the sampled value.
     */
    private fun sampleDistribution(dist: Distribution?): Double {
        if (dist == null) return 0.0
        if (dist.isImmediate) return 0.0
        val samples = dist.sample(1, random)
        return if (samples.isNotEmpty()) samples[0].coerceAtLeast(0.0) else 0.0
    }

    /**
     * Update host statistics before state change.
     */
    private fun updateHostStats(hostIdx: Int) {
        val h = hostIdx - 1  // Convert to 0-based
        if (h < 0 || h >= nhosts) return

        val now = Sim.time()
        val elapsed = now - lastHostUpdateTime[h]
        if (elapsed > 0 && warmupDone) {
            // For utilization, track busy server count.
            // For INF scheduling, don't cap - utilization can exceed 1.0
            // For finite servers (PS, FCFS, etc.), cap at the number of servers
            val strategy = hostSchedStrategy[h]
            val busyServers = if (isINFScheduling(strategy)) {
                hostBusy[h]
            } else {
                minOf(hostBusy[h], hostServers[h])
            }
            hostBusyTime[h] += busyServers.toDouble() * elapsed
        }
        lastHostUpdateTime[h] = now
    }

    /**
     * Update task statistics before state change.
     */
    private fun updateTaskStats(taskIdx: Int) {
        val t = taskIdx - tshift - 1  // Convert to 0-based
        if (t < 0 || t >= ntasks) return

        val now = Sim.time()
        val elapsed = now - lastTaskUpdateTime[t]
        if (elapsed > 0 && warmupDone) {
            taskUtilTime[t] += taskBusy[t].toDouble() * elapsed
        }
        lastTaskUpdateTime[t] = now
    }

    /**
     * Update activity statistics before state change.
     */
    private fun updateActivityStats(activityIdx: Int) {
        val a = activityIdx - ashift - 1  // Convert to 0-based
        if (a < 0 || a >= nacts) return

        val now = Sim.time()
        val elapsed = now - lastActivityUpdateTime[a]
        if (elapsed > 0 && warmupDone) {
            activityBusyTime[a] += activityBusy[a].toDouble() * elapsed
        }
        lastActivityUpdateTime[a] = now
    }

    /**
     * Update task queue length statistics before state change.
     * Call this before incrementing/decrementing taskQueueLength[t].
     */
    private fun updateTaskQueueLengthStats(t: Int) {
        if (t < 0 || t >= ntasks) return

        val now = Sim.time()
        val elapsed = now - lastTaskQueueUpdateTime[t]
        if (elapsed > 0 && warmupDone) {
            taskQueueLengthTime[t] += taskQueueLength[t].toDouble() * elapsed
        }
        lastTaskQueueUpdateTime[t] = now
    }

    /**
     * Update host queue length statistics before state change.
     * Call this before incrementing/decrementing hostQueueLength[h].
     */
    private fun updateHostQueueLengthStats(h: Int) {
        if (h < 0 || h >= nhosts) return

        val now = Sim.time()
        val elapsed = now - lastHostQueueUpdateTime[h]
        if (elapsed > 0 && warmupDone) {
            hostQueueLengthTime[h] += hostQueueLength[h].toDouble() * elapsed
        }
        lastHostQueueUpdateTime[h] = now
    }

    /**
     * Check if simulation should stop based on event count.
     */
    private fun checkEventCountStop() {
        if (!warmupDone && totalEventCount >= warmupEventThreshold) {
            warmupDone = true
            warmupEndTime = Sim.time()
            resetStatistics()
        }

        if (totalEventCount >= maxEvents) {
            Sim.stop()
        }
    }

    /**
     * Reset statistics after warmup period.
     */
    private fun resetStatistics() {
        for (h in 0 until nhosts) {
            hostBusyTime[h] = 0.0
            hostCompletions[h] = 0
            lastHostUpdateTime[h] = Sim.time()
            // Reset host queue length time accumulator (but keep current queue length)
            hostQueueLengthTime[h] = 0.0
            lastHostQueueUpdateTime[h] = Sim.time()
        }
        for (t in 0 until ntasks) {
            taskUtilTime[t] = 0.0
            taskCompletions[t] = 0.0
            lastTaskUpdateTime[t] = Sim.time()
            // Reset task queue length time accumulator (but keep current queue length)
            taskQueueLengthTime[t] = 0.0
            lastTaskQueueUpdateTime[t] = Sim.time()
        }
        for (e in 0 until nentries) {
            entryThroughput[e] = 0.0
            entryResponseTimeTally[e].init()
        }
        for (a in 0 until nacts) {
            activityThroughput[a] = 0.0
            activityServiceTimeTally[a].init()
            activityBusyTime[a] = 0.0
            lastActivityUpdateTime[a] = Sim.time()
        }
        // Clear activity start times from warmup period
        activityStartTimes.clear()
    }

    /**
     * Create a new request and register it.
     */
    private fun createRequest(
        entryIdx: Int,
        callType: CallType,
        callerRequestId: Long,
        callerTaskIdx: Int,
        callMean: Double = 1.0
    ): LNRequest {
        val requestId = nextRequestId++
        val request = LNRequest(
            requestId = requestId,
            entryIdx = entryIdx,
            callType = callType,
            arrivalTime = Sim.time(),
            taskArrivalTime = Sim.time(),
            callerRequestId = callerRequestId,
            callerTaskIdx = callerTaskIdx,
            currentActivityIdx = -1,
            parentForkId = -1
        )
        activeRequests[requestId] = request
        return request
    }

    /**
     * Request arrives at a task's entry.
     */
    private fun arriveAtTask(tidx: Int, request: LNRequest) {
        val t = tidx - tshift - 1
        if (t < 0 || t >= ntasks) return

        updateTaskStats(tidx)
        // Update queue length stats before incrementing
        updateTaskQueueLengthStats(t)
        taskQueueLength[t]++  // Job arrives at task
        request.taskArrivalTime = Sim.time()

        val strategy = taskSchedStrategy[t]
        val mult = taskMultiplicity[t]

        // INF scheduling - unlimited task instances
        if (isINFScheduling(strategy)) {
            taskBusy[t]++
            startRequestProcessing(tidx, request)
            return
        }

        // Standard scheduling (FCFS, LCFS, PS)
        // Note: For tasks, PS behavior is approximated by standard queueing since
        // task service is defined by activity execution, not a fixed service time
        if (taskBusy[t] < mult) {
            // Task instance available - start processing
            taskBusy[t]++
            startRequestProcessing(tidx, request)
        } else {
            // Queue the request using priority queue with appropriate comparator
            taskQueues[t].add(request)
        }
    }

    /**
     * Start processing a request at an entry.
     */
    private fun startRequestProcessing(tidx: Int, request: LNRequest) {
        val entryIdx = request.entryIdx

        // Find the bound activity for this entry
        val boundActivityIdx = entryBoundActivity[entryIdx]
        if (boundActivityIdx != null && boundActivityIdx > 0) {
            request.currentActivityIdx = boundActivityIdx
            startActivityExecution(request, boundActivityIdx)
        } else {
            // No bound activity found - complete the request immediately
            completeRequest(request)
        }
    }

    /**
     * Start executing an activity.
     */
    private fun startActivityExecution(request: LNRequest, aidx: Int) {
        request.currentActivityIdx = aidx

        // Track activity start for utilization and response time
        val a = aidx - ashift - 1
        if (a >= 0 && a < nacts) {
            updateActivityStats(aidx)
            activityBusy[a]++
            activityStartTimes[Pair(request.requestId, aidx)] = Sim.time()
        }

        // Get the task and host for this activity
        val tidx = lsn.parent.get(0, aidx).toInt()
        val hidx = lsn.parent.get(0, tidx).toInt()

        // Get host demand distribution for this activity
        @Suppress("DEPRECATION")
        val hostDemandDist = lsn.hostdem[aidx]
        val hostDemand = sampleDistribution(hostDemandDist)

        if (hostDemand > 0) {
            // Queue for host processing
            queueForHost(hidx, request, aidx, hostDemand)
        } else {
            // No host demand - proceed to calls immediately
            executeActivityCalls(request, aidx)
        }
    }

    /**
     * Queue request for host (processor) service.
     */
    private fun queueForHost(hidx: Int, request: LNRequest, aidx: Int, serviceTime: Double) {
        val h = hidx - 1
        if (h < 0 || h >= nhosts) {
            // Invalid host - proceed without host processing
            executeActivityCalls(request, aidx)
            return
        }

        updateHostStats(hidx)

        val hostRequest = HostRequest(request, aidx, serviceTime, Sim.time())
        val strategy = hostSchedStrategy[h]

        when {
            isINFScheduling(strategy) -> {
                // INF: Infinite server - start immediately, no queuing
                hostBusy[h]++
                ActivityHostComplete(hostRequest, hidx).schedule(serviceTime)
            }
            isPSScheduling(strategy) -> {
                // PS: Processor Sharing - add to PS jobs and update sharing
                startPSJob(h, hidx, hostRequest)
            }
            else -> {
                // FCFS, LCFS, etc. - standard queuing behavior
                if (hostBusy[h] < hostServers[h]) {
                    // Server available
                    hostBusy[h]++
                    ActivityHostComplete(hostRequest, hidx).schedule(serviceTime)
                } else {
                    // Add to host queue
                    hostQueues[h].add(hostRequest)
                }
            }
        }
    }

    /**
     * Start a PS (Processor Sharing) job at a host.
     */
    private fun startPSJob(h: Int, hidx: Int, hostRequest: HostRequest) {
        val psJob = PSJob(
            hostRequest = hostRequest,
            remainingService = hostRequest.serviceTime,
            lastUpdateTime = Sim.time()
        )
        hostPSJobs[h].add(psJob)
        hostBusy[h]++

        // Schedule or reschedule PS update event
        schedulePSUpdate(h, hidx)
    }

    /**
     * Schedule or reschedule PS update event for a host.
     */
    private fun schedulePSUpdate(h: Int, hidx: Int) {
        // Cancel existing event if any
        hostPSEvent[h]?.cancel()

        val jobs = hostPSJobs[h]
        if (jobs.isEmpty()) {
            hostPSEvent[h] = null
            return
        }

        // Calculate time until next completion
        val numJobs = jobs.size
        val numServers = hostServers[h]
        val sharePerJob = minOf(1.0, numServers.toDouble() / numJobs)

        // Find minimum remaining time considering current shares
        var minCompletionTime = Double.MAX_VALUE
        val now = Sim.time()

        for (job in jobs) {
            // Update remaining service based on elapsed time
            val elapsed = now - job.lastUpdateTime
            val served = elapsed * sharePerJob
            job.remainingService = (job.remainingService - served).coerceAtLeast(0.0)
            job.lastUpdateTime = now

            // Time to complete at current rate
            val timeToComplete = if (sharePerJob > 0) job.remainingService / sharePerJob else Double.MAX_VALUE
            if (timeToComplete < minCompletionTime) {
                minCompletionTime = timeToComplete
            }
        }

        // Schedule next PS update event
        if (minCompletionTime < Double.MAX_VALUE && minCompletionTime > 0) {
            val event = PSUpdateEvent(h, hidx)
            event.schedule(minCompletionTime)
            hostPSEvent[h] = event
        } else if (minCompletionTime == 0.0) {
            // Process immediately
            processPSCompletions(h, hidx)
        }
    }

    /**
     * Process PS completions at a host.
     */
    private fun processPSCompletions(h: Int, hidx: Int) {
        val jobs = hostPSJobs[h]
        val numServers = hostServers[h]
        val numJobs = jobs.size
        val sharePerJob = if (numJobs > 0) minOf(1.0, numServers.toDouble() / numJobs) else 1.0
        val now = Sim.time()

        // Update all jobs and find completed ones
        val completedJobs = mutableListOf<PSJob>()
        for (job in jobs) {
            val elapsed = now - job.lastUpdateTime
            val served = elapsed * sharePerJob
            job.remainingService = (job.remainingService - served).coerceAtLeast(0.0)
            job.lastUpdateTime = now

            if (job.remainingService <= 1e-9) {
                completedJobs.add(job)
            }
        }

        // Process completed jobs
        for (job in completedJobs) {
            jobs.remove(job)
            updateHostStats(hidx)
            hostBusy[h]--
            if (warmupDone) {
                hostCompletions[h]++
            }
            // Execute calls for the completed activity
            executeActivityCalls(job.hostRequest.request, job.hostRequest.activityIdx)
        }

        // Reschedule PS update if there are remaining jobs
        if (jobs.isNotEmpty()) {
            schedulePSUpdate(h, hidx)
        } else {
            hostPSEvent[h] = null
        }
    }

    /**
     * Activity host demand completes - proceed to calls.
     */
    private fun completeActivityHostDemand(hostRequest: HostRequest, hidx: Int) {
        val h = hidx - 1
        if (h < 0 || h >= nhosts) return

        updateHostStats(hidx)
        hostBusy[h]--
        if (warmupDone) {
            hostCompletions[h]++
        }

        // Try to start next queued request at this host
        if (hostQueues[h].isNotEmpty()) {
            val nextHostRequest = hostQueues[h].poll()
            hostBusy[h]++
            ActivityHostComplete(nextHostRequest, hidx).schedule(nextHostRequest.serviceTime)
        }

        // Execute calls for this activity
        executeActivityCalls(hostRequest.request, hostRequest.activityIdx)
    }

    /**
     * Execute calls made by an activity.
     * Sync calls are executed sequentially (one at a time, waiting for reply before next).
     * Async calls are executed in parallel (fire-and-forget).
     */
    private fun executeActivityCalls(request: LNRequest, aidx: Int) {
        val calls = lsn.callsof[aidx]

        if (calls == null || calls.isEmpty()) {
            // No calls - complete the activity
            completeActivity(request, aidx)
            return
        }

        // Build queue of sync calls and execute async calls immediately
        val syncCallQueue = LinkedList<Int>()

        // Process each call
        for (cidx in calls) {
            val callType = lsn.calltype[cidx] ?: continue
            // callpair columns: 1 = calling activity, 2 = target entry
            val targetEntryIdx = lsn.callpair.get(cidx, 2).toInt()  // Column 2 has target entry
            if (targetEntryIdx <= 0) continue

            // Get call mean (number of calls)
            val meanCalls = lsn.callproc_mean[cidx] ?: 1.0
            val numCalls = sampleCallCount(meanCalls)

            for (c in 0 until numCalls) {
                when (callType) {
                    CallType.SYNC -> {
                        // Queue sync calls for sequential execution
                        syncCallQueue.add(targetEntryIdx)
                    }
                    CallType.ASYNC -> {
                        // Execute async calls immediately (fire-and-forget)
                        executeAsyncCall(request, targetEntryIdx, meanCalls)
                    }
                    else -> { /* FWD not implemented */ }
                }
            }
        }

        if (syncCallQueue.isNotEmpty()) {
            // Store the queue and activity for when replies come back
            pendingSyncCallQueue[request.requestId] = syncCallQueue
            syncCallActivity[request.requestId] = aidx
            syncCallInFlight[request.requestId] = false

            // Execute the first sync call
            executeNextSyncCall(request, aidx)
        } else {
            // No sync calls - complete the activity
            completeActivity(request, aidx)
        }
    }

    /**
     * Execute the next sync call in the queue for a request.
     */
    private fun executeNextSyncCall(request: LNRequest, aidx: Int) {
        val queue = pendingSyncCallQueue[request.requestId] ?: return

        if (queue.isEmpty()) {
            // All sync calls completed - complete the activity
            pendingSyncCallQueue.remove(request.requestId)
            syncCallActivity.remove(request.requestId)
            syncCallInFlight.remove(request.requestId)
            completeActivity(request, aidx)
            return
        }

        // Get the next target entry
        val targetEntryIdx = queue.poll()
        syncCallInFlight[request.requestId] = true

        // Execute the sync call
        executeSyncCall(request, aidx, targetEntryIdx, 1.0)
    }

    /**
     * Sample the number of calls to make (for stochastic call counts).
     */
    private fun sampleCallCount(meanCalls: Double): Int {
        // For now, use deterministic rounding
        // Could extend to Poisson or geometric sampling
        return meanCalls.toInt().coerceAtLeast(1)
    }

    /**
     * Execute a synchronous call - caller blocks until reply.
     */
    private fun executeSyncCall(callerRequest: LNRequest, callerActivityIdx: Int, targetEntryIdx: Int, callMean: Double = 1.0) {
        val targetTidx = lsn.parent.get(0, targetEntryIdx).toInt()

        val callRequest = createRequest(
            entryIdx = targetEntryIdx,
            callType = CallType.SYNC,
            callerRequestId = callerRequest.requestId,
            callerTaskIdx = lsn.parent.get(0, callerRequest.entryIdx).toInt(),
            callMean = callMean
        )

        // Arrive at target task
        arriveAtTask(targetTidx, callRequest)
    }

    /**
     * Execute an asynchronous call - fire and forget.
     */
    private fun executeAsyncCall(callerRequest: LNRequest, targetEntryIdx: Int, callMean: Double = 1.0) {
        val targetTidx = lsn.parent.get(0, targetEntryIdx).toInt()

        val callRequest = createRequest(
            entryIdx = targetEntryIdx,
            callType = CallType.ASYNC,
            callerRequestId = -1,  // No reply expected
            callerTaskIdx = -1,
            callMean = callMean
        )

        // Arrive at target task
        arriveAtTask(targetTidx, callRequest)
    }

    /**
     * Handle synchronous call reply.
     * Executes the next sync call in the queue, or completes the activity if all done.
     */
    private fun handleSyncReply(completedRequest: LNRequest) {
        val callerRequestId = completedRequest.callerRequestId
        if (callerRequestId < 0) return

        val callerRequest = activeRequests[callerRequestId] ?: return

        // Mark that the sync call has completed
        syncCallInFlight[callerRequestId] = false

        // Get the activity that made the call
        val callerActivityIdx = syncCallActivity[callerRequestId] ?: return

        // Execute the next sync call or complete the activity
        executeNextSyncCall(callerRequest, callerActivityIdx)
    }

    /**
     * Complete an activity and handle precedence.
     */
    private fun completeActivity(request: LNRequest, aidx: Int) {
        val a = aidx - ashift - 1
        if (a >= 0 && a < nacts) {
            // Update activity utilization stats
            updateActivityStats(aidx)
            activityBusy[a]--

            if (warmupDone) {
                // Count each activity completion
                activityThroughput[a] += 1.0

                // Record activity response time
                val startTimeKey = Pair(request.requestId, aidx)
                val startTime = activityStartTimes.remove(startTimeKey)
                if (startTime != null) {
                    val responseTime = Sim.time() - startTime
                    activityServiceTimeTally[a].add(responseTime)
                }
            } else {
                // Remove start time record if still in warmup
                activityStartTimes.remove(Pair(request.requestId, aidx))
            }
        }

        totalEventCount++
        checkEventCountStop()

        // Get post-precedence type
        val postType = lsn.actposttype.get(0, aidx).toInt()
        val successors = activitySuccessors[aidx]

        when (postType) {
            ActivityPrecedenceType.ID_POST_SEQ -> {
                // Sequential: find successor and execute
                if (successors != null && successors.isNotEmpty()) {
                    val successor = successors[0]
                    moveToSuccessorOrJoin(request, successor)
                } else {
                    // No successor - check if this is the reply activity
                    checkAndCompleteRequest(request, aidx)
                }
            }
            ActivityPrecedenceType.ID_POST_AND -> {
                // AND-fork: start all parallel branches
                if (successors != null && successors.size > 1) {
                    handleAndFork(request, aidx, successors)
                } else if (successors != null && successors.size == 1) {
                    // Single successor - use moveToSuccessorOrJoin to detect join points
                    // This handles the case where actposttype is incorrectly set on fork branches
                    moveToSuccessorOrJoin(request, successors[0])
                } else {
                    checkAndCompleteRequest(request, aidx)
                }
            }
            ActivityPrecedenceType.ID_POST_OR -> {
                // OR-fork: probabilistically choose one branch
                if (successors != null && successors.isNotEmpty()) {
                    val chosenSuccessor = chooseOrForkBranch(aidx, successors)
                    startActivityExecution(request, chosenSuccessor)
                } else {
                    checkAndCompleteRequest(request, aidx)
                }
            }
            else -> {
                // No precedence or unknown - detect fork from successor count and probabilities
                if (successors != null && successors.size > 1) {
                    // Check if this is an OR-fork (has probabilities registered)
                    val probs = orForkProbs[aidx]
                    if (probs != null && probs.isNotEmpty()) {
                        // OR-fork detected from probabilities - choose one branch probabilistically
                        val chosenSuccessor = chooseOrForkBranch(aidx, successors)
                        startActivityExecution(request, chosenSuccessor)
                    } else {
                        // Check if edge weights indicate an OR-fork (probabilities not all 1.0)
                        val weights = successors.map { lsn.graph.get(aidx, it) }
                        val hasNonUniformWeights = weights.any { it != 1.0 && it > 0.0 }

                        if (hasNonUniformWeights) {
                            // OR-fork with probabilities in graph - populate orForkProbs and choose
                            val newProbs = mutableListOf<Pair<Int, Double>>()
                            var totalProb = 0.0
                            for (succIdx in successors) {
                                val prob = lsn.graph.get(aidx, succIdx)
                                newProbs.add(Pair(succIdx, prob))
                                totalProb += prob
                            }
                            // Normalize if needed
                            if (totalProb > 0 && totalProb != 1.0) {
                                for (i in newProbs.indices) {
                                    newProbs[i] = Pair(newProbs[i].first, newProbs[i].second / totalProb)
                                }
                            }
                            orForkProbs[aidx] = newProbs
                            val chosenSuccessor = chooseOrForkBranch(aidx, successors)
                            startActivityExecution(request, chosenSuccessor)
                        } else {
                            // Multiple successors with uniform weights - could be AND-fork
                            // For now, just start first successor
                            moveToSuccessorOrJoin(request, successors[0])
                        }
                    }
                } else if (successors != null && successors.isNotEmpty()) {
                    moveToSuccessorOrJoin(request, successors[0])
                } else {
                    checkAndCompleteRequest(request, aidx)
                }
            }
        }
    }

    /**
     * Move to successor activity, handling AND-join detection.
     * If this request is part of an AND-fork branch and the successor is a join point
     * (has multiple incoming edges or has PRE_AND type), signal branch completion.
     */
    private fun moveToSuccessorOrJoin(request: LNRequest, successorIdx: Int) {
        // Check if this is a branch request reaching a join point
        if (request.parentForkId >= 0) {
            val successorPreType = lsn.actpretype.get(0, successorIdx).toInt()

            // Join point detection: marked as PRE_AND
            if (successorPreType == ActivityPrecedenceType.ID_PRE_AND) {
                // This is a branch request reaching a join point - signal branch completion
                handleAndJoinBranchComplete(request)
                return
            }
        }

        // Normal successor transition
        startActivityExecution(request, successorIdx)
    }

    /**
     * Handle AND-fork: start parallel branches.
     */
    // Debug flag - set to true to trace fork/join execution
    private val DEBUG_FORK_JOIN = false

    private fun handleAndFork(request: LNRequest, forkActivityIdx: Int, successors: List<Int>) {
        val parentForkId = nextRequestId++

        val joinActivityIdx = findAndJoinActivity(forkActivityIdx, successors)

        if (DEBUG_FORK_JOIN) {
            val forkName = lsn.hashnames.get(forkActivityIdx)
            val joinName = if (joinActivityIdx > 0) lsn.hashnames.get(joinActivityIdx) else "none"
            println("[FORK] time=${Sim.time()} $forkName -> ${successors.map { lsn.hashnames.get(it) }}, join=$joinName")
        }

        // Create AND-join state
        pendingAndJoins[parentForkId] = AndJoinState(
            parentRequest = request,
            joinActivityIdx = joinActivityIdx,
            pendingBranches = successors.size,
            branchCompletionTimes = mutableListOf()
        )

        // Start each branch
        for (succIdx in successors) {
            val branchRequest = request.copy(
                requestId = nextRequestId++,
                currentActivityIdx = succIdx,
                parentForkId = parentForkId
            )
            activeRequests[branchRequest.requestId] = branchRequest
            startActivityExecution(branchRequest, succIdx)
        }
    }

    /**
     * Find the AND-join activity for an AND-fork.
     */
    private fun findAndJoinActivity(forkActivityIdx: Int, successors: List<Int>): Int {
        // The AND-join activity is the one where all branches converge
        // It has PRE_AND precedence type
        for (a in 1..nacts) {
            val aidx = ashift + a
            val preType = lsn.actpretype.get(0, aidx).toInt()
            if (preType == ActivityPrecedenceType.ID_PRE_AND) {
                // Check if all fork successors eventually lead here
                // For simplicity, return the first PRE_AND activity found
                return aidx
            }
        }
        return -1
    }

    /**
     * Handle AND-join: check if all branches are complete.
     */
    private fun handleAndJoinBranchComplete(request: LNRequest) {
        val parentForkId = request.parentForkId
        if (parentForkId < 0) return

        val joinState = pendingAndJoins[parentForkId] ?: return

        val prevPending = joinState.pendingBranches
        joinState.pendingBranches--
        joinState.branchCompletionTimes.add(Sim.time())

        if (DEBUG_FORK_JOIN) {
            val actName = lsn.hashnames.get(request.currentActivityIdx)
            println("[JOIN-BRANCH] time=${Sim.time()} $actName completed, pending=${prevPending}->${joinState.pendingBranches}")
        }

        // Clean up branch request
        activeRequests.remove(request.requestId)

        if (joinState.pendingBranches <= 0) {
            // All branches complete - continue with join activity
            pendingAndJoins.remove(parentForkId)

            if (DEBUG_FORK_JOIN) {
                val joinName = if (joinState.joinActivityIdx > 0) lsn.hashnames.get(joinState.joinActivityIdx) else "none"
                println("[JOIN-COMPLETE] time=${Sim.time()} all branches done, starting join $joinName")
            }

            val parentRequest = joinState.parentRequest
            if (joinState.joinActivityIdx > 0) {
                startActivityExecution(parentRequest, joinState.joinActivityIdx)
            } else {
                // No join activity - complete the request
                checkAndCompleteRequest(parentRequest, parentRequest.currentActivityIdx)
            }
        }
    }

    /**
     * Choose an OR-fork branch probabilistically.
     */
    private fun chooseOrForkBranch(aidx: Int, successors: List<Int>): Int {
        val probs = orForkProbs[aidx]
        if (probs != null && probs.isNotEmpty()) {
            val u = random.nextDouble()
            var cumProb = 0.0
            for ((succIdx, prob) in probs) {
                cumProb += prob
                if (u <= cumProb) {
                    return succIdx
                }
            }
            return probs.last().first
        }
        // Uniform random if no probabilities defined
        return successors[random.nextInt(successors.size)]
    }

    /**
     * Check if request should complete (reached reply activity or no successors).
     */
    private fun checkAndCompleteRequest(request: LNRequest, aidx: Int) {
        // Check if this is an AND-fork branch
        if (request.parentForkId >= 0) {
            handleAndJoinBranchComplete(request)
            return
        }

        // Check if this is the reply activity for the entry
        val entryIdx = request.entryIdx
        val replyActivity = entryReplyActivity[entryIdx]

        if (replyActivity == null || replyActivity == aidx) {
            // This is the reply activity or no reply activity defined - complete
            completeRequest(request)
        } else {
            // Not at reply activity yet - this might be an error in the model
            // For robustness, complete anyway
            completeRequest(request)
        }
    }

    /**
     * Complete a request and handle reply for sync calls.
     */
    private fun completeRequest(request: LNRequest) {
        val entryIdx = request.entryIdx
        val e = entryIdx - eshift - 1

        val responseTime = Sim.time() - request.taskArrivalTime

        // Record entry statistics
        if (e >= 0 && e < nentries && warmupDone) {
            entryThroughput[e] += 1.0
            entryResponseTimeTally[e].add(responseTime)
        }

        // Release task instance
        val tidx = lsn.parent.get(0, entryIdx).toInt()
        val t = tidx - tshift - 1
        if (t >= 0 && t < ntasks) {
            updateTaskStats(tidx)
            // Update queue length stats before decrementing
            updateTaskQueueLengthStats(t)
            taskQueueLength[t]--  // Job leaves task
            taskBusy[t]--
            if (warmupDone) {
                taskCompletions[t] += 1.0
            }

            // Try to start next queued request at this task
            if (taskQueues[t].isNotEmpty()) {
                val nextRequest = taskQueues[t].poll()
                taskBusy[t]++
                startRequestProcessing(tidx, nextRequest)
            }
        }

        // Handle reply for sync call
        if (request.callType == CallType.SYNC && request.callerRequestId >= 0) {
            handleSyncReply(request)
        }

        // If this request originated from a REF task, schedule next think time
        if (request.callerRequestId < 0 && request.callerTaskIdx < 0) {
            // This is a REF task request - schedule next cycle
            val sched = lsn.sched[tidx]
            if (sched == SchedStrategy.REF) {
                scheduleNextRefTaskRequest(tidx, entryIdx)
            }
        }

        // Clean up
        activeRequests.remove(request.requestId)
    }

    /**
     * Extract results into LNDESResult.
     */
    fun getLNDESResult(): LNDESResult {
        val result = LNDESResult()
        result.lsn = lsn

        val simTime = Sim.time() - warmupEndTime
        if (simTime <= 0) {
            return result
        }

        // Initialize result matrices
        result.QLN = Matrix(1, lsn.nidx + 1)
        result.ULN = Matrix(1, lsn.nidx + 1)
        result.RLN = Matrix(1, lsn.nidx + 1)
        result.WLN = Matrix(1, lsn.nidx + 1)
        result.TLN = Matrix(1, lsn.nidx + 1)
        result.ALN = Matrix(1, lsn.nidx + 1)

        // Host metrics (indices 1 to nhosts)
        // LQNS reports aggregate processor utilization (total busy time / sim time),
        // which can exceed 1.0 for multi-server processors.
        for (h in 0 until nhosts) {
            val hidx = h + 1
            val util = hostBusyTime[h] / simTime
            val tput = hostCompletions[h].toDouble() / simTime
            result.ULN.set(0, hidx, util)
            result.TLN.set(0, hidx, tput)
        }

        // Task metrics (indices tshift+1 to tshift+ntasks)
        // In LQN convention, task utilization = sum of activity utilizations on that task
        // This is computed after activity metrics are calculated (below)
        for (t in 0 until ntasks) {
            val tidx = tshift + t + 1
            // Initialize task utilization to 0 (will be summed from activities)
            result.ULN.set(0, tidx, 0.0)
            val tput = taskCompletions[t] / simTime
            result.TLN.set(0, tidx, tput)
            // Finalize time-weighted queue length and compute average
            updateTaskQueueLengthStats(t)
            val avgQLen = taskQueueLengthTime[t] / simTime
            result.QLN.set(0, tidx, avgQLen)
        }

        // Entry metrics (indices eshift+1 to eshift+nentries)
        // Initialize entry utilization to 0 (will be summed from activities below)
        for (e in 0 until nentries) {
            val eidx = eshift + e + 1
            val tput = entryThroughput[e] / simTime
            val respT = if (entryResponseTimeTally[e].numberObs() > 0) {
                entryResponseTimeTally[e].average()
            } else {
                0.0
            }
            result.TLN.set(0, eidx, tput)
            result.RLN.set(0, eidx, respT)
            result.ALN.set(0, eidx, tput)  // Arrival rate equals throughput in steady state
            // Queue length for entry using Little's Law: Q = λ * R
            val qLen = tput * respT
            result.QLN.set(0, eidx, qLen)
            // Initialize entry utilization to 0
            result.ULN.set(0, eidx, 0.0)
        }

        // Activity metrics (indices ashift+1 to ashift+nacts)
        // Two-pass algorithm:
        // Pass 1: Compute processor utilization contributions and accumulate to task/entry
        // Pass 2: Set activity utilization normalized to task utilization (LQNS convention)

        // Store processor utilization contributions for pass 2
        val activityProcUtil = DoubleArray(nacts)

        // Pass 1: Compute throughput, response time, queue length, and accumulate processor utilization
        for (a in 0 until nacts) {
            val aidx = ashift + a + 1
            val tput = activityThroughput[a] / simTime
            result.TLN.set(0, aidx, tput)

            // Processor utilization contribution: throughput × mean host demand
            val meanHostDem = lsn.hostdem_mean[aidx] ?: 0.0
            val procUtil = tput * meanHostDem
            activityProcUtil[a] = procUtil

            // Add activity's processor utilization to parent task utilization
            val tidx = lsn.parent.get(0, aidx).toInt()
            if (tidx > 0 && tidx <= lsn.nidx) {
                val currentTaskUtil = result.ULN.get(0, tidx)
                result.ULN.set(0, tidx, currentTaskUtil + procUtil)

                // Add activity utilization to parent entry utilization
                // For tasks with a single entry, all activities belong to that entry
                val taskEntries = lsn.entriesof[tidx]
                if (taskEntries != null && taskEntries.isNotEmpty()) {
                    // Distribute activity utilization among entries
                    // For simplicity, if task has one entry, add to that entry
                    // For multiple entries, add to the entry whose bound activity starts this activity's chain
                    val parentEntryIdx = findParentEntry(aidx, tidx, taskEntries)
                    if (parentEntryIdx > 0) {
                        val currentEntryUtil = result.ULN.get(0, parentEntryIdx)
                        result.ULN.set(0, parentEntryIdx, currentEntryUtil + procUtil)
                    }
                }
            }

            // Activity response time: average from tally
            val respT = if (activityServiceTimeTally[a].numberObs() > 0) {
                activityServiceTimeTally[a].average()
            } else {
                0.0
            }
            result.RLN.set(0, aidx, respT)
            // Queue length for activity using Little's Law: Q = λ * R
            val qLen = tput * respT
            result.QLN.set(0, aidx, qLen)
        }

        // Pass 2: Set activity utilization normalized to task utilization (LQNS convention)
        // Activity utilization = activity's execution time / task's busy time
        // When activity is always busy during task execution, this equals 1.0
        for (a in 0 until nacts) {
            val aidx = ashift + a + 1
            val procUtil = activityProcUtil[a]
            val tidx = lsn.parent.get(0, aidx).toInt()
            val taskUtil = if (tidx > 0 && tidx <= lsn.nidx) result.ULN.get(0, tidx) else 0.0

            val activityUtil = if (taskUtil > 0.0 && taskUtil < 1.0) {
                // Normalize activity's processor utilization by task utilization
                Math.min(1.0, procUtil / taskUtil)
            } else if (procUtil > 0.0) {
                // Task is fully utilized or activity has demand: activity is always busy
                1.0
            } else {
                0.0
            }
            result.ULN.set(0, aidx, activityUtil)
        }

        return result
    }

    // ==================== Event Classes ====================

    /**
     * REF task think time completion - generates new request.
     */
    private inner class RefTaskThinkComplete(
        private val tidx: Int,
        private val entryIdx: Int
    ) : Event() {
        override fun actions() {
            // Create new request to entry
            val request = createRequest(
                entryIdx = entryIdx,
                callType = CallType.ASYNC,  // REF task requests are not blocking calls
                callerRequestId = -1,
                callerTaskIdx = -1
            )

            // Arrive at the task
            arriveAtTask(tidx, request)
        }
    }

    /**
     * REF task request completion - schedule next think time.
     * This is called when a request from a REF task completes.
     */
    private fun scheduleNextRefTaskRequest(tidx: Int, entryIdx: Int) {
        @Suppress("DEPRECATION")
        val thinkDist = lsn.think[tidx]
        val thinkTime = sampleDistribution(thinkDist)
        RefTaskThinkComplete(tidx, entryIdx).schedule(thinkTime)
    }

    /**
     * Activity host demand completes.
     */
    private inner class ActivityHostComplete(
        private val hostRequest: HostRequest,
        private val hidx: Int
    ) : Event() {
        override fun actions() {
            completeActivityHostDemand(hostRequest, hidx)
        }
    }

    /**
     * PS (Processor Sharing) update event - checks for completions and reschedules.
     */
    private inner class PSUpdateEvent(
        private val h: Int,
        private val hidx: Int
    ) : Event() {
        override fun actions() {
            processPSCompletions(h, hidx)
        }
    }
}
