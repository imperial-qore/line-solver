/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.des.analyzers

import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.ProcessType
import jline.lang.constant.SchedStrategy
import jline.lang.nodes.Queue
import jline.solvers.SolverOptions
import jline.solvers.des.DESResult
import jline.solvers.des.SolverDES
import jline.solvers.des.handlers.solver_ssj
import jline.solvers.des.handlers.solver_ssj_transient
import jline.streaming.Collector

/**
 * @brief DES analyzer that selects and executes the appropriate analysis method.
 *
 * @details This function serves as the main entry point for DES analysis, validating
 * the network structure and dispatching to the appropriate backend solver (SSJ).
 *
 * @section analyzer_validation Network Validation
 * The analyzer validates that the network is suitable for DES simulation:
 * - Open networks: Must have Source and Sink nodes
 * - Closed networks: Must have service nodes with finite job populations
 * - Mixed networks: Combination of open and closed classes
 * - Petri nets: Must have Place and/or Transition nodes
 *
 * @section analyzer_scheduling Supported Scheduling
 * Queue nodes must use one of the supported scheduling strategies:
 * - FCFS, HOL, FCFSPRIO (FIFO variants)
 * - LCFS, LCFSPR, LCFSPI (LIFO variants with priority)
 * - PS, DPS, GPS (processor sharing variants)
 * - SIRO, SJF, LJF, SEPT, LEPT (other disciplines)
 * - EXT (external arrivals)
 *
 * @section analyzer_distributions Supported Distributions
 * Service and arrival processes must be one of:
 * - EXP, ERLANG, HYPEREXP (exponential family)
 * - PH, APH, COXIAN (phase-type)
 * - IMMEDIATE, DISABLED (special cases)
 * - REPLAYER (trace-driven)
 *
 * @section analyzer_modes Analysis Modes
 * - **Steady-state**: When options.timespan is null or infinite
 * - **Transient**: When options.timespan specifies a finite time horizon
 *
 * @param sn Network structure representing the queueing network model
 * @param options Solver configuration options including method selection
 * @param solverDES DES solver instance to use for analysis
 * @return DES analysis results including method used and runtime
 * @throws RuntimeException If the network is not supported or analysis fails
 *
 * @see solver_ssj Steady-state SSJ backend
 * @see solver_ssj_transient Transient SSJ backend
 */
fun solver_des_analyzer(sn: NetworkStruct, options: SolverOptions, solverDES: SolverDES): DESResult {
    val Tstart = System.nanoTime()

    var res = DESResult()
    val requestedMethod = options.method
    var backendMethod = requestedMethod

    // Validate that this is a Jackson network (open network with exponential service and FCFS queues)
    val isValidJackson = validateJacksonNetwork(sn)
    if (!isValidJackson) {
        throw RuntimeException("solver_des_analyzer: Currently only Jackson queueing networks (possibly multiclass) are supported. " +
                "The model must have: Source, one or more Queues with FCFS and exponential service, Sink, and probabilistic routing.")
    }

    // Check if this is transient analysis (finite timespan)
    val isTransient = options.timespan != null &&
                      options.timespan.size >= 2 &&
                      java.lang.Double.isFinite(options.timespan[1])

    // Map "default" to "ssj" for backend selection
    if (backendMethod == "default") {
        backendMethod = "ssj"
    }

    // Get streaming collector from solver (may be null if not streaming)
    val stream: Collector? = solverDES.getStream()

    // Pick analysis back-end
    when (backendMethod) {
        "ssj" -> {
            res = if (isTransient) {
                solver_ssj_transient(sn, options, stream)
            } else {
                solver_ssj(sn, options, stream)
            }
        }
        else -> {
            throw RuntimeException("solver_des_analyzer:UnknownMethod - Unknown analysis method: ${options.method}")
        }
    }

    // Preserve original requested method name in result
    res.method = requestedMethod
    res.runtime = (System.nanoTime() - Tstart).toDouble() / 1000000000.0
    return res
}

/**
 * @brief Validates that the network structure represents a queueing network suitable for DES.
 *
 * @details This function checks the network topology and configuration to ensure
 * it can be simulated by the DES solver. The validation covers node types,
 * scheduling strategies, and service time distributions.
 *
 * @section validate_networks Supported Network Types
 * - **Open networks**: Source nodes generating arrivals, Queue/Delay nodes, Sink nodes
 * - **Closed networks**: Queue/Delay nodes with fixed populations, no Source/Sink required
 * - **Mixed networks**: Both open and closed classes in the same model
 * - **Petri nets**: Place and Transition nodes for stochastic Petri net modeling
 *
 * @section validate_rules Validation Rules
 * 1. Must have at least one service node (Queue/Delay) or Petri net node (Place/Transition)
 * 2. Open classes require both Source and Sink nodes
 * 3. Queue nodes must use supported scheduling strategies
 * 4. Delay nodes must use INF or EXT scheduling
 * 5. Service processes must be from the supported distribution family
 *
 * @param sn The network structure to validate
 * @return true if valid network for DES, false otherwise
 */
private fun validateJacksonNetwork(sn: NetworkStruct): Boolean {
    var hasSource = false
    var hasSink = false
    var hasServiceNode = false  // Queue or Delay
    var hasPlace = false
    var hasTransition = false
    var hasCache = false

    // Check node types
    for (nodeType in sn.nodetype) {
        when (nodeType) {
            NodeType.Source -> hasSource = true
            NodeType.Sink -> hasSink = true
            NodeType.Queue -> hasServiceNode = true
            NodeType.Delay -> hasServiceNode = true
            NodeType.Place -> hasPlace = true
            NodeType.Transition -> hasTransition = true
            NodeType.Cache -> hasCache = true
            else -> { /* Other types allowed but not required */ }
        }
    }

    // Check if there are open classes (njobs = Inf)
    val hasOpenClasses = (0 until sn.nclasses).any { k ->
        java.lang.Double.isInfinite(sn.njobs.get(k))
    }

    // Must have service nodes OR Place/Transition nodes (Petri net) OR Cache nodes
    val hasPetriNet = hasPlace || hasTransition
    if (!hasServiceNode && !hasPetriNet && !hasCache) {
        return false
    }

    // Open classes require Source and Sink (unless it's a pure closed Petri net)
    if (hasOpenClasses && (!hasSource || !hasSink)) {
        return false
    }

    // Check that queues use FCFS/HOL scheduling and delays use INF
    for ((station, sched) in sn.sched) {
        val stationIdx = station.stationIdx
        val nodeType = sn.nodetype[sn.stationToNode[stationIdx].toInt()]
        when (nodeType) {
            NodeType.Queue -> {
                // FCFS, HOL (priority FCFS), FCFSPRIO, LCFS variants, PS variants, SJF/LJF, SRPT, SEPT/LEPT,
                // PSJF, FB/LAS, LRPT, EDF, and EXT are supported
                val supportedStrategies = setOf(
                    SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRIO,
                    SchedStrategy.LCFS, SchedStrategy.LCFSPR, SchedStrategy.LCFSPI,
                    SchedStrategy.LCFSPRIO, SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO,
                    SchedStrategy.EXT,
                    SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
                    SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO,
                    SchedStrategy.SIRO,
                    SchedStrategy.SJF, SchedStrategy.LJF,
                    SchedStrategy.LEPT, SchedStrategy.SEPT,
                    SchedStrategy.SRPT, SchedStrategy.SRPTPRIO,
                    SchedStrategy.PSJF, SchedStrategy.FB, SchedStrategy.LRPT,  // Size-based preemptive policies
                    SchedStrategy.POLLING,  // Polling scheduling (GATED, EXHAUSTIVE, KLIMITED)
                    SchedStrategy.EDF  // Earliest Deadline First (preemptive)
                )
                if (sched !in supportedStrategies) {
                    return false
                }
            }
            NodeType.Delay -> {
                if (sched != SchedStrategy.INF && sched != SchedStrategy.EXT) {
                    return false
                }
            }
            else -> { /* Other node types don't need scheduling validation */ }
        }
    }

    // Check that service processes are exponential or phase-type
    val allowedProcTypes = setOf(
        ProcessType.DISABLED,
        ProcessType.EXP,
        ProcessType.PH,
        ProcessType.APH,
        ProcessType.HYPEREXP,
        ProcessType.COXIAN,
        ProcessType.COX2,
        ProcessType.ERLANG,
        ProcessType.MAP,
        ProcessType.MMPP2,
        ProcessType.BMAP,
        ProcessType.MMAP,
        ProcessType.ME,
        ProcessType.RAP,
        ProcessType.IMMEDIATE,
        ProcessType.REPLAYER,
        // Additional continuous distributions
        ProcessType.DET,
        ProcessType.UNIFORM,
        ProcessType.GAMMA,
        ProcessType.PARETO,
        // ProcessType.WEIBULL, // TODO: Weibull support disabled - parameter mapping needs investigation
        ProcessType.LOGNORMAL
    )
    for ((_, procMap) in sn.procid) {
        for ((_, procType) in procMap) {
            if (procType !in allowedProcTypes) {
                return false
            }
        }
    }

    // Validate setup/delayoff compatibility
    for (station in sn.stations) {
        if (station is Queue && station.isDelayOffEnabled()) {
            val stationIdx = station.stationIdx
            val sched = sn.sched[station]

            // Setup/delayoff is incompatible with Processor Sharing (PS) variants
            // PS scheduling assumes all jobs share the server simultaneously, while
            // setup/delayoff assumes servers can be turned on/off which conflicts with PS semantics
            val psVariants = setOf(
                SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
                SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO
            )

            if (sched in psVariants) {
                throw RuntimeException(
                    "solver_des_analyzer: Setup/delayoff is incompatible with Processor Sharing (PS) scheduling. " +
                    "Station '${station.name}' uses ${sched} which does not support server on/off transitions."
                )
            }

            // Validate setup and delayoff distributions for all job classes
            for (jobClass in sn.jobclasses) {
                val setupDist = station.getSetupTime(jobClass)
                val delayoffDist = station.getDelayOffTime(jobClass)

                // Setup distribution validation
                if (setupDist != null && !setupDist.isDisabled()) {
                    // Get process type from distribution name
                    // Currently only exponential distributions are fully supported
                    val distName = setupDist.getName().uppercase()
                    if (distName != "EXP" && distName != "EXPONENTIAL") {
                        throw RuntimeException(
                            "solver_des_analyzer: Setup time distribution for station '${station.name}' class '${jobClass.name}' " +
                            "must be Exponential. Distribution '${setupDist.getName()}' is not yet supported."
                        )
                    }
                }

                // Delayoff distribution validation
                if (delayoffDist != null && !delayoffDist.isDisabled()) {
                    val distName = delayoffDist.getName().uppercase()
                    if (distName != "EXP" && distName != "EXPONENTIAL") {
                        throw RuntimeException(
                            "solver_des_analyzer: Delayoff time distribution for station '${station.name}' class '${jobClass.name}' " +
                            "must be Exponential. Distribution '${delayoffDist.getName()}' is not yet supported."
                        )
                    }
                }
            }
        }
    }

    return true
}
