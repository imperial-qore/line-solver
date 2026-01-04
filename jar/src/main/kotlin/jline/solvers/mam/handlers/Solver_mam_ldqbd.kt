/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam.handlers

import jline.io.line_error
import jline.io.mfilename
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.SchedStrategy
import jline.solvers.SolverOptions
import jline.solvers.mam.MAMResult
import jline.util.matrix.Matrix
import jline.api.mam.ldqbd
import jline.api.mam.LdqbdOptions
import jline.api.mam.map_pie
import jline.api.mam.map_mean
import jline.api.mam.map_scale
import kotlin.math.min

/**
 * Solver for single-class closed queueing networks using Level-Dependent QBD.
 *
 * Uses Level-Dependent Quasi-Birth-Death (LD-QBD) process to compute
 * exact performance metrics for single-class closed queueing networks
 * consisting of a Delay (infinite server) and a Queue (FCFS).
 *
 * The LD-QBD approach models the system where:
 *   - Level n = number of jobs at the queue (0 <= n <= N)
 *   - Jobs at delay = N - n
 *   - Transition rates depend on the current level
 *
 * Supports PH-type service distributions (Exp, Erlang, HyperExp, etc.)
 */
fun solver_mam_ldqbd(sn: NetworkStruct, options: SolverOptions): MAMResult {
    val M = sn.nstations
    val K = sn.nclasses

    // Check: single-class closed network
    if (K != 1) {
        line_error(mfilename(object {}), "LDQBD method requires a single-class model.")
        return createEmptyResult(M, K)
    }

    val N = sn.njobs[0, 0].toInt()
    if (!java.lang.Double.isFinite(N.toDouble()) || N <= 0) {
        line_error(mfilename(object {}), "LDQBD method requires a closed model with finite population.")
        return createEmptyResult(M, K)
    }

    // Check: must have exactly one delay and one queue
    var nDelay = 0
    var nQueue = 0
    var delayIdx = -1
    var queueIdx = -1

    for (i in 0 until M) {
        val sched = sn.sched[sn.stations[i]]
        if (sched == SchedStrategy.INF) {
            nDelay++
            delayIdx = i
        } else if (sched == SchedStrategy.FCFS) {
            nQueue++
            queueIdx = i
        }
    }

    if (nDelay != 1 || nQueue != 1 || M != 2) {
        line_error(mfilename(object {}), "LDQBD method requires exactly one Delay and one Queue station.")
        return createEmptyResult(M, K)
    }

    // Get service parameters
    val rates = sn.rates
    val nservers = sn.nservers

    // Delay station rate (infinite server)
    val lambda_d = rates[delayIdx, 0]

    // Get station and job class objects for proc access
    val queueStation = sn.stations[queueIdx]
    val jobClass = sn.jobclasses[0]

    // Queue service process
    val PH_queue = sn.proc[queueStation]?.get(jobClass)
        ?: throw RuntimeException("No service process for queue station")
    val nServers = nservers[queueIdx, 0].toInt()

    // Check if queue service is exponential (1x1 matrix) or PH
    val D0 = PH_queue[0]
        ?: throw RuntimeException("No D0 matrix in service process")
    val isExponential = D0.numRows == 1 && D0.numCols == 1
    val mu: Double
    val nPhases: Int

    if (isExponential) {
        mu = -D0[0, 0]
        nPhases = 1
    } else {
        nPhases = D0.numRows
        mu = 1.0 / map_mean(PH_queue)
    }

    // Build LD-QBD matrices
    val Q0 = mutableListOf<Matrix>()
    val Q1 = mutableListOf<Matrix>()
    val Q2 = mutableListOf<Matrix>()

    if (isExponential) {
        // Exponential service case (scalar matrices)
        for (n in 0 until N) {
            // Arrival rate: (N-n) * lambda_d
            val arrRate = Matrix(1, 1)
            arrRate[0, 0] = (N - n) * lambda_d
            Q0.add(arrRate)
        }

        for (n in 0..N) {
            val arrivalRate = (N - n) * lambda_d
            val departureRate = if (n > 0) min(n, nServers) * mu else 0.0
            val local = Matrix(1, 1)
            local[0, 0] = -(arrivalRate + departureRate)
            Q1.add(local)
        }

        for (n in 1..N) {
            // Service rate: min(n, nServers) * mu
            val depRate = Matrix(1, 1)
            depRate[0, 0] = min(n, nServers) * mu
            Q2.add(depRate)
        }
    } else {
        // PH-type service case (matrix-valued)
        val D1 = PH_queue[1]
            ?: throw RuntimeException("No D1 matrix in PH service process")
        val alpha = map_pie(PH_queue)

        // Level 0 -> 1: arrival starts service in some phase
        Q0.add(alpha.scale(N * lambda_d))

        // Level n -> n+1 (n >= 1): arrival, preserve phase
        for (n in 1 until N) {
            Q0.add(Matrix.eye(nPhases).scale((N - n) * lambda_d))
        }

        // Level 0: 1x1 (no service, only arrivals)
        val Q1_0 = Matrix(1, 1)
        Q1_0[0, 0] = -(N * lambda_d)
        Q1.add(Q1_0)

        // Level n >= 1: phase transitions within level
        for (n in 1..N) {
            val arrivalRate = (N - n) * lambda_d
            val c_n = min(n, nServers)
            val local = D0.scale(c_n.toDouble()).sub(Matrix.eye(nPhases).scale(arrivalRate))
            Q1.add(local)
        }

        // Level 1 -> 0: service completion, go to empty state
        val c_1 = min(1, nServers)
        val Q2_1 = D1.scale(c_1.toDouble()).mult(Matrix.ones(nPhases, 1))
        Q2.add(Q2_1)

        // Level n -> n-1 (n >= 2): service completion, next job starts
        for (n in 2..N) {
            val c_n = min(n, nServers)
            Q2.add(D1.scale(c_n.toDouble()))
        }
    }

    // Solve LD-QBD
    val ldqbdOptions = LdqbdOptions(
        epsilon = options.tol,
        maxIter = options.iter_max,
        verbose = false
    )

    val ldqbdResult = ldqbd(Q0, Q1, Q2, ldqbdOptions)
    val pi_ldqbd = ldqbdResult.pi

    // Compute performance metrics from steady-state distribution
    var mean_queue = 0.0
    for (n in 0..N) {
        mean_queue += n * pi_ldqbd[0, n]
    }

    val mean_delay = N - mean_queue

    // System throughput X = mean_delay * lambda_d (flow balance)
    val X = mean_delay * lambda_d

    // Mean service time at queue
    val mean_service = if (isExponential) 1.0 / mu else map_mean(PH_queue)

    // Queue utilization
    val util_queue: Double
    if (nServers == 1) {
        util_queue = 1.0 - pi_ldqbd[0, 0]
    } else {
        var u = 0.0
        for (n in 1..N) {
            u += (min(n, nServers).toDouble() / nServers) * pi_ldqbd[0, n]
        }
        util_queue = u
    }

    // Response time at queue (using Little's law: R = Q/X)
    val R_queue = if (X > 0) mean_queue / X else 0.0

    // Response time at delay (constant for infinite server)
    val R_delay = 1.0 / lambda_d

    // Populate output matrices
    val result = MAMResult()
    result.QN = Matrix(M, K)
    result.UN = Matrix(M, K)
    result.RN = Matrix(M, K)
    result.TN = Matrix(M, K)
    result.CN = Matrix(1, K)
    result.XN = Matrix(1, K)

    // Delay station metrics
    result.QN[delayIdx, 0] = mean_delay
    result.UN[delayIdx, 0] = mean_delay  // For infinite server, U = Q
    result.RN[delayIdx, 0] = R_delay
    result.TN[delayIdx, 0] = X

    // Queue station metrics
    result.QN[queueIdx, 0] = mean_queue
    result.UN[queueIdx, 0] = util_queue / nServers  // Per-server utilization
    result.RN[queueIdx, 0] = R_queue
    result.TN[queueIdx, 0] = X

    // System-level metrics
    result.XN[0, 0] = X
    result.CN[0, 0] = R_delay + R_queue  // Cycle time

    result.iter = 1  // LDQBD is a direct method
    result.method = "ldqbd"

    return result
}

private fun createEmptyResult(M: Int, K: Int): MAMResult {
    val result = MAMResult()
    result.QN = Matrix(M, K)
    result.UN = Matrix(M, K)
    result.RN = Matrix(M, K)
    result.TN = Matrix(M, K)
    result.CN = Matrix(1, K)
    result.XN = Matrix(1, K)
    result.iter = 0
    result.method = "ldqbd"
    return result
}
