/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.fluid.analyzers

import jline.GlobalConstants
import jline.lang.NetworkStruct
import jline.lang.constant.JobClassType
import jline.lang.constant.NodeType
import jline.lib.butools.queues.fluFluQueue
import jline.solvers.SolverOptions
import jline.solvers.SolverResult
import jline.util.matrix.Matrix
import kotlin.math.min

/**
 * MFQ (Markovian Fluid Queue) analyzer for single-queue open systems.
 *
 * Uses BUTools FluFluQueue function to compute exact steady-state queue length
 * and sojourn time moments for single-queue networks with phase-type arrivals
 * and service distributions.
 *
 * Applicability:
 * - Single-queue open system: Source -> Queue -> Sink
 * - Single-server (c=1) or infinite-server (c=Inf)
 * - Exponential or phase-type distributions
 *
 * References:
 * Horvath G, Telek M, "Sojourn times in fluid queues with independent and
 * dependent input and output processes", PERFORMANCE EVALUATION 79: pp. 160-181, 2014.
 */
class MFQAnalyzer : FluidAnalyzer {
    private var xvecIt: Matrix? = null

    override fun analyze(sn: NetworkStruct, options: SolverOptions, result: SolverResult) {
        val M = sn.nstations
        val K = sn.nclasses

        // Initialize result matrices
        result.QN = Matrix.zeros(M, K)
        result.UN = Matrix.zeros(M, K)
        result.RN = Matrix.zeros(M, K)
        result.TN = Matrix.zeros(M, K)
        result.WN = Matrix.zeros(M, K)
        result.AN = Matrix.zeros(M, K)

        // Validate topology
        val topologyInfo = validateTopology(sn)
        if (!topologyInfo.isValid) {
            throw RuntimeException("MFQ requires single-queue topology: ${topologyInfo.errorMsg}")
        }

        val sourceStation = topologyInfo.sourceStation
        val queueStation = topologyInfo.queueStation

        // Process each open class
        for (k in 0 until K) {
            if (sn.njobs[0, k] == Double.POSITIVE_INFINITY) {
                // Open class
                val classResult = analyzeOpenClass(sn, k, sourceStation, queueStation, options)

                result.QN[queueStation, k] = classResult.meanQueueLength
                result.RN[queueStation, k] = classResult.meanResponseTime
                result.TN[queueStation, k] = classResult.throughput
                result.UN[queueStation, k] = classResult.utilization

                // Throughput at source equals queue throughput
                result.TN[sourceStation, k] = classResult.throughput
            }
        }

        // Create dummy state vector for compatibility
        xvecIt = Matrix.zeros(1, 1)
    }

    override fun getXVecIt(): Matrix? = xvecIt

    /**
     * Validates that the network has a single-queue topology suitable for MFQ.
     */
    private fun validateTopology(sn: NetworkStruct): TopologyInfo {
        val info = TopologyInfo()

        // Check for open classes
        var hasOpenClass = false
        for (k in 0 until sn.nclasses) {
            if (sn.njobs[0, k] == Double.POSITIVE_INFINITY) {
                hasOpenClass = true
                break
            }
        }
        if (!hasOpenClass) {
            info.errorMsg = "Not an open model - all classes are closed"
            return info
        }

        // Find Source, Queue, and Sink nodes
        var sourceNode = -1
        var queueNode = -1
        var sinkNode = -1
        var sourceCount = 0
        var queueCount = 0
        var sinkCount = 0

        for (i in 0 until sn.nnodes) {
            when (sn.nodetype[i]) {
                NodeType.Source -> {
                    sourceNode = i
                    sourceCount++
                }
                NodeType.Queue -> {
                    queueNode = i
                    queueCount++
                }
                NodeType.Sink -> {
                    sinkNode = i
                    sinkCount++
                }
                else -> {}
            }
        }

        // Validate counts
        when {
            sourceCount == 0 -> {
                info.errorMsg = "No source node found"
                return info
            }
            sourceCount > 1 -> {
                info.errorMsg = "Multiple source nodes found ($sourceCount)"
                return info
            }
            sinkCount == 0 -> {
                info.errorMsg = "No sink node found"
                return info
            }
            sinkCount > 1 -> {
                info.errorMsg = "Multiple sink nodes found ($sinkCount)"
                return info
            }
            queueCount == 0 -> {
                info.errorMsg = "No queue node found"
                return info
            }
            queueCount > 1 -> {
                info.errorMsg = "Multiple queue nodes found ($queueCount) - MFQ supports single queue only"
                return info
            }
        }

        // Get station indices
        val sourceStation = sn.nodeToStation[sourceNode].toInt()
        val queueStation = sn.nodeToStation[queueNode].toInt()

        // Check server count (single-server or infinite-server)
        val c = sn.nservers[queueStation, 0]
        if (c > 1 && c < Double.POSITIVE_INFINITY) {
            info.errorMsg = "Multi-server queue (c=$c) not supported - MFQ requires c=1 or c=Inf"
            return info
        }

        info.isValid = true
        info.sourceNode = sourceNode
        info.queueNode = queueNode
        info.sinkNode = sinkNode
        info.sourceStation = sourceStation
        info.queueStation = queueStation
        return info
    }

    /**
     * Analyzes an open class using FluFluQueue.
     */
    private fun analyzeOpenClass(
        sn: NetworkStruct,
        classIdx: Int,
        sourceStation: Int,
        queueStation: Int,
        options: SolverOptions
    ): ClassResult {
        // Extract arrival rate
        val lambda = sn.rates[sourceStation, classIdx]
        if (lambda.isNaN() || lambda <= 0) {
            throw RuntimeException("No valid arrival rate for class $classIdx at source")
        }

        // Extract service rate (mean)
        val mu = sn.rates[queueStation, classIdx]
        if (mu.isNaN() || mu <= 0) {
            throw RuntimeException("No valid service rate for class $classIdx at queue")
        }

        // Check stability
        val rho = lambda / mu
        if (rho >= 1.0) {
            // Unstable system
            return ClassResult(
                meanQueueLength = Double.POSITIVE_INFINITY,
                meanResponseTime = Double.POSITIVE_INFINITY,
                throughput = lambda,
                utilization = 1.0
            )
        }

        // Extract process matrices
        val arrivalParams = extractArrivalProcess(sn, sourceStation, classIdx, lambda)
        val serviceParams = extractServiceProcess(sn, queueStation, classIdx, mu)

        // Check if simple M/M/1 case
        val isSimpleExponential = arrivalParams.isSimple && serviceParams.isSimple

        val flMoms: DoubleArray
        val stMoms: DoubleArray

        if (isSimpleExponential) {
            // Use analytical M/M/1 formulas
            val meanL = rho / (1 - rho)
            val varL = rho / ((1 - rho) * (1 - rho))
            val meanW = 1.0 / (mu - lambda)
            val varW = 1.0 / ((mu - lambda) * (mu - lambda))

            flMoms = doubleArrayOf(meanL, meanL * meanL + varL)
            stMoms = doubleArrayOf(meanW, meanW * meanW + varW)
        } else {
            // Use FluFluQueue for general case
            val prec = maxOf(options.tol, 1e-14)
            val srv0stop = true  // Work-conserving behavior

            val fluResult = fluFluQueue(
                Qin = arrivalParams.Q,
                Rin = arrivalParams.R,
                Qout = serviceParams.Q,
                Rout = serviceParams.R,
                srv0stop = srv0stop,
                numFluidMoments = 2,
                numSojournMoments = 2,
                prec = prec
            )

            flMoms = fluResult.fluidMoments ?: throw RuntimeException("FluFluQueue did not return fluid moments")
            stMoms = fluResult.sojournMoments ?: throw RuntimeException("FluFluQueue did not return sojourn moments")
        }

        // Map moments to LINE metrics
        val meanQueueLength = flMoms[0]
        val meanResponseTime = stMoms[0]

        // Throughput from Little's Law: X = L / W
        val throughput = if (meanResponseTime > GlobalConstants.FineTol) {
            meanQueueLength / meanResponseTime
        } else {
            lambda
        }

        // Utilization: U = lambda * E[S]
        val utilization = min(1.0, lambda / mu)

        return ClassResult(
            meanQueueLength = meanQueueLength,
            meanResponseTime = meanResponseTime,
            throughput = throughput,
            utilization = utilization
        )
    }

    /**
     * Extracts arrival process parameters (Q, R matrices) from network structure.
     */
    private fun extractArrivalProcess(
        sn: NetworkStruct,
        stationIdx: Int,
        classIdx: Int,
        lambda: Double
    ): ProcessParams {
        // Get the station and job class objects
        val station = sn.stations[stationIdx]
        val jobClass = sn.jobclasses[classIdx]

        // Check if proc data is available
        val proc = sn.proc[station]?.get(jobClass)

        if (proc == null || proc.isEmpty()) {
            // Simple Poisson arrival (1-state)
            return ProcessParams(
                Q = Matrix.zeros(1, 1),
                R = Matrix.singleton(lambda),
                isSimple = true
            )
        }

        // Check number of phases
        val nPhases = proc[0].numRows
        if (nPhases == 1) {
            // Single state: Poisson arrival
            return ProcessParams(
                Q = Matrix.zeros(1, 1),
                R = Matrix.singleton(lambda),
                isSimple = true
            )
        }

        // Multi-phase: Extract D0, D1 and convert to (Q, R) format
        // D0 + D1 = Q (generator), sum of D1 row = R (fluid rate)
        val D0 = proc[0]
        val D1 = proc[1]

        val Q = D0.add(1.0, D1)
        val R = Matrix.zeros(nPhases, nPhases)
        for (i in 0 until nPhases) {
            var rowSum = 0.0
            for (j in 0 until D1.numCols) {
                rowSum += D1[i, j]
            }
            R[i, i] = rowSum
        }

        return ProcessParams(Q = Q, R = R, isSimple = false)
    }

    /**
     * Extracts service process parameters (Q, R matrices) from network structure.
     */
    private fun extractServiceProcess(
        sn: NetworkStruct,
        stationIdx: Int,
        classIdx: Int,
        mu: Double
    ): ProcessParams {
        // Get the station and job class objects
        val station = sn.stations[stationIdx]
        val jobClass = sn.jobclasses[classIdx]

        // Check if proc data is available
        val proc = sn.proc[station]?.get(jobClass)

        if (proc == null || proc.isEmpty()) {
            // Simple exponential service (1-phase)
            return ProcessParams(
                Q = Matrix.zeros(1, 1),
                R = Matrix.singleton(mu),
                isSimple = true
            )
        }

        // Check number of phases
        val nPhases = proc[0].numRows
        if (nPhases == 1) {
            // Single phase: exponential service
            return ProcessParams(
                Q = Matrix.zeros(1, 1),
                R = Matrix.singleton(mu),
                isSimple = true
            )
        }

        // Multi-phase: Extract D0, D1 and convert to (Q, R) format
        val D0 = proc[0]
        val D1 = proc[1]

        val Q = D0.add(1.0, D1)
        val R = Matrix.zeros(nPhases, nPhases)
        for (i in 0 until nPhases) {
            var rowSum = 0.0
            for (j in 0 until D1.numCols) {
                rowSum += D1[i, j]
            }
            R[i, i] = rowSum
        }

        return ProcessParams(Q = Q, R = R, isSimple = false)
    }

    /**
     * Internal data class for topology validation results.
     */
    private data class TopologyInfo(
        var isValid: Boolean = false,
        var errorMsg: String = "",
        var sourceNode: Int = -1,
        var queueNode: Int = -1,
        var sinkNode: Int = -1,
        var sourceStation: Int = -1,
        var queueStation: Int = -1
    )

    /**
     * Internal data class for process parameters.
     */
    private data class ProcessParams(
        val Q: Matrix,
        val R: Matrix,
        val isSimple: Boolean
    )

    /**
     * Internal data class for class analysis results.
     */
    private data class ClassResult(
        val meanQueueLength: Double,
        val meanResponseTime: Double,
        val throughput: Double,
        val utilization: Double
    )
}
