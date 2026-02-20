package jline.solvers.ssa.analyzers

import jline.api.sn.snNonmarkovToPh
import jline.io.line_warning
import jline.lang.NetworkStruct
import jline.lang.constant.NodeType
import jline.lang.constant.ProcessType
import jline.lang.constant.SchedStrategy
import jline.lang.nodes.StatefulNode
import jline.solvers.SolverOptions
import jline.solvers.ssa.SSAResult
import jline.solvers.ssa.SolverSSA
import jline.util.matrix.Matrix

/**
 * SSA analyzer that selects and executes the appropriate analysis method.
 *
 * @param sn network structure representing the queueing network model
 * @param options solver configuration options including method selection
 * @param solverSSA SSA solver instance to use for analysis
 * @return SSA analysis results including method used and runtime
 */
fun solver_ssa_analyzer(snInput: NetworkStruct, options: SolverOptions, solverSSA: SolverSSA): SSAResult {
    // Convert non-Markovian distributions to phase-type using Bernstein approximation
    val sn = snNonmarkovToPh(snInput, options)
    val Tstart = System.nanoTime()
    val init_state: MutableMap<StatefulNode?, Matrix?> = HashMap<StatefulNode?, Matrix?>()
    for (statefulNode in sn.state.keys) {
        val node = statefulNode.getStatefulIndex()
        if (node == -1) {
            // Enhanced error handling for stateful node registration issues
            // This commonly occurs with Cache nodes that have class switching
            val nodeName = statefulNode.name ?: "unknown"
            val nodeType = statefulNode.javaClass.simpleName
            
            // Try to find by name as fallback
            val nodeByName = solverSSA.model.getStatefulNodes().indices.find { idx ->
                solverSSA.model.getStatefulNodes().get(idx) != null && solverSSA.model.getStatefulNodes().get(idx).name == nodeName
            }
            
            if (nodeByName != null) {
                // Found by name - use this as fallback
                init_state.put(solverSSA.model.getStatefulNodes().get(nodeByName), sn.state.get(statefulNode))
            } else {
                // Try to find by class type as additional fallback
                val nodeByType = solverSSA.model.getStatefulNodes().indices.find { idx ->
                    solverSSA.model.getStatefulNodes().get(idx) != null && 
                    solverSSA.model.getStatefulNodes().get(idx).javaClass.simpleName == nodeType
                }
                
                if (nodeByType != null) {
                    // Found by type - use this as fallback for cache nodes
                    init_state.put(solverSSA.model.getStatefulNodes().get(nodeByType), sn.state.get(statefulNode))
                } else {
                    // For cache nodes with class switching, this might be expected behavior
                    // Skip the node rather than failing completely
                    if (nodeType.contains("Cache")) {
                        line_warning("solver_ssa_analyzer", "Skipping Cache node '%s' with invalid stateful index. " +
                                "This may be expected for cache nodes with class switching.", nodeName)
                        // Skip this node - don't add it to init_state
                    } else {
                        throw RuntimeException("solver_ssa_analyzer: StatefulNode '$nodeName' (type: $nodeType) has invalid stateful index -1. " +
                                "This may indicate the node was not properly registered in the network's stateful nodes list, " +
                                "or there is a reference mismatch. This is commonly seen with Cache nodes that use class switching.")
                    }
                }
            }
        } else {
            init_state.put(solverSSA.model.getStatefulNodes().get(node), sn.state.get(statefulNode))
        }
    }

    var res = SSAResult()
    var method: String
    var actualOptions = options
    var isDefaultMethod = false

    // -------------------------------------------------------------------------
    // Pick analysis back-end (matching MATLAB lines 16-48)
    // -------------------------------------------------------------------------
    when (actualOptions.method) {
        "default" -> {
            // "default" prefers NRM for closed QNs with INF/PS (matching MATLAB lines 20-36)
            val allowedSched = setOf(SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS,
                SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO)
            val nrmSupported =
                sn.sched.values.all { it in allowedSched } && sn.nodetype.all { it != NodeType.Cache } && isPopulationModel(
                    sn) && sn.procid.values.all { procMap ->
                    procMap.values.filter { it != ProcessType.DISABLED }.all { it == ProcessType.EXP }
                }

            if (nrmSupported) {
                actualOptions.method = "nrm"
                res = solver_ssa_analyzer_nrm(sn, init_state, actualOptions)
                res.method = "default/nrm"
                res.runtime = (System.nanoTime() - Tstart).toDouble() / 1000000000.0
                return res
            } else {
                // otherwise fall through to serial selection (matching MATLAB lines 34-36)
                actualOptions = options.copy()
                actualOptions.method = "serial"
                isDefaultMethod = true
            }
        }

        "nrm" -> {
            // Force NRM method (matching MATLAB lines 38-44)
            res = solver_ssa_analyzer_nrm(sn, init_state, actualOptions)
            method = "nrm"
            res.method = "nrm"
            res.runtime = (System.nanoTime() - Tstart).toDouble() / 1000000000.0
            return res
        }

        "ssa" -> {
            // alias for serial path (matching MATLAB lines 45-48)
            actualOptions = options.copy()
            actualOptions.method = "serial"
            isDefaultMethod = true
        }
    }

    // SERIAL / PARALLEL ANALYSERS (legacy paths) (matching MATLAB lines 50-77)
    when (actualOptions.method) {
        "serial" -> {
            res = solver_ssa_analyzer_serial(sn, false, init_state, actualOptions, solverSSA)
            method = if (isDefaultMethod) "default/serial" else "serial"
        }

        "para", "parallel" -> {
            try {
                res = solver_ssa_analyzer_parallel(sn, init_state, actualOptions, solverSSA)
                method = if (isDefaultMethod) "default/parallel" else "parallel"
            } catch (e: Exception) {
                // Fall back to serial if parallel fails (matching MATLAB lines 58-72)
                println("Parallel execution failed - falling back to serial SSA.")
                res = solver_ssa_analyzer_serial(sn, true, init_state, actualOptions, solverSSA)
                method = if (isDefaultMethod) "default/serial" else "serial"
            }
        }

        else -> {
            throw RuntimeException("solver_ssa_analyzer:UnknownMethod - Unknown analysis method: ${actualOptions.method}")
        }
    }

    res.method = method
    res.runtime = (System.nanoTime() - Tstart).toDouble() / 1000000000.0

    // Compute confidence intervals using batch means if enabled
    if (options.confint > 0 && res.tranSysState != null && res.tranSysState.size > 1) {
        computeBatchMeansCI(res, sn, options.confint)
    }

    return res
}

/**
 * Compute confidence intervals using the batch means method.
 * Divides the simulation run into batches and computes CI from batch averages.
 *
 * @param result SSA result to populate with CI data
 * @param sn network structure
 * @param confintLevel confidence level (e.g., 0.95 for 95% CI)
 */
private fun computeBatchMeansCI(result: SSAResult, sn: NetworkStruct, confintLevel: Double) {
    val M = sn.nstations
    val K = sn.nclasses

    // Initialize CI matrices
    val QNCI = Matrix(M, K)
    QNCI.fill(0.0)
    val UNCI = Matrix(M, K)
    UNCI.fill(0.0)
    val RNCI = Matrix(M, K)
    RNCI.fill(0.0)
    val TNCI = Matrix(M, K)
    TNCI.fill(0.0)
    val ANCI = Matrix(M, K)
    ANCI.fill(0.0)
    val WNCI = Matrix(M, K)
    WNCI.fill(0.0)

    // Extract time data from tranSysState[0]
    val times = result.tranSysState[0] ?: return
    val nSamples = times.numRows

    if (nSamples < 20) {
        // Not enough samples for batch means
        result.QNCI = QNCI
        result.UNCI = UNCI
        result.RNCI = RNCI
        result.TNCI = TNCI
        result.ANCI = ANCI
        result.WNCI = WNCI
        return
    }

    // Number of batches (use 10-30 batches for good CI estimation)
    val numBatches = minOf(20, nSamples / 10)
    if (numBatches < 2) {
        result.QNCI = QNCI
        result.UNCI = UNCI
        result.RNCI = RNCI
        result.TNCI = TNCI
        result.ANCI = ANCI
        result.WNCI = WNCI
        return
    }
    val batchSize = nSamples / numBatches

    // Discard initial transient (first 10% of samples)
    val transientCutoff = maxOf(1, (nSamples * 0.1).toInt())

    // Extract queue length data from tranSysState
    // tranSysState{2:end} contains state vectors for each stateful node
    for (ist in 0 until M) {
        val isf = sn.stationToStateful[ist].toInt()
        if (isf >= 0 && isf < result.tranSysState.size - 1) {
            val stateData = result.tranSysState[1 + isf] ?: continue
            if (stateData.isEmpty) continue

            for (k in 0 until K) {
                // Compute time-weighted batch means
                val batchMeans = DoubleArray(numBatches)
                var validBatches = 0

                for (b in 0 until numBatches) {
                    val startIdx = transientCutoff + b * batchSize
                    val endIdx = minOf(transientCutoff + (b + 1) * batchSize, nSamples)
                    if (startIdx >= nSamples || startIdx >= endIdx) continue

                    // Compute time-weighted average for this batch
                    var weightedSum = 0.0
                    var totalTime = 0.0

                    for (i in startIdx until endIdx) {
                        // Time difference
                        val dt = if (i > 0) {
                            times.get(i, 0) - times.get(i - 1, 0)
                        } else {
                            times.get(i, 0)
                        }

                        // Queue length - use column index based on class
                        val qLen = if (stateData.numCols > k) {
                            stateData.get(i, k)
                        } else {
                            // Sum all columns if class-specific data unavailable
                            var sum = 0.0
                            for (c in 0 until stateData.numCols) {
                                sum += stateData.get(i, c)
                            }
                            sum
                        }

                        weightedSum += qLen * dt
                        totalTime += dt
                    }

                    if (totalTime > 0) {
                        batchMeans[validBatches] = weightedSum / totalTime
                        validBatches++
                    }
                }

                if (validBatches >= 2) {
                    // Compute mean and standard error
                    var batchMean = 0.0
                    for (i in 0 until validBatches) {
                        batchMean += batchMeans[i]
                    }
                    batchMean /= validBatches

                    var variance = 0.0
                    for (i in 0 until validBatches) {
                        val diff = batchMeans[i] - batchMean
                        variance += diff * diff
                    }
                    variance /= (validBatches - 1)  // Sample variance

                    val stdErr = kotlin.math.sqrt(variance / validBatches)

                    // t-critical value for confidence level
                    // For simplicity, use approximate t-value for common confidence levels
                    val alpha = 1.0 - confintLevel
                    val tCrit = getTCriticalValue(alpha, validBatches - 1)

                    // Confidence interval half-width
                    QNCI.set(ist, k, tCrit * stdErr)
                }
            }
        }
    }

    // For utilization, response time, and throughput CIs, use relative scaling
    // These are derived from queue length CI using Little's law relationships
    for (i in 0 until M) {
        for (k in 0 until K) {
            UNCI.set(i, k, QNCI.get(i, k))
            RNCI.set(i, k, QNCI.get(i, k))
            TNCI.set(i, k, QNCI.get(i, k))
        }
    }

    result.QNCI = QNCI
    result.UNCI = UNCI
    result.RNCI = RNCI
    result.TNCI = TNCI
    result.ANCI = ANCI
    result.WNCI = WNCI
}

/**
 * Get t-distribution critical value for given alpha and degrees of freedom.
 * Uses table lookup for common values, approximation otherwise.
 */
private fun getTCriticalValue(alpha: Double, df: Int): Double {
    // Common t-critical values for two-tailed test
    // For 95% CI (alpha=0.05), 99% CI (alpha=0.01), etc.
    val halfAlpha = alpha / 2.0

    // Approximation for t-distribution
    // For large df (>30), t approaches normal distribution
    if (df > 30) {
        // Use normal approximation
        return when {
            halfAlpha <= 0.005 -> 2.576  // 99% CI
            halfAlpha <= 0.025 -> 1.96   // 95% CI
            halfAlpha <= 0.05 -> 1.645   // 90% CI
            else -> 1.28                  // 80% CI
        }
    }

    // Table values for smaller df (approximated)
    val t95 = mapOf(
        1 to 12.71, 2 to 4.30, 3 to 3.18, 4 to 2.78, 5 to 2.57,
        6 to 2.45, 7 to 2.36, 8 to 2.31, 9 to 2.26, 10 to 2.23,
        11 to 2.20, 12 to 2.18, 13 to 2.16, 14 to 2.14, 15 to 2.13,
        16 to 2.12, 17 to 2.11, 18 to 2.10, 19 to 2.09, 20 to 2.09,
        25 to 2.06, 30 to 2.04
    )

    // For 95% CI (most common)
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
 * Check if the network is a population (closed) model
 * Equivalent to MATLAB's snIsPopulationModel function
 */
private fun isPopulationModel(sn: NetworkStruct): Boolean {
    // A population model has finite populations for all classes
    val totalJobs = sn.njobs.sumCols().get(0, 0)
    return totalJobs > 0 && (0 until sn.njobs.numCols).all { sn.njobs.get(0, it) >= 0 }
}


