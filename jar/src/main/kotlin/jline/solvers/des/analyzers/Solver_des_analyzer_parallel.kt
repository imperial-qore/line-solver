/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.des.analyzers

import jline.lang.NetworkStruct
import jline.solvers.des.DESOptions
import jline.solvers.des.DESResult
import jline.solvers.des.SolverDES
import jline.solvers.des.handlers.solver_ssj
import jline.util.RandomManager
import jline.util.matrix.Matrix
import org.apache.commons.math3.distribution.TDistribution
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.TimeUnit

/**
 * Parallel replication analyzer for DES solver.
 *
 * Runs multiple independent simulation replications in parallel and aggregates
 * results using cross-replication statistics for confidence intervals.
 *
 * @param sn Network structure
 * @param options DES solver options (including replications and numThreads)
 * @param solverDES DES solver instance with thread pool
 * @return Aggregated DESResult with cross-replication confidence intervals
 */
fun solver_des_analyzer_parallel(
    sn: NetworkStruct,
    options: DESOptions,
    solverDES: SolverDES
): DESResult {
    val Tstart = System.nanoTime()
    val numReplications = options.replications
    val numThreads = solverDES.numThreads

    val M = sn.nstations
    val K = sn.nclasses

    // Thread-safe storage for per-replication results
    val QNs = ConcurrentHashMap<Int, Matrix>()
    val UNs = ConcurrentHashMap<Int, Matrix>()
    val RNs = ConcurrentHashMap<Int, Matrix>()
    val TNs = ConcurrentHashMap<Int, Matrix>()
    val ANs = ConcurrentHashMap<Int, Matrix>()
    val CNs = ConcurrentHashMap<Int, Matrix>()
    val XNs = ConcurrentHashMap<Int, Matrix>()

    // Impatience metrics storage
    val renegedCounts = ConcurrentHashMap<Int, Matrix>()
    val avgRenegingWaits = ConcurrentHashMap<Int, Matrix>()
    val renegingRates = ConcurrentHashMap<Int, Matrix>()
    val balkedCounts = ConcurrentHashMap<Int, Matrix>()
    val balkingProbs = ConcurrentHashMap<Int, Matrix>()
    val retriedCounts = ConcurrentHashMap<Int, Matrix>()
    val retrialDroppedCounts = ConcurrentHashMap<Int, Matrix>()
    val avgOrbitSizes = ConcurrentHashMap<Int, Matrix>()

    // Divide samples across replications (ensure each replication gets reasonable amount)
    val samplesPerReplication = Math.max(1000, (options.samples + numReplications - 1) / numReplications)

    // Submit replication tasks
    for (repIdx in 0 until numReplications) {
        val replicationIndex = repIdx
        solverDES.threadPool.submit(Runnable {
            // Create per-replication RNG stream for reproducibility
            val repRandom = RandomManager.getParallelRandom("DES_PARALLEL", replicationIndex)

            // Create per-replication options copy
            val repOptions = options.copy()
            repOptions.samples = samplesPerReplication
            // Use deterministic per-replication seed based on master seed
            repOptions.seed = options.seed + replicationIndex

            // Run single replication (no streaming in parallel mode)
            val repResult = solver_ssj(sn, repOptions, null)

            // Store results
            QNs[replicationIndex] = repResult.QN
            UNs[replicationIndex] = repResult.UN
            RNs[replicationIndex] = repResult.RN
            TNs[replicationIndex] = repResult.TN
            if (repResult.AN != null) {
                ANs[replicationIndex] = repResult.AN
            } else {
                ANs[replicationIndex] = Matrix(M, K)
            }
            CNs[replicationIndex] = repResult.CN
            XNs[replicationIndex] = repResult.XN

            // Store impatience metrics if available
            if (repResult.renegedCustomers != null) {
                renegedCounts[replicationIndex] = repResult.renegedCustomers
            }
            if (repResult.avgRenegingWaitTime != null) {
                avgRenegingWaits[replicationIndex] = repResult.avgRenegingWaitTime
            }
            if (repResult.renegingRate != null) {
                renegingRates[replicationIndex] = repResult.renegingRate
            }
            if (repResult.balkedCustomers != null) {
                balkedCounts[replicationIndex] = repResult.balkedCustomers
            }
            if (repResult.balkingProbability != null) {
                balkingProbs[replicationIndex] = repResult.balkingProbability
            }
            if (repResult.retriedCustomers != null) {
                retriedCounts[replicationIndex] = repResult.retriedCustomers
            }
            if (repResult.retrialDropped != null) {
                retrialDroppedCounts[replicationIndex] = repResult.retrialDropped
            }
            if (repResult.avgOrbitSize != null) {
                avgOrbitSizes[replicationIndex] = repResult.avgOrbitSize
            }
        })
    }

    // Wait for all replications to complete
    solverDES.threadPool.shutdown()
    try {
        solverDES.threadPool.awaitTermination(600, TimeUnit.SECONDS)
    } catch (e: InterruptedException) {
        throw RuntimeException("Parallel DES execution interrupted", e)
    }

    // Aggregate results - compute grand mean across all replications
    val QN = averageMatrices(QNs.values, numReplications)
    val UN = averageMatrices(UNs.values, numReplications)
    val RN = averageMatrices(RNs.values, numReplications)
    val TN = averageMatrices(TNs.values, numReplications)
    val AN = averageMatrices(ANs.values, numReplications)
    val CN = averageMatrices(CNs.values, numReplications)
    val XN = averageMatrices(XNs.values, numReplications)

    // Compute cross-replication confidence intervals
    val alpha = 1.0 - options.confint
    val tCrit = if (numReplications > 1) {
        val tDist = TDistribution((numReplications - 1).toDouble())
        tDist.inverseCumulativeProbability(1.0 - alpha / 2.0)
    } else {
        0.0  // No CI with single replication
    }

    val QNCI = computeCrossReplicationCI(QNs.values.toList(), QN, tCrit, numReplications)
    val UNCI = computeCrossReplicationCI(UNs.values.toList(), UN, tCrit, numReplications)
    val RNCI = computeCrossReplicationCI(RNs.values.toList(), RN, tCrit, numReplications)
    val TNCI = computeCrossReplicationCI(TNs.values.toList(), TN, tCrit, numReplications)
    val ANCI = computeCrossReplicationCI(ANs.values.toList(), AN, tCrit, numReplications)
    val WNCI = RNCI.copy()  // Residence time CI same as response time

    // Compute relative precision
    val QNRelPrec = computeRelativePrecision(QNCI, QN)
    val UNRelPrec = computeRelativePrecision(UNCI, UN)
    val RNRelPrec = computeRelativePrecision(RNCI, RN)
    val TNRelPrec = computeRelativePrecision(TNCI, TN)

    // Build result
    val result = DESResult()
    result.QN = QN
    result.UN = UN
    result.RN = RN
    result.TN = TN
    result.AN = AN
    result.CN = CN
    result.XN = XN
    result.sn = sn
    result.method = "parallel"
    result.runtime = (System.nanoTime() - Tstart).toDouble() / 1_000_000_000.0

    // Set CI data
    result.QNCI = QNCI
    result.UNCI = UNCI
    result.RNCI = RNCI
    result.TNCI = TNCI
    result.ANCI = ANCI
    result.WNCI = WNCI
    result.QNRelPrec = QNRelPrec
    result.UNRelPrec = UNRelPrec
    result.RNRelPrec = RNRelPrec
    result.TNRelPrec = TNRelPrec

    // Convergence: consider converged if all relative precisions < tolerance
    result.converged = checkConvergence(QNRelPrec, UNRelPrec, RNRelPrec, TNRelPrec, options.cnvgtol)
    result.stoppingReason = if (result.converged) "convergence" else "max_replications"
    result.convergenceBatches = numReplications

    // Aggregate impatience metrics
    if (renegedCounts.isNotEmpty()) {
        result.renegedCustomers = averageMatrices(renegedCounts.values, numReplications)
    }
    if (avgRenegingWaits.isNotEmpty()) {
        result.avgRenegingWaitTime = averageMatrices(avgRenegingWaits.values, numReplications)
    }
    if (renegingRates.isNotEmpty()) {
        result.renegingRate = averageMatrices(renegingRates.values, numReplications)
    }
    if (balkedCounts.isNotEmpty()) {
        result.balkedCustomers = averageMatrices(balkedCounts.values, numReplications)
    }
    if (balkingProbs.isNotEmpty()) {
        result.balkingProbability = averageMatrices(balkingProbs.values, numReplications)
    }
    if (retriedCounts.isNotEmpty()) {
        result.retriedCustomers = averageMatrices(retriedCounts.values, numReplications)
    }
    if (retrialDroppedCounts.isNotEmpty()) {
        result.retrialDropped = averageMatrices(retrialDroppedCounts.values, numReplications)
    }
    if (avgOrbitSizes.isNotEmpty()) {
        result.avgOrbitSize = averageMatrices(avgOrbitSizes.values, numReplications)
    }

    return result
}

/**
 * Average a collection of matrices element-wise.
 */
private fun averageMatrices(matrices: Collection<Matrix>, count: Int): Matrix {
    if (matrices.isEmpty()) {
        throw IllegalArgumentException("Cannot average empty collection")
    }
    val first = matrices.first()
    val result = Matrix(first.numRows, first.numCols)
    result.zero()
    for (m in matrices) {
        result.addEq(m)
    }
    return Matrix.scaleMult(result, 1.0 / count.toDouble())
}

/**
 * Compute cross-replication confidence interval half-widths.
 *
 * Uses the formula: CI = t_{α/2, n-1} × √(s²/n)
 * where s² is the sample variance across replications.
 */
private fun computeCrossReplicationCI(
    replicationResults: List<Matrix>,
    grandMean: Matrix,
    tCrit: Double,
    numReplications: Int
): Matrix {
    val rows = grandMean.numRows
    val cols = grandMean.numCols
    val ciHalfWidth = Matrix(rows, cols)

    if (numReplications < 2) {
        ciHalfWidth.fill(0.0)
        return ciHalfWidth
    }

    for (i in 0 until rows) {
        for (j in 0 until cols) {
            val mean = grandMean.get(i, j)
            var sumSquaredDiff = 0.0
            for (rep in replicationResults) {
                val diff = rep.get(i, j) - mean
                sumSquaredDiff += diff * diff
            }
            val variance = sumSquaredDiff / (numReplications - 1)
            val stdError = kotlin.math.sqrt(variance / numReplications)
            ciHalfWidth.set(i, j, tCrit * stdError)
        }
    }

    return ciHalfWidth
}

/**
 * Compute relative precision (CI half-width / mean).
 */
private fun computeRelativePrecision(ciHalfWidth: Matrix, mean: Matrix): Matrix {
    val result = Matrix(mean.numRows, mean.numCols)
    for (i in 0 until mean.numRows) {
        for (j in 0 until mean.numCols) {
            val m = mean.get(i, j)
            val ci = ciHalfWidth.get(i, j)
            if (m > 0) {
                result.set(i, j, ci / m)
            } else {
                result.set(i, j, 0.0)
            }
        }
    }
    return result
}

/**
 * Check if all metrics have converged (relative precision < tolerance).
 */
private fun checkConvergence(
    qnRelPrec: Matrix,
    unRelPrec: Matrix,
    rnRelPrec: Matrix,
    tnRelPrec: Matrix,
    tolerance: Double
): Boolean {
    val allMatrices = listOf(qnRelPrec, unRelPrec, rnRelPrec, tnRelPrec)
    for (m in allMatrices) {
        for (i in 0 until m.numRows) {
            for (j in 0 until m.numCols) {
                val value = m.get(i, j)
                // Skip zero/NaN values (no data) and check if any exceed tolerance
                if (value > 0 && value.isFinite() && value > tolerance) {
                    return false
                }
            }
        }
    }
    return true
}
