package jline.solvers.mam

import jline.api.fj.extractFJParams
import jline.api.fj.isFJ
import jline.lang.NetworkStruct
import jline.lib.fjcodes.FJPercentileResult
import jline.lib.fjcodes.mainFJ
import jline.util.matrix.Matrix

/**
 * Fork-Join solver integration for SolverMAM
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Result from Fork-Join analysis
 */
data class MAMFJResult(
    val QN: Matrix,                             // Queue length
    val UN: Matrix,                             // Utilization
    val RN: Matrix,                             // Response time
    val TN: Matrix,                             // Throughput
    val XN: Matrix,                             // System throughput
    val percentileResults: List<FJPercentileResult>  // Percentile data per class
)

/**
 * Solve Fork-Join network using FJ_codes algorithm
 *
 * Automatically detects FJ topology and computes response time percentiles
 * using the algorithm from "Beyond the Mean in Fork-Join Queues" (IFIP Performance 2015).
 *
 * @param sn Network structure
 * @param percentiles Percentile levels to compute (0-1 scale, e.g., [0.50, 0.90, 0.95, 0.99])
 * @param C Accuracy parameter (default: 100)
 * @param tMode T-matrix method: "NARE" or "Sylvest" (default: "NARE")
 * @return MAMFJResult with performance metrics and percentiles
 * @throws IllegalArgumentException if network is not valid FJ topology
 */
fun solverMAMFJ(
    sn: NetworkStruct,
    percentiles: DoubleArray = doubleArrayOf(0.50, 0.90, 0.95, 0.99),
    C: Int = 100,
    tMode: String = "NARE"
): MAMFJResult {
    // Validate FJ topology
    val (isValid, fjInfo) = isFJ(sn)
    if (!isValid || fjInfo == null) {
        throw IllegalArgumentException(
            "Network is not a valid Fork-Join topology. " +
            "Required: Source → Fork → K Queues → Join → Sink"
        )
    }

    // Extract FJ parameters per class
    val (arrivals, services) = extractFJParams(sn, fjInfo)

    // Compute percentiles for each class
    val percentileResults = mutableListOf<FJPercentileResult>()

    for (r in 0 until sn.nclasses) {
        try {
            val result = mainFJ(
                arrival = arrivals[r],
                service = services[r],
                pers = percentiles,
                K = fjInfo.K,
                C = C,
                tMode = tMode
            )
            percentileResults.add(result)
        } catch (e: Exception) {
            throw RuntimeException(
                "FJ_codes failed for class $r: ${e.message}", e
            )
        }
    }

    // Compute mean response times from percentiles (via trapezoidal integration)
    val meanRT = DoubleArray(sn.nclasses)
    for (r in 0 until sn.nclasses) {
        val result = percentileResults[r]
        // Simple trapezoidal integration over percentile CDF
        var sum = 0.0
        for (i in 0 until result.RTp.size - 1) {
            val dp = result.percentiles[i + 1] - result.percentiles[i]
            sum += (result.RTp[i] + result.RTp[i + 1]) / 2.0 * dp / 100.0
        }
        meanRT[r] = sum
    }

    // Compute performance metrics using Little's Law
    val QN = Matrix(sn.nstations, sn.nclasses)
    val RN = Matrix(sn.nstations, sn.nclasses)
    val TN = Matrix(sn.nstations, sn.nclasses)
    val UN = Matrix(sn.nstations, sn.nclasses)
    val XN = Matrix(1, sn.nclasses)

    for (r in 0 until sn.nclasses) {
        val lambda = arrivals[r].lambda
        val mu = services[r].mu

        // System throughput = arrival rate (stable system)
        XN.set(0, r, lambda)

        // Distribute metrics across parallel queues
        for (qIdx in fjInfo.queueIndices) {
            val station = sn.nodeToStation.get(qIdx).toInt()

            // Mean response time per queue (assumes equal split)
            RN.set(station, r, meanRT[r] / fjInfo.K)

            // Throughput per queue
            TN.set(station, r, lambda / fjInfo.K)

            // Queue length via Little's Law: Q = λ * R
            QN.set(station, r, lambda * meanRT[r] / fjInfo.K)

            // Utilization: U = λ / (K * μ)
            UN.set(station, r, lambda / (fjInfo.K * mu))
        }

        // Fork, Join, Source, Sink get zero metrics (instantaneous)
        // Already initialized to zero by Matrix constructor
    }

    return MAMFJResult(QN, UN, RN, TN, XN, percentileResults)
}

/**
 * Interpolate percentile from stored results
 *
 * @param storedPercentiles Stored percentile levels (0-100 scale)
 * @param storedValues Stored percentile values
 * @param targetPercentile Target percentile to interpolate (0-100 scale)
 * @return Interpolated value
 */
fun interpolatePercentile(
    storedPercentiles: DoubleArray,
    storedValues: DoubleArray,
    targetPercentile: Double
): Double {
    // Find bracketing indices
    for (i in 0 until storedPercentiles.size - 1) {
        if (targetPercentile >= storedPercentiles[i] &&
            targetPercentile <= storedPercentiles[i + 1]) {
            // Linear interpolation
            val t = (targetPercentile - storedPercentiles[i]) /
                    (storedPercentiles[i + 1] - storedPercentiles[i])
            return storedValues[i] + t * (storedValues[i + 1] - storedValues[i])
        }
    }

    // Extrapolation (outside range)
    if (targetPercentile < storedPercentiles[0]) {
        return storedValues[0]
    } else {
        return storedValues[storedValues.size - 1]
    }
}
