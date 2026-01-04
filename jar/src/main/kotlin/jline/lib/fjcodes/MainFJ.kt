package jline.lib.fjcodes

import jline.util.matrix.Matrix
import kotlin.math.ln

/**
 * Main orchestrator for FJ_codes Fork-Join percentile analysis
 *
 * Implements the approximation method from "Beyond the Mean in Fork-Join Queues:
 * Efficient Approximation for Response-Time Tails" (IFIP Performance 2015).
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Result for a specific K value
 */
data class FJPercentileResult(
    val K: Int,                     // Number of parallel queues
    val percentiles: DoubleArray,   // Percentile levels (0-100 scale)
    val RTp: DoubleArray            // Response time percentile values
)

/**
 * Compute response time percentiles for K-node Fork-Join queue
 *
 * For K=1: Uses exact MAP/PH/1 analysis
 * For K=2: Uses approximation from Section 4 of the paper
 * For K>2: Uses logarithmic extrapolation from K=1 and K=2 results
 *
 * @param arrival Arrival process (lambda, lambda0, lambda1)
 * @param service Service process (mu, ST, St, tau_st)
 * @param pers Percentile levels to compute (e.g., [0.50, 0.90, 0.95, 0.99])
 * @param K Array of K values (number of parallel queues) to analyze
 * @param Cs Array of C values (accuracy parameter, default [100])
 * @param tMode T-matrix computation method ("NARE" or "Sylvest")
 * @return List of FJPercentileResult, one for each K value
 */
fun mainFJ(
    arrival: FJArrival,
    service: FJService,
    pers: DoubleArray,
    K: IntArray,
    Cs: IntArray = intArrayOf(100),
    tMode: String = "NARE"
): List<FJPercentileResult> {
    // Check stability condition
    val load = arrival.lambda / service.mu
    if (load >= 1.0) {
        throw IllegalArgumentException(
            "System not stable: mean arrival rate ${arrival.lambda} >= " +
            "mean service rate ${service.mu} (load = $load)"
        )
    }

    // Use the last C value from Cs array
    val C = Cs.last()

    // Compute percentiles for K=1 (exact)
    val percentileRT_1 = returnRT1(arrival, service, pers)

    // Compute percentiles for K=2 (approximation)
    val percentileRT_2 = returnRT2(arrival, service, pers, C, tMode)

    // Predict percentiles for requested K values using logarithmic extrapolation
    val results = mutableListOf<FJPercentileResult>()

    for (k in K) {
        val percentilesScaled = DoubleArray(pers.size) { pers[it] * 100.0 }
        val RTp = DoubleArray(pers.size)

        for (p in pers.indices) {
            if (k == 1) {
                // Exact result for K=1
                RTp[p] = percentileRT_1.get(p, 1)
            } else if (k == 2) {
                // Exact result for K=2
                RTp[p] = percentileRT_2.get(p, 1)
            } else {
                // Logarithmic extrapolation for K > 2
                // RT_K = RT_1 + (RT_2 - RT_1) * log(K) / log(2)
                val rt1 = percentileRT_1.get(p, 1)
                val rt2 = percentileRT_2.get(p, 1)
                RTp[p] = rt1 + (rt2 - rt1) * ln(k.toDouble()) / ln(2.0)
            }
        }

        results.add(FJPercentileResult(k, percentilesScaled, RTp))
    }

    return results
}

/**
 * Convenience overload for single K value
 */
fun mainFJ(
    arrival: FJArrival,
    service: FJService,
    pers: DoubleArray,
    K: Int,
    C: Int = 100,
    tMode: String = "NARE"
): FJPercentileResult {
    return mainFJ(arrival, service, pers, intArrayOf(K), intArrayOf(C), tMode).first()
}
