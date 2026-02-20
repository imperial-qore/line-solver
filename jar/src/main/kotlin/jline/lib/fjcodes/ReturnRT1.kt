package jline.lib.fjcodes

import jline.io.line_warning
import jline.lib.butools.MMAPPH1FCFS
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Compute response time percentiles for K=1 (single queue case)
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Compute response time percentiles for K=1 Fork-Join queue
 *
 * For the K=1 case, this reduces to a standard MAP/PH/1 queue.
 * Uses MMAPPH1FCFS to solve the queueing system.
 *
 * @param arrival Arrival process (lambda0, lambda1)
 * @param service Service process
 * @param pers Array of percentile levels to compute
 * @return Matrix with percentile results [percentile_level, percentile_value]
 */
fun returnRT1(arrival: FJArrival, service: FJService, pers: DoubleArray): Matrix {
    // For K=1, the Fork-Join queue reduces to a single MAP/PH/1 queue

    // Use MMAPPH1FCFS to solve the MAP/PH/1 queue
    // Prepare input for MMAPPH1FCFS
    val D = MatrixCell(2)
    D[0] = arrival.lambda0
    D[1] = arrival.lambda1

    val sigma = HashMap<Int?, Matrix>()
    val S = HashMap<Int?, Matrix>()

    // Service process: exit vector times initial probability
    // service.St is exit rate vector (column), service.tau_st is initial probability (column)
    // For PH representation: sigma is initial probability (row), S is sub-generator
    sigma[0] = service.tau_st.transpose()  // Convert column to row
    S[0] = service.ST

    try {
        // Call MMAPPH1FCFS to get sojourn time distribution
        val result = MMAPPH1FCFS(
            D = D,
            sigma = sigma,
            S = S,
            numOfQLMoms = null,
            numOfQLProbs = null,
            numOfSTMoms = null,
            stDistr = null,
            stDistrME = false,
            stDistrPH = true,
            prec = 1e-14,
            classes_ = null
        )

        // Extract sojourn time (response time) PH distribution from result
        // MMAPPH1FCFS returns stDistrPH_alpha and stDistrPH_A for the sojourn time
        val stDistrPH_alpha = result["stDistrPH_alpha"]
        val stDistrPH_A = result["stDistrPH_A"]

        if (stDistrPH_alpha == null || stDistrPH_A == null) {
            throw RuntimeException("MMAPPH1FCFS did not return sojourn time PH distribution")
        }

        // Get the distribution for class 0 (single class case)
        val res_alpha: Matrix? = stDistrPH_alpha[0]
        val Smat: Matrix? = stDistrPH_A[0]

        if (res_alpha == null || Smat == null) {
            throw RuntimeException("MMAPPH1FCFS returned null sojourn time distribution for class 0")
        }



        // Compute percentiles from PH sojourn time distribution
        return returnPer(res_alpha, Smat, pers)

    } catch (e: Exception) {
        line_warning("returnRT1", "MMAPPH1FCFS failed, using simplified calculation: %s", e.message)

        // Fallback: use service distribution only (ignoring queueing)
        val res_alpha = service.tau_st.transpose()
        val Smat = service.ST
        return returnPer(res_alpha, Smat, pers)
    }
}
