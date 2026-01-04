package jline.lib.fjcodes

import jline.util.matrix.Matrix

/**
 * Build PH representation for 2-node Fork-Join job
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Build service representation for a 2-node FJ job
 *
 * This function constructs the Phase-Type (PH) representation for the service
 * time of a 2-node Fork-Join job, which is the maximum of two independent
 * service times.
 *
 * @param service Single subtask service process
 * @return FJServiceH with PH representation for 2-node job
 */
fun build_Service_h(service: FJService): FJServiceH {
    val dim_single = service.tau_st.length()
    val phases_single = build_index(dim_single, 1)

    // Possible service phases for the 2-node FJ job
    // Each row represents a state: [longest queue phases, shortest queue phases]
    val totalStates = dim_single * dim_single
    val service_phases = Matrix(totalStates, 2 * dim_single)

    var k = 0
    for (i in 0 until dim_single) {
        for (j in 0 until dim_single) {
            // [longest, shortest]
            for (col in 0 until dim_single) {
                service_phases.set(k, col, phases_single.get(i, col))
                service_phases.set(k, dim_single + col, phases_single.get(j, col))
            }
            k++
        }
    }

    // PH representation for the service time of a 2-node FJ job
    // beta = kron(tau_st, tau_st)
    val beta = service.tau_st.kron(service.tau_st).transpose()  // Row vector

    // S = kronsum(ST, ST) = ST \otimes I + I \otimes ST
    val S = kronsum(service.ST, service.ST)

    return FJServiceH(service_phases, beta, S)
}
