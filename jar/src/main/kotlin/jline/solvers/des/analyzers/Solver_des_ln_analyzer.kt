/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.des.analyzers

import jline.lang.layered.LayeredNetworkStruct
import jline.solvers.SolverOptions
import jline.solvers.des.LNDESResult
import jline.solvers.des.SolverDES
import jline.solvers.des.handlers.solver_ssj_ln

/**
 * LN DES analyzer - validates LQN structure and dispatches to SSJ backend.
 *
 * @param lsn the LayeredNetwork structure
 * @param options solver options
 * @param solverDES the parent solver instance
 * @return LNDESResult containing simulation results
 */
fun solver_des_ln_analyzer(
    lsn: LayeredNetworkStruct,
    options: SolverOptions,
    solverDES: SolverDES
): LNDESResult {
    val Tstart = System.nanoTime()

    // Validate LQN structure
    validateLQNStructure(lsn)

    // Map method name
    var backendMethod = options.method
    if (backendMethod == null || backendMethod == "default") {
        backendMethod = "ssj"
    }

    // Route to appropriate backend
    val result = when (backendMethod) {
        "ssj" -> solver_ssj_ln(lsn, options)
        else -> throw RuntimeException("solver_des_ln_analyzer: Unknown method '$backendMethod'")
    }

    // Preserve original method name and calculate runtime
    result.method = options.method ?: "default"
    result.runtime = (System.nanoTime() - Tstart).toDouble() / 1e9

    return result
}

/**
 * Validates the LayeredNetwork structure for DES simulation.
 *
 * @param lsn the LayeredNetwork structure to validate
 * @throws RuntimeException if validation fails
 */
private fun validateLQNStructure(lsn: LayeredNetworkStruct) {
    // Check basic structure requirements
    if (lsn.nhosts == 0) {
        throw RuntimeException("solver_des_ln_analyzer: LQN must have at least one Host/Processor")
    }
    if (lsn.ntasks == 0) {
        throw RuntimeException("solver_des_ln_analyzer: LQN must have at least one Task")
    }

    // Validate all tasks have at least one entry
    for (t in 1..lsn.ntasks) {
        val tidx = lsn.tshift + t
        val entries = lsn.entriesof[tidx]
        if (entries == null || entries.isEmpty()) {
            val taskName = lsn.names[tidx] ?: "Task$tidx"
            throw RuntimeException("solver_des_ln_analyzer: Task '$taskName' has no entries")
        }
    }

    // Validate all entries have a bound activity
    for (e in 1..lsn.nentries) {
        val eidx = lsn.eshift + e
        // Check that the entry has at least one activity bound to it
        val tidx = lsn.parent.get(0, eidx).toInt()
        val activities = lsn.actsof[tidx]
        if (activities == null || activities.isEmpty()) {
            val entryName = lsn.names[eidx] ?: "Entry$eidx"
            throw RuntimeException("solver_des_ln_analyzer: Entry '$entryName' has no activities in its task")
        }
    }
}
