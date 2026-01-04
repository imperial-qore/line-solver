/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam.handlers

import jline.io.line_warning
import jline.io.mfilename
import jline.lang.NetworkStruct
import jline.solvers.SolverOptions
import jline.solvers.mam.MAMOptions
import jline.solvers.mam.MAMResult
import jline.util.matrix.Matrix

/**
 * RCAT-based solver for SolverMAM.
 *
 * Uses RCAT (Reversed Compound Agent Theorem) to find product-form
 * solutions for queueing networks.
 *
 * Methods:
 *   'inap'     - Iterative Numerical Approximation Procedure (default, fast)
 *   'inapplus' - Improved INAP with weighted rates (no normalization)
 *   'exact'    - Optimization-based search using autocat (not yet available)
 */
fun solver_mam_ag(sn: NetworkStruct, options: SolverOptions): MAMResult {
    val M = sn.nstations
    val K = sn.nclasses

    val maxStates = if (options is MAMOptions) options.maxStates else 100

    val tol = if (options.iter_tol > 0) options.iter_tol else 1e-6
    val maxiter = if (options.iter_max > 0) options.iter_max else 1000

    val rcat = solver_mam_build_ag(sn, maxStates)

    val numProcesses = rcat.N.size
    val numActions = rcat.actionMap.size

    if (numProcesses == 0 || (numActions == 0 && numProcesses > 1)) {
        line_warning(mfilename(object {}), "Network could not be mapped to RCAT format.")
        return createEmptyMAMResult(M, K)
    }

    var method = options.method
    if (method == "default") {
        method = "inap"
    }

    val inapResult = when (method) {
        "inap" -> {
            solver_mam_inap(rcat, tol, maxiter, options.seed, "inap")
        }
        "inapplus" -> {
            solver_mam_inap(rcat, tol, maxiter, options.seed, "inapplus")
        }
        "exact" -> {
            line_warning(mfilename(object {}), "AutoCAT (exact method) is not yet available in JAR. Falling back to INAP.")
            solver_mam_inap(rcat, tol, maxiter, options.seed, "inap")
        }
        else -> {
            line_warning(mfilename(object {}), "Unknown method: $method. Using INAP.")
            solver_mam_inap(rcat, tol, maxiter, options.seed, "inap")
        }
    }

    val metrics = solver_mam_metrics(sn, inapResult, rcat)

    val result = MAMResult()
    result.QN = metrics.QN
    result.UN = metrics.UN
    result.RN = metrics.RN
    result.TN = metrics.TN
    result.CN = metrics.CN
    result.XN = metrics.XN
    result.iter = inapResult.iter
    result.method = method
    result.actionRates = inapResult.x
    result.equilibrium = inapResult.pi
    result.generators = inapResult.Q

    return result
}

private fun createEmptyMAMResult(M: Int, K: Int): MAMResult {
    val result = MAMResult()
    result.QN = Matrix(M, K)
    result.UN = Matrix(M, K)
    result.RN = Matrix(M, K)
    result.TN = Matrix(M, K)
    result.CN = Matrix(1, K)
    result.XN = Matrix(1, K)
    result.iter = 0
    result.method = "inap"
    return result
}
