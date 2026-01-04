package jline.lib.fjcodes

import jline.util.matrix.Matrix

/**
 * Compute T-matrix using either Sylvester or NARE method
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Result of computeT
 */
data class ComputeTResult(
    val T: Matrix,
    val S: Matrix,
    val A_jump: Matrix,
    val S_Arr: Matrix,
    val sum_Ajump: Matrix
)

/**
 * Compute T-matrix using specified method
 *
 * The T-matrix is computed using either the Sylvester equation approach
 * or the NARE (Nonsymmetric Algebraic Riccati Equation) method.
 *
 * @param arrival Arrival process (lambda0, lambda1)
 * @param services Service process for single subtask
 * @param service_h Service representation for 2-node job
 * @param C Capacity parameter
 * @param tMode Method to use: "NARE" (default) or "Sylvest"
 * @return ComputeTResult with T, S, A_jump, S_Arr, sum_Ajump
 */
fun computeT(
    arrival: FJArrival,
    services: FJService,
    service_h: FJServiceH,
    C: Int,
    tMode: String = "NARE"
): ComputeTResult {
    // Build S and A_jump matrices
    val saResult = build_SA(services, service_h, C)
    val S = saResult.S
    val A_jump = saResult.A_jump

    val d0 = arrival.lambda0.getNumRows()
    val S_Arr = S.kron(Matrix.eye(d0))
    val A_jump_Arr = A_jump.kron(Matrix.eye(d0))

    // Compute T-matrix using specified method
    val T = if (tMode.contains("Sylvest", ignoreCase = true)) {
        // Sylvester equation approach
        computeT_Sylvester(arrival.lambda0, arrival.lambda1, S_Arr, A_jump)
    } else {
        // NARE method (default)
        computeT_NARE(arrival.lambda0, arrival.lambda1, S_Arr, A_jump)
    }

    // Compute sum of A_jump_Arr rows
    val sum_Ajump = Matrix(A_jump_Arr.getNumRows(), 1)
    for (i in 0 until A_jump_Arr.getNumRows()) {
        var sum = 0.0
        for (j in 0 until A_jump_Arr.getNumCols()) {
            sum += A_jump_Arr.get(i, j)
        }
        sum_Ajump.set(i, 0, sum)
    }

    return ComputeTResult(T, S, A_jump, S_Arr, sum_Ajump)
}
