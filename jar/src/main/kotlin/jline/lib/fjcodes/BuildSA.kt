package jline.lib.fjcodes

import jline.util.matrix.Matrix

/**
 * Build state transition matrices for FJ queue
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Result of build_SA
 */
data class SAResult(
    val S: Matrix,      // State transition matrix (phase changes)
    val A_jump: Matrix  // Jump matrix (completions)
)

/**
 * Build state transition matrices S and A_jump
 *
 * Constructs the transition rate matrices for the FJ queue:
 * - S: Phase changes within each level (no job completions)
 * - A_jump: Transitions due to job completions
 *
 * @param services Service process for single subtask
 * @param service_h Service representation for 2-node job
 * @param C Capacity parameter
 * @return SAResult with S and A_jump matrices
 */
fun build_SA(services: FJService, service_h: FJServiceH, C: Int): SAResult {
    val dim = service_h.beta.length()
    val m = services.tau_st.length()
    val dim_C = C + 1
    val newdim = dim_C * dim

    val S = Matrix(newdim, newdim)
    val A_jump = Matrix(newdim, newdim)

    // Phase changes within each level c in C (service_h.S)
    // No job completes service
    for (row in 0 until dim_C) {
        setSubMatrix(S, row * dim, row * dim, service_h.S)
    }

    // When a job completes service and a new job starts service
    // A = -sum(ST, 2) * tau_st
    val A = Matrix(m, m)
    for (i in 0 until m) {
        var rowSum = 0.0
        for (j in 0 until m) {
            rowSum += services.ST.get(i, j)
        }
        for (j in 0 until m) {
            A.set(i, j, -rowSum * services.tau_st.get(j, 0))
        }
    }

    // When the job in the longer queue completes service
    val S_Cminus1 = Matrix(dim, dim)

    for (row in 0 until dim) {
        val countvect = getRowAsArray(service_h.service_phases, row)

        for (i in 0 until m) {  // The starting phase
            if (countvect[i] > 0.0) {
                for (j in 0 until m) {  // The resulting phase
                    val tovect = countvect.copyOf()
                    tovect[i] = tovect[i] - 1.0  // Jump from phase i to phase j
                    tovect[j] = tovect[j] + 1.0

                    val col = vectmatch(tovect, service_h.service_phases) - 1  // Convert to 0-based

                    if (col >= 0) {
                        // The job in the longer queue completes service
                        val currentVal = S_Cminus1.get(row, col)
                        S_Cminus1.set(row, col, currentVal + countvect[i] * A.get(i, j))
                    }
                }
            }
        }
    }

    // Form c to c-1 (queue length difference decreases)
    for (c in C downTo 1) {
        val sourceRow = (C - c) * dim
        val targetCol = (C - c + 1) * dim
        setSubMatrix(S, sourceRow, targetCol, S_Cminus1)
    }

    // When the job in the shorter queue completes service
    val A_Cplus1 = Matrix(dim, dim)

    for (row in 0 until dim) {
        val countvect = getRowAsArray(service_h.service_phases, row)

        for (i in m until 2 * m) {  // The starting phase (shorter queue)
            if (countvect[i] > 0.0) {
                for (j in m until 2 * m) {  // The resulting phase
                    val tovect = countvect.copyOf()
                    tovect[i] = tovect[i] - 1.0  // Jump from phase i to phase j
                    tovect[j] = tovect[j] + 1.0

                    val col = vectmatch(tovect, service_h.service_phases) - 1  // Convert to 0-based

                    if (col >= 0) {
                        // The job in the shorter queue completes service
                        val currentVal = A_Cplus1.get(row, col)
                        A_Cplus1.set(row, col, currentVal + A.get(i - m, j - m))
                    }
                }
            }
        }
    }

    // Set A_jump blocks (queue length difference increases)
    for (c in (C - 1) downTo 1) {
        val sourceRow = (C - c) * dim
        val targetCol = (C - c - 1) * dim
        setSubMatrix(A_jump, sourceRow, targetCol, A_Cplus1)
    }

    setSubMatrix(A_jump, 0, 0, A_Cplus1)

    // Both queues have the same length
    // The one that didn't finish becomes the longer queue
    val A_last = Matrix(dim, dim)

    for (row in 0 until dim) {
        val countvect = getRowAsArray(service_h.service_phases, row)

        for (k in 0..1) {
            val startIdx = k * m
            val endIdx = (k + 1) * m

            for (i in startIdx until endIdx) {  // The starting phase
                if (countvect[i] > 0.0) {
                    // Extract tovect from the other queue
                    val otherStart = (1 - k) * m
                    val tovect = DoubleArray(2 * m)

                    // Copy the phases of the queue that didn't complete
                    for (idx in 0 until m) {
                        tovect[idx] = countvect[otherStart + idx]
                    }

                    // The new job starts in phase j
                    for (j in 0 until m) {
                        val temp_tovector = tovect.copyOf()
                        temp_tovector[m + j] = 1.0

                        val col = vectmatch(temp_tovector, service_h.service_phases) - 1  // Convert to 0-based

                        if (col >= 0) {
                            val currentVal = A_last.get(row, col)
                            A_last.set(row, col, currentVal + countvect[i] * A.get(i - startIdx, j))
                        }
                    }
                }
            }
        }
    }

    setSubMatrix(A_jump, C * dim, (C - 1) * dim, A_last)

    return SAResult(S, A_jump)
}
