package jline.lib.fjcodes

import jline.util.matrix.Matrix

/**
 * Generate service time Phase-Type representation for Fork-Join queue
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Result of generateService
 */
data class GenerateServiceResult(
    val T: Matrix,          // Service time PH sub-generator
    val newdim: Int,        // Dimension of busy states
    val dim_notbusy: Int    // Dimension of not-all-busy states
)

/**
 * Generate Phase-Type representation for service time
 *
 * Constructs the PH representation for the service time distribution
 * considering both busy and not-all-busy states.
 *
 * @param services Service process for single subtask
 * @param service_h Service representation for 2-node job
 * @param C Capacity parameter
 * @param S State transition matrix
 * @return GenerateServiceResult with service time PH representation
 */
fun generateService(
    services: FJService,
    service_h: FJServiceH,
    C: Int,
    S: Matrix
): GenerateServiceResult {
    val dim = service_h.beta.length()
    val m = services.tau_st.length()

    val indexes_notbusy = build_index(m, 1)
    val dim_NB = indexes_notbusy.getNumRows()

    val dim_C = C + 1

    val newdim = dim_C * dim
    val dim_notbusy = dim_C * dim_NB

    val T = Matrix(newdim + dim_notbusy, newdim + dim_notbusy)
    val t = Matrix(newdim + dim_notbusy, 1)

    // Set exit rates for not-busy states
    for (i in 0 until dim_NB) {
        t.set(newdim + dim_notbusy - dim_NB + i, 0, services.St.get(i, 0))
    }

    // Copy S into upper-left block
    setSubMatrix(T, 0, 0, S)

    // S_long: job in shorter queue completes service
    val S_long = Matrix(dim, dim_NB)

    for (row in 0 until dim) {
        val countvect = getRowAsArray(service_h.service_phases, row)

        for (i in m until 2 * m) {  // Starting phase (shorter queue)
            if (countvect[i] > 0.0) {
                // Extract tovect (first m elements)
                val tovect = DoubleArray(m)
                for (k in 0 until m) {
                    tovect[k] = countvect[k]
                }

                val col = vectmatch(tovect, indexes_notbusy) - 1  // Convert to 0-based

                if (col >= 0) {
                    val currentVal = S_long.get(row, col)
                    S_long.set(row, col, currentVal + countvect[i] * services.St.get(i - m, 0))
                }
            }
        }
    }

    // Set S_long blocks
    for (row in 0 until dim_C - 1) {
        setSubMatrix(T, row * dim, newdim + row * dim_NB, S_long)
    }

    // S_last: all queues have the same length
    val S_last = Matrix(dim, dim_NB)

    for (row in 0 until dim) {
        val countvect = getRowAsArray(service_h.service_phases, row)

        for (k in 0..1) {
            val startIdx = k * m
            val endIdx = (k + 1) * m

            for (i in startIdx until endIdx) {
                if (countvect[i] > 0.0) {
                    // Extract tovect from the other queue
                    val otherStart = (1 - k) * m
                    val tovect = DoubleArray(m)
                    for (j in 0 until m) {
                        tovect[j] = countvect[otherStart + j]
                    }

                    val col = vectmatch(tovect, indexes_notbusy) - 1  // Convert to 0-based

                    if (col >= 0) {
                        val currentVal = S_last.get(row, col)
                        S_last.set(row, col, currentVal + countvect[i] * services.St.get(i - startIdx, 0))
                    }
                }
            }
        }
    }

    setSubMatrix(T, (dim_C - 1) * dim, newdim + (dim_C - 1) * dim_NB, S_last)

    // Only one subtask in service
    for (row in 0 until dim_C) {
        setSubMatrix(T, newdim + row * dim_NB, newdim + row * dim_NB, services.ST)
    }

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

    // Set A blocks
    for (row in 0 until dim_C - 1) {
        setSubMatrix(T, newdim + row * dim_NB, newdim + (row + 1) * dim_NB, A)
    }

    // Normalize diagonal for not-busy states
    for (row in newdim until newdim + dim_notbusy) {
        T.set(row, row, 0.0)
        var rowSum = 0.0
        for (col in 0 until newdim + dim_notbusy) {
            rowSum += T.get(row, col)
        }
        rowSum += t.get(row, 0)
        T.set(row, row, -rowSum)
    }

    return GenerateServiceResult(T, newdim, dim_notbusy)
}
