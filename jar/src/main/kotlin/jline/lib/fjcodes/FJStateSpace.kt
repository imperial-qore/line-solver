package jline.lib.fjcodes

import jline.util.matrix.Matrix

/**
 * FJ_codes state space construction functions
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Service representation for FJ_codes
 */
data class FJService(
    val mu: Double,           // Service rate
    val ST: Matrix,           // PH sub-generator
    val St: Matrix,           // PH exit rate vector (column)
    val tau_st: Matrix,       // PH initial probability (column)
    val SerChoice: Int        // Service type: 1=Exponential, other=PH
)

/**
 * Service representation for 2-node FJ job
 */
data class FJServiceH(
    val service_phases: Matrix,  // Possible phase combinations (rows are states)
    val beta: Matrix,             // Initial probability for 2-node job (row)
    val S: Matrix                 // PH sub-generator for 2-node job
)

/**
 * Construct state space matrix for not-all-busy states
 *
 * @param C Capacity parameter
 * @param services Service process
 * @param service_h 2-node service representation
 * @return State transition matrix for not-all-busy states
 */
fun constructNotAllBusy(C: Int, services: FJService, service_h: FJServiceH): Matrix {
    val dim = service_h.beta.length()
    val m = services.tau_st.length()

    val indexes_notbusy = build_index(m, 1)
    val dim_NB = indexes_notbusy.getNumRows()
    val dim_C = C + 1
    val dim_notbusy = (dim_C - 1) * dim_NB + 1

    val S_notallbusy = Matrix(dim_notbusy, dim_notbusy)

    // From not-busy to not-busy
    for (row in 0 until dim_C - 1) {
        val startRow = row * dim_NB
        val endRow = (row + 1) * dim_NB
        setSubMatrix(S_notallbusy, startRow, startRow, services.ST)
    }

    // A = -sum(ST, 2) * tau_st  (row sum times initial prob)
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

    // Transitions between levels
    for (row in 0 until dim_C - 2) {
        val startRow = row * dim_NB
        val endRow = (row + 1) * dim_NB
        val startCol = (row + 1) * dim_NB
        setSubMatrix(S_notallbusy, startRow, startCol, A)
    }

    // Last transition to absorbing state
    val startRow = (dim_C - 2) * dim_NB
    val endRow = (dim_C - 1) * dim_NB
    val startCol = (dim_C - 1) * dim_NB
    for (i in startRow until endRow) {
        var rowSum = 0.0
        for (j in 0 until m) {
            rowSum += services.ST.get(i - startRow, j)
        }
        S_notallbusy.set(i, startCol, -rowSum)
    }

    // Normalize diagonal
    for (row in 0 until dim_notbusy) {
        var rowSum = 0.0
        for (col in 0 until dim_notbusy) {
            if (row != col) {
                rowSum += S_notallbusy.get(row, col)
            }
        }
        S_notallbusy.set(row, row, -rowSum)
    }

    return S_notallbusy
}

/**
 * Result of constructSRK
 */
data class SRKResult(
    val Se: Matrix,
    val Sestar: Matrix,
    val R0: Matrix,
    val Ke: Matrix,
    val Kc: Matrix,
    val S_long: Matrix,
    val S_last: Matrix
)

/**
 * Construct extended state space matrices (Se, Sestar, R0, Ke, Kc)
 *
 * @param C Capacity parameter
 * @param services Service process
 * @param service_h 2-node service representation
 * @param S Busy state transition matrix
 * @return SRKResult with all constructed matrices
 */
fun constructSRK(C: Int, services: FJService, service_h: FJServiceH, S: Matrix): SRKResult {
    val dim = service_h.beta.length()
    val m = services.tau_st.length()

    val indexes_notbusy = build_index(m, 1)
    val dim_NB = indexes_notbusy.getNumRows()
    val dim_C = C + 1
    val newdim = dim_C * dim
    val dim_notbusy = (dim_C - 1) * dim_NB + 1

    var Se = Matrix(newdim + dim_notbusy, newdim + dim_notbusy)
    var Sestar = Matrix(newdim + dim_notbusy, newdim + dim_notbusy)

    // Copy S into upper-left block
    setSubMatrix(Se, 0, 0, S)

    // From busy to not-busy: when the job in the shorter queue completes service
    val S_long = Matrix(dim, dim_NB)

    for (row in 0 until dim) {
        val countvect = getRowAsArray(service_h.service_phases, row)

        for (i in m until 2 * m) {  // Starting phase
            if (countvect[i] > 0) {
                // Extract tovect (first m elements)
                val tovect = DoubleArray(m)
                for (k in 0 until m) {
                    tovect[k] = countvect[k]
                }

                val col = vectmatch(tovect, indexes_notbusy) - 1  // Convert to 0-based

                if (col >= 0) {
                    // The job in the shorter queue completes service
                    val currentVal = S_long.get(row, col)
                    S_long.set(row, col, currentVal + countvect[i] * services.St.get(i - m, 0))
                }
            }
        }
    }

    // Set S_long blocks in Se
    setSubMatrix(Se, 0, newdim, S_long)
    for (row in 1 until dim_C - 1) {
        setSubMatrix(Se, row * dim, newdim + (row - 1) * dim_NB, S_long)
    }

    // When the 2 queues have the same length
    val S_last = Matrix(dim, dim_NB)

    for (row in 0 until dim) {
        val countvect = getRowAsArray(service_h.service_phases, row)

        for (k in 0..1) {
            val startIdx = k * m
            val endIdx = (k + 1) * m

            for (i in startIdx until endIdx) {
                if (countvect[i] > 0) {
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

    setSubMatrix(Se, (dim_C - 1) * dim, newdim + (dim_C - 2) * dim_NB, S_last)

    // Sestar gets transitions from busy to not-busy
    setSubMatrix(Sestar, 0, newdim, Matrix.getSubMatrix(Se, 0, newdim, newdim, newdim + dim_notbusy))

    // From not-busy to not-busy (same as constructNotAllBusy logic)
    for (row in 0 until dim_C - 1) {
        val start = newdim + row * dim_NB
        setSubMatrix(Se, start, start, services.ST)
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

    for (row in 0 until dim_C - 2) {
        val startRow = newdim + row * dim_NB
        val startCol = newdim + (row + 1) * dim_NB
        setSubMatrix(Se, startRow, startCol, A)
    }

    // Last block
    val lastBlockRow = newdim + (dim_C - 2) * dim_NB
    val lastBlockCol = newdim + (dim_C - 1) * dim_NB
    for (i in 0 until dim_NB) {
        var rowSum = 0.0
        for (j in 0 until m) {
            rowSum += services.ST.get(i, j)
        }
        Se.set(lastBlockRow + i, lastBlockCol, -rowSum)
    }

    // Normalize diagonal for not-busy part
    for (row in newdim until newdim + dim_notbusy) {
        var rowSum = 0.0
        for (col in 0 until newdim + dim_notbusy) {
            if (row != col) {
                rowSum += Se.get(row, col)
            }
        }
        Se.set(row, row, -rowSum)
    }

    // Construct R0 matrix: from not-busy to busy
    val R0 = Matrix(newdim + dim_notbusy, newdim + dim_notbusy)
    val R_NB = Matrix(dim_NB, dim)

    for (row in 0 until dim_NB) {
        val countvect = getRowAsArray(indexes_notbusy, row)

        for (i in 0 until m) {
            // tovect_short has a 1 at position i
            val tovect = DoubleArray(2 * m)
            for (j in 0 until m) {
                tovect[j] = countvect[j]
            }
            tovect[m + i] = 1.0

            val col = vectmatch(tovect, service_h.service_phases) - 1  // Convert to 0-based

            if (col >= 0) {
                val currentVal = R_NB.get(row, col)
                R_NB.set(row, col, currentVal + services.tau_st.get(i, 0))
            }
        }
    }

    for (row in 0 until dim_C - 1) {
        setSubMatrix(R0, newdim + row * dim_NB, row * dim, R_NB)
    }

    // Last row of R0
    for (i in 0 until dim) {
        R0.set(newdim + dim_notbusy - 1, newdim - dim + i, service_h.beta.get(0, i))
    }

    // Construct Ke and Kc
    val Ke = Matrix(newdim, newdim + dim_notbusy)
    for (i in 0 until newdim) {
        Ke.set(i, i, 1.0)
    }

    val Kc = Matrix(newdim + dim_notbusy, newdim)
    for (i in 0 until newdim) {
        Kc.set(i, i, 1.0)
    }

    return SRKResult(Se, Sestar, R0, Ke, Kc, S_long, S_last)
}
