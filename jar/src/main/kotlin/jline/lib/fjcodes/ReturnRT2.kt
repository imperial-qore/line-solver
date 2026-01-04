package jline.lib.fjcodes

import jline.util.matrix.Matrix

/**
 * Compute response time percentiles for K=2 (two parallel queues)
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Compute response time percentiles for K=2 Fork-Join queue
 *
 * Combines waiting time and service time distributions to compute
 * response time percentiles for a 2-node Fork-Join system.
 *
 * @param arrival Arrival process (lambda0, lambda1)
 * @param service Service process for single subtask
 * @param pers Array of percentile levels to compute
 * @param C Capacity parameter (accuracy control)
 * @param tMode T-matrix computation method ("NARE" or "Sylvest")
 * @return Matrix with percentile results [percentile_level, percentile_value]
 */
fun returnRT2(
    arrival: FJArrival,
    service: FJService,
    pers: DoubleArray,
    C: Int,
    tMode: String = "NARE"
): Matrix {
    // Build 2-node service representation
    val service_h = build_Service_h(service)

    // Compute T-matrix
    val computeTResult = computeT(arrival, service, service_h, C, tMode)
    val T = computeTResult.T
    val S = computeTResult.S
    val A_jump = computeTResult.A_jump
    val S_Arr = computeTResult.S_Arr
    val sum_Ajump = computeTResult.sum_Ajump

    val mWait = A_jump.getNumRows()

    // phi = sum(T - S_Arr, 2) (row sums)
    val phi = Matrix(T.getNumRows(), 1)
    for (i in 0 until T.getNumRows()) {
        var sum = 0.0
        for (j in 0 until T.getNumCols()) {
            sum += T.get(i, j) - S_Arr.get(i, j)
        }
        phi.set(i, 0, sum)
    }

    // Compute steady-state distribution
    val piResult = computePi(T, arrival, service, service_h, C, S, A_jump)
    var pi0 = piResult.pi0
    val En1 = piResult.En1

    // Compute waiting time distribution
    val waitResult = returnWait(En1, pi0, T, phi, sum_Ajump)
    val wait_alpha = waitResult.wait_alpha
    val wait_Smat = waitResult.wait_Smat
    val prob_wait = waitResult.prob_wait
    val alfa = waitResult.alfa


    // Generate service time representation
    val serviceResult = generateService(service, service_h, C, S)
    var ST = serviceResult.T
    val dim = serviceResult.newdim
    val dim_notbusy = serviceResult.dim_notbusy

    // Kronecker product with arrival process
    ST = ST.kron(Matrix.eye(arrival.lambda0.getNumRows()))

    // Normalize pi0
    pi0 = pi0.scale(1.0 / pi0.elementSum())

    val ma = arrival.lambda0.getNumRows()
    val dim_ma = ma * dim
    val dim_notbusy_ma = ma * dim_notbusy
    val dim_service = dim_ma + dim_notbusy_ma

    // Starting state of service for a job in not-all-busy period
    val notbusy_start = Matrix(1, dim_service)
    for (i in 0 until dim_ma) {
        notbusy_start.set(0, i, (1.0 - prob_wait) * pi0.get(0, i))
    }

    // Starting state of service for a job in all-busy
    val busy_start = Matrix(1, dim_service)

    // TS = T - S_Arr
    val TS = T.add(-1.0, S_Arr)

    // alfa * TS / sum(alfa * TS)
    val alfaTS = alfa.mult(TS)
    val sumAlfaTS = alfaTS.elementSum()

    for (i in 0 until dim_ma) {
        busy_start.set(0, i, prob_wait * alfaTS.get(0, i) / sumAlfaTS)
    }

    // PH representation of response time
    val Tr = ST.getNumRows()
    val Sc = wait_Smat.getNumCols()

    // TS_sum = sum(TS, 2) * ones(1, ma * mWait)
    val TS_sum = Matrix(TS.getNumRows(), ma * mWait)
    for (i in 0 until TS.getNumRows()) {
        var rowSum = 0.0
        for (j in 0 until TS.getNumCols()) {
            rowSum += TS.get(i, j)
        }
        for (j in 0 until ma * mWait) {
            TS_sum.set(i, j, rowSum)
        }
    }

    // TS = TS ./ TS_sum
    val TS_normalized = Matrix(TS.getNumRows(), TS.getNumCols())
    for (i in 0 until TS.getNumRows()) {
        for (j in 0 until TS.getNumCols()) {
            if (TS_sum.get(i, 0) > 1e-12) {
                TS_normalized.set(i, j, TS.get(i, j) / TS_sum.get(i, 0))
            }
        }
    }

    // stat_service_phase = busy_start / (-ST)
    val negST = ST.scale(-1.0)
    val stat_service_phase = busy_start.rightMatrixDivide(negST)

    // Find non-zero phases
    // Use a small negative threshold to match MATLAB's behavior with numerical precision
    // MATLAB stat_service_phase can have values as small as 1e-14 that are still counted
    val busy_nz = BooleanArray(Tr)
    var m_tr_ST = 0

    for (i in 0 until stat_service_phase.length()) {
        val val_i = stat_service_phase.get(0, i)
        // Allow tiny negative values from numerical error (MATLAB includes phases ~1e-14)
        if (val_i > -1e-15) {
            busy_nz[i] = true
            m_tr_ST++
        }
    }


    // tr_start_state = -sum(ST, 2)' * diag(stat_service_phase)
    val tr_start_state_full = Matrix(1, Tr)
    for (i in 0 until Tr) {
        var rowSum = 0.0
        for (j in 0 until ST.getNumCols()) {
            rowSum += ST.get(i, j)
        }
        tr_start_state_full.set(0, i, -rowSum * stat_service_phase.get(0, i))
    }

    // Build tr_ST (time-reversed service time matrix)
    val tr_ST_full = Matrix(Tr, Tr)
    for (i in 0 until Tr) {
        for (j in 0 until Tr) {
            if (stat_service_phase.get(0, i) > 1e-12) {
                tr_ST_full.set(i, j,
                    ST.get(j, i) * stat_service_phase.get(0, j) / stat_service_phase.get(0, i))
            }
        }
    }

    // tr_ST_exit = sum(-tr_ST, 2)
    val tr_ST_exit_full = Matrix(Tr, 1)
    for (i in 0 until Tr) {
        var sum = 0.0
        for (j in 0 until Tr) {
            sum += tr_ST_full.get(i, j)
        }
        tr_ST_exit_full.set(i, 0, -sum)
    }

    // tr_ST_exit_mat = tr_ST_exit(1:Sc) * ones(1, Sc)
    val tr_ST_exit_mat = Matrix(Sc, Sc)
    for (i in 0 until Sc) {
        for (j in 0 until Sc) {
            tr_ST_exit_mat.set(i, j, tr_ST_exit_full.get(i, 0))
        }
    }

    // TS2 = (TS')' * diag(alfa)
    val TS2 = Matrix(ma * mWait, TS_normalized.getNumCols())
    for (i in 0 until ma * mWait) {
        for (j in 0 until TS_normalized.getNumCols()) {
            TS2.set(i, j, TS_normalized.get(j, i) * alfa.get(0, i))
        }
    }

    // tr_TS_ind2 = sum(TS2, 2)
    val tr_TS_ind2 = Matrix(TS2.getNumRows(), 1)
    for (i in 0 until TS2.getNumRows()) {
        var sum = 0.0
        for (j in 0 until TS2.getNumCols()) {
            sum += TS2.get(i, j)
        }
        tr_TS_ind2.set(i, 0, if (kotlin.math.abs(sum) < 1e-12) 1.0 else sum)
    }

    // tr_TS_ind2 = TS2 ./ (tr_TS_ind2 * ones(1, ma * mWait))
    val tr_TS_ind2_normalized = Matrix(TS2.getNumRows(), Sc)
    for (i in 0 until TS2.getNumRows()) {
        for (j in 0 until Sc) {
            tr_TS_ind2_normalized.set(i, j, TS2.get(i, j) / tr_TS_ind2.get(i, 0))
        }
    }

    // tildeP = [tr_ST_exit_mat .* tr_TS_ind2; zeros(Tr - Sc, Sc)]
    val tildeP = Matrix(Tr, Sc)
    for (i in 0 until Sc) {
        for (j in 0 until Sc) {
            tildeP.set(i, j, tr_ST_exit_mat.get(i, j) * tr_TS_ind2_normalized.get(i, j))
        }
    }

    // Extract non-zero phases
    val tr_ST = Matrix(m_tr_ST, m_tr_ST)
    val tr_start_state = Matrix(1, m_tr_ST)
    val tildeP_reduced = Matrix(m_tr_ST, Sc)

    var idx = 0
    val indexMap = IntArray(Tr)
    for (i in 0 until Tr) {
        if (busy_nz[i]) {
            indexMap[i] = idx
            tr_start_state.set(0, idx, tr_start_state_full.get(0, i))
            idx++
        }
    }

    idx = 0
    for (i in 0 until Tr) {
        if (busy_nz[i]) {
            var idx2 = 0
            for (j in 0 until Tr) {
                if (busy_nz[j]) {
                    tr_ST.set(idx, idx2, tr_ST_full.get(i, j))
                    idx2++
                }
            }
            for (j in 0 until Sc) {
                tildeP_reduced.set(idx, j, tildeP.get(i, j))
            }
            idx++
        }
    }

    // gamma_res = [notbusy_start, tr_start_state, zeros(1, Sc)]
    val gamma_res = Matrix(1, dim_service + m_tr_ST + Sc)
    for (i in 0 until dim_service) {
        gamma_res.set(0, i, notbusy_start.get(0, i))
    }
    for (i in 0 until m_tr_ST) {
        gamma_res.set(0, dim_service + i, tr_start_state.get(0, i))
    }

    // C_res block matrix
    val C_res = Matrix(Tr + m_tr_ST + Sc, Tr + m_tr_ST + Sc)
    setSubMatrix(C_res, 0, 0, ST)
    setSubMatrix(C_res, Tr, Tr, tr_ST)
    setSubMatrix(C_res, Tr, Tr + m_tr_ST, tildeP_reduced)
    setSubMatrix(C_res, Tr + m_tr_ST, Tr + m_tr_ST, wait_Smat)

    // Compute response time percentiles
    return returnPer(gamma_res, C_res, pers)
}

/**
 * Compute element-wise sum of matrix
 */
private fun Matrix.elementSum(): Double {
    var sum = 0.0
    for (i in 0 until getNumRows()) {
        for (j in 0 until getNumCols()) {
            sum += get(i, j)
        }
    }
    return sum
}
