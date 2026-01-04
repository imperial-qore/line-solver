/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.nc.analyzers

import jline.api.lossn.lossn_erlangfp
import jline.lang.NetworkStruct
import jline.solvers.SolverOptions
import jline.solvers.nc.NCResult
import jline.util.matrix.Matrix

/**
 * Analyzes open loss networks with FCR using Erlang fixed-point approximation.
 *
 * This analyzer handles open queueing networks with a single multiclass
 * Delay node inside a Finite Capacity Region (FCR) with DROP policy.
 *
 * @param sn Network structure
 * @param options Solver options
 * @return NCResult containing performance metrics
 */
fun solver_nc_lossn_analyzer(sn: NetworkStruct, options: SolverOptions): NCResult {
    val Tstart = System.nanoTime()
    val K = sn.nclasses  // number of classes
    val M = sn.nstations

    // 1. Extract arrival rates from Source
    val nu = DoubleArray(K)
    for (r in 0..<K) {
        val sourceIdx = sn.refstat.get(r).toInt()
        nu[r] = sn.rates.get(sourceIdx, r)  // arrival rate
    }

    // 2. Find delay station in FCR and extract constraints
    val regionMatrix = sn.region.get(0)
    var delayIdx = -1
    for (i in 0..<M) {
        var hasConstraint = false
        // Check per-class constraints
        for (r in 0..<K) {
            if (regionMatrix.get(i, r) >= 0) {
                hasConstraint = true
                break
            }
        }
        // Check global constraint
        if (regionMatrix.get(i, K) >= 0) {
            hasConstraint = true
        }
        if (hasConstraint) {
            delayIdx = i
            break
        }
    }

    val globalMax = regionMatrix.get(delayIdx, K)
    val classMax = DoubleArray(K)
    for (r in 0..<K) {
        classMax[r] = regionMatrix.get(delayIdx, r)
    }

    // Handle unbounded constraints (replace -1 with large value for Erlang)
    val globalMaxVal = if (globalMax < 0) 1e6 else globalMax
    for (r in 0..<K) {
        if (classMax[r] < 0) {
            classMax[r] = 1e6
        }
    }

    // 3. Build A matrix (J x K) where J = K+1 links
    // Link 0: global constraint (all classes contribute)
    // Links 1..K: per-class constraints (only class r contributes to link r)
    val J = K + 1
    val A = Matrix(J, K)
    A.fill(0.0)
    // Global link: all classes contribute
    for (r in 0..<K) {
        A.set(0, r, 1.0)
    }
    // Per-class links
    for (r in 0..<K) {
        A.set(r + 1, r, 1.0)
    }

    // 4. Build C vector (J x 1)
    val C_vec = Matrix(J, 1)
    C_vec.set(0, 0, globalMaxVal)
    for (r in 0..<K) {
        C_vec.set(r + 1, 0, classMax[r])
    }

    // 5. Call lossn_erlangfp
    val lossnResult = lossn_erlangfp(Matrix(nu), A, C_vec)
    val QLen = lossnResult.qLen  // effective throughput after loss
    val niter = lossnResult.niter

    // 6. Convert to standard outputs
    val Q = Matrix(M, K)
    Q.fill(0.0)
    val U = Matrix(M, K)
    U.fill(0.0)
    val T = Matrix(M, K)
    T.fill(0.0)
    val R = Matrix(M, K)
    R.fill(0.0)

    // At delay node: QLen is effective throughput (after loss)
    for (r in 0..<K) {
        val effTput = QLen.get(r)
        T.set(delayIdx, r, effTput)  // effective throughput = arrival rate * (1-loss)
        val mu_r = sn.rates.get(delayIdx, r)  // service rate at delay
        Q.set(delayIdx, r, effTput / mu_r)  // Little's law: Q = X * S
        R.set(delayIdx, r, 1.0 / mu_r)  // response time = service time (infinite server)
        U.set(delayIdx, r, effTput / mu_r)  // "utilization" for delay
    }

    // System throughput per class (row vector 1 x K)
    val X = Matrix(1, K)
    for (r in 0..<K) {
        X.set(0, r, QLen.get(r))
    }

    // Cycle time not applicable for open networks
    val CN = Matrix(1, K)
    CN.fill(0.0)

    // Build result
    val res = NCResult()
    res.QN = Q
    res.UN = U
    res.RN = R
    res.TN = T
    res.XN = X
    res.CN = CN
    res.lG = Double.NaN  // no normalizing constant for loss networks
    res.it = niter
    res.iter = niter
    res.method = "erlangfp"
    res.runtime = (System.nanoTime() - Tstart).toDouble() / 1000000000.0

    return res
}
