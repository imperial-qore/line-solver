package jline.api.fj

import jline.lang.NetworkStruct
import jline.lang.constant.ProcessType
import jline.lib.fjcodes.FJArrival
import jline.lib.fjcodes.FJService
import jline.util.matrix.Matrix

/**
 * Convert LINE distributions to FJ_codes format
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */

/**
 * Extract Fork-Join parameters from network structure
 *
 * Extracts arrival and service processes for each class from the network,
 * converting LINE distributions to FJ_codes format.
 *
 * @param sn Network structure
 * @param fjInfo Fork-Join topology information
 * @return Pair of (arrivals per class, services per class)
 */
fun extractFJParams(
    sn: NetworkStruct,
    fjInfo: FJInfo
): Pair<List<FJArrival>, List<FJService>> {
    val arrivals = mutableListOf<FJArrival>()
    val services = mutableListOf<FJService>()

    // For each class, extract arrival and service parameters
    for (r in 0 until sn.nclasses) {
        // Get arrival process from Source
        val sourceNode = sn.nodes[fjInfo.sourceIdx]
        val sourceStation = sn.stations[sn.nodeToStation.get(fjInfo.sourceIdx).toInt()]

        // Extract arrival rate
        val arrivalProc = sn.proc[sourceStation]?.get(sn.jobclasses[r])
        val lambda = sn.rates.get(fjInfo.sourceIdx, r)

        // For MVP: assume exponential arrivals
        // TODO: Support MAP arrivals
        val arrival = FJArrival(
            lambda = lambda,
            lambda0 = Matrix(1, 1).apply { set(0, 0, -lambda) },
            lambda1 = Matrix(1, 1).apply { set(0, 0, lambda) },
            ArrChoice = 1  // 1 = Exponential
        )
        arrivals.add(arrival)

        // Get service process from first queue
        val firstQueueIdx = fjInfo.queueIndices[0]
        val queueStation = sn.stations[sn.nodeToStation.get(firstQueueIdx).toInt()]

        // Extract service parameters
        val serviceProc = sn.proc[queueStation]?.get(sn.jobclasses[r])
        val mu_val = sn.mu[queueStation]?.get(sn.jobclasses[r])?.get(0, 0) ?: 1.0
        val procType = sn.procid[queueStation]?.get(sn.jobclasses[r])

        // For MVP: assume exponential service
        // TODO: Support PH service distributions
        val service = FJService(
            mu = mu_val,
            ST = Matrix(1, 1).apply { set(0, 0, -mu_val) },
            St = Matrix(1, 1).apply { set(0, 0, mu_val) },
            tau_st = Matrix(1, 1).apply { set(0, 0, 1.0) },
            SerChoice = 1  // 1 = Exponential
        )
        services.add(service)

        // Validate stability: lambda < mu
        if (lambda >= mu_val) {
            throw IllegalStateException(
                "Class $r unstable: arrival rate $lambda >= service rate $mu_val"
            )
        }
    }

    return Pair(arrivals, services)
}

/**
 * Convert LINE MAP to FJ arrival format
 *
 * @param D0 D0 matrix (transitions without arrivals)
 * @param D1 D1 matrix (transitions with arrivals)
 * @return FJArrival
 */
fun convertToFJArrival(D0: Matrix, D1: Matrix): FJArrival {
    // Compute arrival rate
    val pi = Matrix(1, D0.getNumRows())
    // TODO: Compute stationary distribution properly
    // For now, use uniform distribution
    val n = D0.getNumRows()
    for (i in 0 until n) {
        pi.set(0, i, 1.0 / n)
    }

    val lambda = pi.mult(D1).elementSum()

    return FJArrival(
        lambda = lambda,
        lambda0 = D0,
        lambda1 = D1,
        ArrChoice = if (D0.getNumRows() == 1) 1 else 2  // 1=Exp, 2=MAP
    )
}

/**
 * Convert LINE PH to FJ service format
 *
 * @param S Sub-generator matrix
 * @param s Exit rate vector
 * @param tau Initial probability vector
 * @return FJService
 */
fun convertToFJService(S: Matrix, s: Matrix, tau: Matrix): FJService {
    // Compute mean service rate
    val pi = Matrix(1, S.getNumRows())
    // TODO: Compute stationary distribution properly
    val n = S.getNumRows()
    for (i in 0 until n) {
        pi.set(0, i, 1.0 / n)
    }

    val mu = 1.0 / pi.mult(S.inv().scale(-1.0)).mult(Matrix.ones(n, 1)).get(0, 0)

    return FJService(
        mu = mu,
        ST = S,
        St = s,
        tau_st = tau,
        SerChoice = if (S.getNumRows() == 1) 1 else 2  // 1=Exp, 2=PH
    )
}

/**
 * Helper: compute element sum
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
