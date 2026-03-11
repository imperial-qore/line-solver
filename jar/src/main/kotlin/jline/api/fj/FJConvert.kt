package jline.api.fj

import jline.api.mam.map_lambda
import jline.api.mam.map_pie
import jline.lang.NetworkStruct
import jline.lang.constant.ProcessType
import jline.lib.fjcodes.FJArrival
import jline.lib.fjcodes.FJService
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

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
 * converting LINE distributions to FJ_codes format. Supports:
 * - Exponential arrivals (single-phase MAP)
 * - General MAP arrivals (multi-phase)
 * - Exponential service (single-phase PH)
 * - General PH/APH/Erlang/HyperExp/Coxian service (multi-phase)
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
        val sourceStation = sn.stations[sn.nodeToStation.get(fjInfo.sourceIdx).toInt()]
        val arrivalProc = sn.proc[sourceStation]?.get(sn.jobclasses[r])
        val lambda = sn.rates.get(fjInfo.sourceIdx, r)

        // Build arrival from process representation
        val arrival: FJArrival
        if (arrivalProc != null && arrivalProc.size() >= 2 && arrivalProc.get(0).getNumRows() > 1) {
            // Multi-phase MAP arrival: use D0, D1 directly
            arrival = convertToFJArrival(arrivalProc.get(0), arrivalProc.get(1))
        } else {
            // Exponential arrival (single-phase)
            arrival = FJArrival(
                lambda = lambda,
                lambda0 = Matrix(1, 1).apply { set(0, 0, -lambda) },
                lambda1 = Matrix(1, 1).apply { set(0, 0, lambda) },
                ArrChoice = 1  // 1 = Exponential
            )
        }
        arrivals.add(arrival)

        // Get service process from first queue
        val firstQueueIdx = fjInfo.queueIndices[0]
        val queueStation = sn.stations[sn.nodeToStation.get(firstQueueIdx).toInt()]
        val serviceProc = sn.proc[queueStation]?.get(sn.jobclasses[r])
        val procType = sn.procid[queueStation]?.get(sn.jobclasses[r])
        val mu_val = sn.mu[queueStation]?.get(sn.jobclasses[r])?.get(0, 0) ?: 1.0

        // Build service from process representation
        val service: FJService
        if (serviceProc != null && serviceProc.size() >= 2 &&
            procType != null && procType != ProcessType.EXP && serviceProc.get(0).getNumRows() > 1) {
            // Multi-phase PH service: extract sub-generator, exit rates, and initial vector
            service = convertToFJServiceFromProc(serviceProc)
        } else {
            // Exponential service (single-phase)
            service = FJService(
                mu = mu_val,
                ST = Matrix(1, 1).apply { set(0, 0, -mu_val) },
                St = Matrix(1, 1).apply { set(0, 0, mu_val) },
                tau_st = Matrix(1, 1).apply { set(0, 0, 1.0) },
                SerChoice = 1  // 1 = Exponential
            )
        }
        services.add(service)

        // Validate stability: lambda < mean service rate
        val meanServiceRate = if (serviceProc != null && serviceProc.size() >= 2 && !serviceProc.get(0).hasNaN()) {
            map_lambda(serviceProc.get(0), serviceProc.get(1))
        } else {
            mu_val
        }
        if (lambda >= meanServiceRate) {
            throw IllegalStateException(
                "Class $r unstable: arrival rate $lambda >= service rate $meanServiceRate"
            )
        }
    }

    return Pair(arrivals, services)
}

/**
 * Convert LINE MAP to FJ arrival format
 *
 * Computes the stationary distribution of the underlying Markov chain
 * to determine the correct arrival rate: lambda = pi * D1 * e.
 *
 * @param D0 D0 matrix (transitions without arrivals)
 * @param D1 D1 matrix (transitions with arrivals)
 * @return FJArrival
 */
fun convertToFJArrival(D0: Matrix, D1: Matrix): FJArrival {
    val n = D0.getNumRows()

    // Compute arrival rate using MAP stationary distribution
    // pi is the stationary distribution of D0+D1, lambda = pi * D1 * e
    val lambda = map_lambda(D0, D1)

    return FJArrival(
        lambda = lambda,
        lambda0 = D0,
        lambda1 = D1,
        ArrChoice = if (n == 1) 1 else 2  // 1=Exp, 2=MAP
    )
}

/**
 * Convert LINE PH to FJ service format
 *
 * Uses the PH initial vector (pie from MAP representation) as the
 * initial probability for phase selection. Mean service rate is
 * computed as mu = 1 / (tau * (-S)^{-1} * e).
 *
 * @param S Sub-generator matrix (negative diagonal = phase rates)
 * @param s Exit rate vector (s = -S * e for absorbing PH)
 * @param tau Initial probability vector
 * @return FJService
 */
fun convertToFJService(S: Matrix, s: Matrix, tau: Matrix): FJService {
    val n = S.getNumRows()

    // Compute mean service rate: mu = 1 / (tau * (-S)^{-1} * e)
    val negSinv = S.scale(-1.0).inv()
    val ones = Matrix.ones(n, 1)
    val meanServiceTime = tau.mult(negSinv).mult(ones).get(0, 0)
    val mu = 1.0 / meanServiceTime

    return FJService(
        mu = mu,
        ST = S,
        St = s,
        tau_st = tau,
        SerChoice = if (n == 1) 1 else 2  // 1=Exp, 2=PH
    )
}

/**
 * Convert LINE process cell {D0, D1} to FJ service format.
 *
 * Extracts the sub-generator S=D0, exit rate vector s=-D0*e (or D1*e),
 * and initial probability vector tau from the MAP/PH representation.
 *
 * @param proc MatrixCell containing {D0, D1}
 * @return FJService
 */
fun convertToFJServiceFromProc(proc: MatrixCell): FJService {
    val D0 = proc.get(0)  // Sub-generator S
    val D1 = proc.get(1)  // Completion matrix
    val n = D0.getNumRows()

    // Exit rate vector: s = -S * e (row sums of -D0)
    val ones = Matrix.ones(n, 1)
    val s = D0.scale(-1.0).mult(ones)

    // Initial probability vector: use MAP pie (embedded stationary distribution)
    val tau = map_pie(D0, D1)

    // Mean service rate
    val negSinv = D0.scale(-1.0).inv()
    val meanServiceTime = tau.mult(negSinv).mult(ones).get(0, 0)
    val mu = 1.0 / meanServiceTime

    return FJService(
        mu = mu,
        ST = D0,
        St = s,
        tau_st = tau,
        SerChoice = if (n == 1) 1 else 2  // 1=Exp, 2=PH
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
