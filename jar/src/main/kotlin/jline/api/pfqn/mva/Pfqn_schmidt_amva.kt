/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.mva

import jline.io.Ret
import jline.lang.constant.SchedStrategy
import jline.util.PopulationLattice.hashpop
import jline.util.PopulationLattice.pprod
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.CombinatoricsUtils
import kotlin.math.max
import kotlin.math.min
import kotlin.math.pow

/**
 * Schmidt MVA algorithm for multi-class FCFS queueing networks.
 *
 * Computes performance metrics for product-form queueing networks using the
 * Schmidt MVA algorithm. This is an exact algorithm for closed queueing networks
 * with class-dependent FCFS scheduling.
 *
 * Reference: R. Schmidt, "An approximate MVA algorithm for exponential, class-dependent
 * multiple server stations," Performance Evaluation, vol. 29, no. 4, pp. 245-254, 1997.
 *
 * @param rates Service rate matrix (M x R) - rates per station and class
 * @param N Population vector (1 x R)
 * @param S Number of servers at each station (M x 1)
 * @param v Visit ratio matrix (M x R)
 * @param sched List of scheduling strategies for each station
 * @return pfqnAMVASchmidt result containing throughput, queue length, utilization, response time, and probabilities
 */
fun pfqn_schmidt(
    rates: Matrix,
    N: Matrix,
    S: Matrix,
    v: Matrix,
    sched: List<SchedStrategy>
): Ret.pfqnAMVASchmidt {
    val M = rates.numRows
    val R = rates.numCols

    // Convert rates to service demands (D = 1/rate)
    val D = Matrix(M, R)
    for (i in 0 until M) {
        for (r in 0 until R) {
            if (rates[i, r] == 0.0) {
                D[i, r] = 0.0
            } else {
                D[i, r] = 1.0 / rates[i, r]
            }
        }
    }

    val closedClasses = IntArray(R) { it }

    val XN = Matrix.zeros(1, R)
    val UN = Matrix.zeros(M, R)
    val CN = Matrix.zeros(M, R)
    val QN = Matrix.zeros(M, R)
    val PN = mutableListOf<MutableMap<Any, Double>>()

    for (i in 0 until M) {
        PN.add(mutableMapOf())
    }

    val C = closedClasses.size
    val Nc = N

    // Calculate products for hashing
    val prods = Matrix.zeros(1, C)
    for (r in 0 until C) {
        var prod = 1.0
        for (i in 0 until r) {
            prod *= (Nc[i] + 1)
        }
        prods[0, r] = prod
    }

    var kvec = pprod(Nc)
    val totalStates = prod(Nc).toInt()

    // Initialize data structures
    val L = mutableListOf<Matrix>()
    val Pc = mutableListOf<Matrix?>()

    for (i in 0 until M) {
        L.add(Matrix.zeros(R, totalStates))
        when (sched[i]) {
            SchedStrategy.INF -> Pc.add(null)
            SchedStrategy.PS -> {
                val singleServerPS = S[i] == 1.0
                if (singleServerPS) {
                    Pc.add(null)
                } else {
                    Pc.add(Matrix.zeros(S[i].toInt(), totalStates))
                }
            }
            SchedStrategy.FCFS -> {
                var classIndependent = true
                for (r in 1 until R) {
                    if (D[i, r] != D[i, 0]) {
                        classIndependent = false
                        break
                    }
                }
                val isSingleServer = S[i] == 1.0

                if (classIndependent) {
                    if (isSingleServer) {
                        Pc.add(null)
                    } else {
                        Pc.add(Matrix.zeros(S[i].toInt(), totalStates))
                    }
                } else {
                    Pc.add(Matrix.zeros(totalStates, totalStates))
                }
            }
            else -> Pc.add(null)
        }
    }

    val x = Array(M) { Matrix.zeros(C, totalStates) }
    val w = Array(M) { Matrix.zeros(C, totalStates) }

    // Initialize base case
    for (ist in 0 until M) {
        Pc[ist]?.set(0, hashpop(kvec, Nc), 1.0)
    }

    kvec = pprod(kvec, Nc)
    var hkvec = hashpop(kvec, Nc)

    while (allGE(kvec, 0) && allLE(kvec, Nc)) {
        val kprods = Matrix.zeros(1, C)
        for (r in 0 until C) {
            // kprods[r] = prod(kvec(0:r-1) + 1), matches MATLAB: prod(kvec(1:r-1)+1)
            var prod = 1.0
            for (i in 0 until r) {
                prod *= (kvec[i] + 1)
            }
            kprods[r] = prod
        }

        for (i in 0 until M) {
            for (c in 0 until C) {
                val ns = S[i]
                hkvec = hashpop(kvec, Nc)

                val kvec_c = oner(kvec, c)
                val hkvec_c = hashpop(kvec_c, Nc)

                if (kvec[c] == 0.0) continue

                when (sched[i]) {
                    SchedStrategy.INF -> {
                        w[i][c, hkvec] = D[i, c]
                    }
                    SchedStrategy.FCFS -> {
                        var classIndependent = true
                        for (r in 1 until R) {
                            if (D[i, r] != D[i, 0]) {
                                classIndependent = false
                                break
                            }
                        }

                        if (!classIndependent) {
                            if (ns == 1.0) {
                                var totalQueueLengthKvecC = 0.0
                                for (r in 0 until R) {
                                    totalQueueLengthKvecC += L[i][r, hkvec_c]
                                }
                                w[i][c, hkvec] = D[i, c] * (1 + totalQueueLengthKvecC)
                            } else {
                                var waitTime = 0.0
                                var nvec = pprod(kvec)
                                while (nvec != null && nvec.value() >= 0) {
                                    if (nvec[c] > 0) {
                                        val nvec_c = oner(nvec, c)
                                        val hnvec_c = hashpop(nvec_c, Nc, C, prods)

                                        val Bcn = getBcn(D, i, c, nvec, C, ns)
                                        val prob = Pc[i]!![hnvec_c, hkvec_c]
                                        waitTime += (Bcn * prob)
                                    }
                                    nvec = pprod(nvec, kvec)
                                }
                                w[i][c, hkvec] = waitTime
                            }
                        }
                        // Fall through to PS for class-independent case
                        if (classIndependent) {
                            if (ns == 1.0) {
                                var totalQueueLengthKvecC = 0.0
                                for (r in 0 until R) {
                                    totalQueueLengthKvecC += L[i][r, hkvec_c]
                                }
                                w[i][c, hkvec] = D[i, c] * (1 + totalQueueLengthKvecC)
                            } else {
                                var totalQueueLengthKvecC = 0.0
                                for (r in 0 until R) {
                                    totalQueueLengthKvecC += L[i][r, hkvec_c]
                                }
                                var weightedQueueLength = 0.0
                                for (j in 0 until (ns - 1).toInt()) {
                                    weightedQueueLength += ((ns - j - 1) * Pc[i]!![j, hkvec_c])
                                }
                                val totalWaitTime = (D[i, c] / ns) * (1 + totalQueueLengthKvecC + weightedQueueLength)
                                w[i][c, hkvec] = totalWaitTime
                            }
                        }
                    }
                    SchedStrategy.PS -> {
                        if (ns == 1.0) {
                            var totalQueueLengthKvecC = 0.0
                            for (r in 0 until R) {
                                totalQueueLengthKvecC += L[i][r, hkvec_c]
                            }
                            w[i][c, hkvec] = D[i, c] * (1 + totalQueueLengthKvecC)
                        } else {
                            var totalQueueLengthKvecC = 0.0
                            for (r in 0 until R) {
                                totalQueueLengthKvecC += L[i][r, hkvec_c]
                            }
                            var weightedQueueLength = 0.0
                            for (j in 0 until (ns - 1).toInt()) {
                                weightedQueueLength += ((ns - j - 1) * Pc[i]!![j, hkvec_c])
                            }
                            val totalWaitTime = (D[i, c] / ns) * (1 + totalQueueLengthKvecC + weightedQueueLength)
                            w[i][c, hkvec] = totalWaitTime
                        }
                    }
                    else -> {}
                }
            }
        }

        // Compute throughputs
        for (c in 0 until C) {
            var denom = 0.0
            for (i in 0 until M) {
                denom += (v[i, c] * w[i][c, hkvec])
            }

            for (i in 0 until M) {
                if (denom > 0) {
                    x[i][c, hkvec] = v[i, c] * kvec[c] / denom
                } else {
                    x[i][c, hkvec] = 0.0
                }
            }
        }

        // Update queue lengths
        for (i in 0 until M) {
            for (c in 0 until C) {
                L[i][c, hkvec] = x[i][c, hkvec] * w[i][c, hkvec]
            }
            val ns = S[i].toInt()
            var K_j = 0
            for (r in 0 until R) {
                if (v[i, r] > 0) {
                    K_j += Nc[r].toInt()
                }
            }

            when (sched[i]) {
                SchedStrategy.FCFS -> {
                    var classIndependent = true
                    for (r in 1 until R) {
                        if (D[i, r] != D[i, 0]) {
                            classIndependent = false
                            break
                        }
                    }

                    if (!classIndependent) {
                        var nvec = pprod(kvec)
                        nvec = pprod(nvec, kvec)
                        var sumOfAllProbs = 0.0
                        while (nvec != null && nvec.value() >= 0) {
                            val hnvec = hashpop(nvec, Nc, C, prods)

                            var prob = 0.0
                            for (r in 0 until C) {
                                if (nvec[r] > 0) {
                                    val hnvec_c = hashpop(oner(nvec, r), Nc, C, prods)
                                    val hkvec_c = hashpop(oner(kvec, r), Nc)
                                    val Bcn = getBcn(D, i, r, nvec, C, ns.toDouble())
                                    val capacity_inv = 1 / nvec.elementSum()
                                    val x_ir = x[i][r, hkvec]
                                    val prob_c = Pc[i]!![hnvec_c, hkvec_c]
                                    val classProb = (Bcn * capacity_inv * x_ir * prob_c)
                                    prob += classProb
                                }
                            }
                            Pc[i]!![hnvec, hkvec] = prob
                            sumOfAllProbs += prob
                            nvec = pprod(nvec, kvec)
                        }
                        Pc[i]!![0, hkvec] = max(1e-12, 1.0 - sumOfAllProbs)
                    }
                    // For class-independent, fall through to PS probability update
                    if (!classIndependent) {
                        // Already handled above
                    } else if (ns > 1) {
                        var meanQueueLength = 0.0
                        for (r in 0 until R) {
                            meanQueueLength += L[i][r, hkvec]
                        }
                        for (n in 1 until ns) {
                            val x1 = CombinatoricsUtils.binomialCoefficient(K_j, n).toDouble()
                            val frac = meanQueueLength / K_j
                            val x2 = frac.pow(n.toDouble())
                            val x3 = (1 - frac).pow((K_j - n).toDouble())
                            val prob = x1 * x2 * x3
                            Pc[i]!![n, hkvec] = prob
                        }
                        var sum1 = 0.0
                        var sum2 = 0.0
                        for (r in 0 until R) {
                            sum1 += D[i, r] * x[i][r, hkvec]
                        }
                        for (n in 0 until ns) {
                            sum2 += (ns - n) * Pc[i]!![n, hkvec]
                        }
                        val term = (sum1 + sum2) / ns
                        Pc[i]!![0, hkvec] = max(1e-12, 1 - term)
                    }
                }
                SchedStrategy.PS -> {
                    if (ns > 1) {
                        var meanQueueLength = 0.0
                        for (r in 0 until R) {
                            meanQueueLength += L[i][r, hkvec]
                        }
                        for (n in 1 until ns) {
                            val x1 = CombinatoricsUtils.binomialCoefficient(K_j, n).toDouble()
                            val frac = meanQueueLength / K_j
                            val x2 = frac.pow(n.toDouble())
                            val x3 = (1 - frac).pow((K_j - n).toDouble())
                            val prob = x1 * x2 * x3
                            Pc[i]!![n, hkvec] = prob
                        }
                        var sum1 = 0.0
                        var sum2 = 0.0
                        for (r in 0 until R) {
                            sum1 += D[i, r] * x[i][r, hkvec]
                        }
                        for (n in 0 until ns) {
                            sum2 += (ns - n) * Pc[i]!![n, hkvec]
                        }
                        val term = (sum1 + sum2) / ns
                        Pc[i]!![0, hkvec] = max(1e-12, 1 - term)
                    }
                }
                else -> {}
            }
        }

        kvec = pprod(kvec, Nc)
    }

    val hkvecFinal = hkvec

    // Extract final results
    for (c in 0 until C) {
        var totalResponseTime = 0.0
        for (i in 0 until M) {
            totalResponseTime += w[i][c, hkvecFinal]
        }
        XN[c] = Nc[c] / totalResponseTime
    }

    for (m in 0 until M) {
        for (c in 0 until C) {
            UN[m, c] = (D[m, c] * XN[c]) / S[m]
        }
    }

    for (m in 0 until M) {
        for (c in 0 until C) {
            CN[m, c] = w[m][c, hkvecFinal]
        }
    }

    for (m in 0 until M) {
        val Lmat = L[m]
        for (c in 0 until C) {
            QN[m, c] = Lmat[c, hkvecFinal]
        }
    }

    @Suppress("UNCHECKED_CAST")
    return Ret.pfqnAMVASchmidt(XN, QN, UN, CN, PN as List<MutableMap<Any, Double>>)
}

/**
 * Extended Schmidt MVA algorithm with queue-aware alpha corrections.
 *
 * A queue-aware version of the Schmidt algorithm that precomputes alpha values
 * for improved accuracy in networks with class-dependent FCFS scheduling.
 *
 * @param rates Service rate matrix (M x R)
 * @param N Population vector (1 x R)
 * @param S Number of servers at each station (M x 1)
 * @param v Visit ratio matrix (M x R)
 * @param sched List of scheduling strategies for each station
 * @return pfqnAMVASchmidt result
 */
fun pfqn_schmidt_ext(
    rates: Matrix,
    N: Matrix,
    S: Matrix,
    v: Matrix,
    sched: List<SchedStrategy>
): Ret.pfqnAMVASchmidt {
    val M = rates.numRows
    val R = rates.numCols

    // Convert rates to service demands (D = 1/rate)
    val D = Matrix(M, R)
    for (i in 0 until M) {
        for (r in 0 until R) {
            if (rates[i, r] == 0.0 || rates[i, r].isNaN()) {
                D[i, r] = 0.0
            } else {
                D[i, r] = 1.0 / rates[i, r]
            }
        }
    }

    // Precompute alphas for FCFS stations
    val alphas = mutableMapOf<Pair<Int, Int>, Matrix>()

    for (i in 0 until M) {
        if (sched[i] == SchedStrategy.FCFS) {
            for (r in 0 until R) {
                val rates_mod = Matrix(M, R + 1)
                val N_mod = Matrix(1, R + 1)
                val v_mod = Matrix(M, R + 1)

                for (k in 0 until R) {
                    N_mod[k] = if (k == r) N[k] - 1 else N[k]
                    for (j in 0 until M) {
                        rates_mod[j, k] = rates[j, k]
                        v_mod[j, k] = v[j, k]
                        rates_mod[j, R] = if (i == j) rates[j, r] else 0.0
                        v_mod[j, R] = if (i == j) 1.0 else 0.0
                    }
                }
                N_mod[0, R] = 1.0

                val result = pfqn_schmidt(rates_mod, N_mod, S, v_mod, sched)
                alphas[Pair(i, r)] = result.U
            }
        }
    }

    val closedClasses = IntArray(R) { it }

    val XN = Matrix.zeros(1, R)
    val UN = Matrix.zeros(M, R)
    val CN = Matrix.zeros(M, R)
    val QN = Matrix.zeros(M, R)
    val PN = mutableListOf<MutableMap<Any, Double>>()

    for (i in 0 until M) {
        PN.add(mutableMapOf())
    }

    val C = closedClasses.size
    val Nc = N

    val prods = Matrix.zeros(1, C)
    for (r in 0 until C) {
        var prod = 1.0
        for (i in 0 until r) {
            prod *= (Nc[i] + 1)
        }
        prods[0, r] = prod
    }

    var kvec = pprod(Nc)
    val totalStates = prod(Nc).toInt()

    val L = mutableListOf<Matrix>()
    val Pc = mutableListOf<Matrix?>()

    for (i in 0 until M) {
        L.add(Matrix.zeros(R, totalStates))
        when (sched[i]) {
            SchedStrategy.INF -> Pc.add(null)
            SchedStrategy.PS -> {
                val singleServerPS = S[i] == 1.0
                if (singleServerPS) {
                    Pc.add(null)
                } else {
                    Pc.add(Matrix.zeros(S[i].toInt(), totalStates))
                }
            }
            SchedStrategy.FCFS -> {
                var classIndependent = true
                for (r in 1 until R) {
                    if (D[i, r] != D[i, 0]) {
                        classIndependent = false
                        break
                    }
                }
                val isSingleServer = S[i] == 1.0

                if (classIndependent) {
                    if (isSingleServer) {
                        Pc.add(null)
                    } else {
                        Pc.add(Matrix.zeros(S[i].toInt(), totalStates))
                    }
                } else {
                    Pc.add(Matrix.zeros(totalStates, totalStates))
                }
            }
            else -> Pc.add(null)
        }
    }

    val x = Array(M) { Matrix.zeros(C, totalStates) }
    val w = Array(M) { Matrix.zeros(C, totalStates) }

    for (ist in 0 until M) {
        Pc[ist]?.set(0, hashpop(kvec, Nc), 1.0)
    }

    kvec = pprod(kvec, Nc)
    var hkvec = hashpop(kvec, Nc)

    while (allGE(kvec, 0) && allLE(kvec, Nc)) {
        val kprods = Matrix.zeros(1, C)
        for (r in 0 until C) {
            // kprods[r] = prod(kvec(0:r-1) + 1), matches MATLAB: prod(kvec(1:r-1)+1)
            var prod = 1.0
            for (i in 0 until r) {
                prod *= (kvec[i] + 1)
            }
            kprods[r] = prod
        }

        for (i in 0 until M) {
            for (c in 0 until C) {
                val ns = S[i]
                hkvec = hashpop(kvec, Nc)

                val kvec_c = oner(kvec, c)
                val hkvec_c = hashpop(kvec_c, Nc)

                if (kvec[c] == 0.0) continue

                when (sched[i]) {
                    SchedStrategy.INF -> {
                        w[i][c, hkvec] = D[i, c]
                    }
                    SchedStrategy.FCFS -> {
                        var classIndependent = true
                        for (r in 1 until R) {
                            if (D[i, r] != D[i, 0]) {
                                classIndependent = false
                                break
                            }
                        }

                        if (!classIndependent) {
                            if (ns == 1.0) {
                                var totalQueueLengthKvecC = 0.0
                                for (r in 0 until R) {
                                    totalQueueLengthKvecC += L[i][r, hkvec_c]
                                }
                                w[i][c, hkvec] = D[i, c] * (1 + totalQueueLengthKvecC)
                            } else {
                                var waitTime = 0.0
                                var nvec = pprod(kvec)
                                while (nvec != null && nvec.value() >= 0) {
                                    if (nvec[c] > 0) {
                                        val nvec_c = oner(nvec, c)
                                        val hnvec_c = hashpop(nvec_c, Nc, C, prods)

                                        // Use extended Bcn with alpha
                                        val Bcn = if (Nc[c] > 1.0) {
                                            getBcnExt(alphas[Pair(i, c)]!!, D, i, c, nvec, C, ns)
                                        } else {
                                            getBcn(D, i, c, nvec, C, ns)
                                        }

                                        val prob = Pc[i]!![hnvec_c, hkvec_c]
                                        waitTime += (Bcn * prob)
                                    }
                                    nvec = pprod(nvec, kvec)
                                }
                                w[i][c, hkvec] = waitTime
                            }
                        }
                        // Handle class-independent case like PS
                        if (classIndependent) {
                            if (ns == 1.0) {
                                var totalQueueLengthKvecC = 0.0
                                for (r in 0 until R) {
                                    totalQueueLengthKvecC += L[i][r, hkvec_c]
                                }
                                w[i][c, hkvec] = D[i, c] * (1 + totalQueueLengthKvecC)
                            } else {
                                var totalQueueLengthKvecC = 0.0
                                for (r in 0 until R) {
                                    totalQueueLengthKvecC += L[i][r, hkvec_c]
                                }
                                var weightedQueueLength = 0.0
                                for (j in 0 until (ns - 1).toInt()) {
                                    weightedQueueLength += ((ns - j - 1) * Pc[i]!![j, hkvec_c])
                                }
                                val totalWaitTime = (D[i, c] / ns) * (1 + totalQueueLengthKvecC + weightedQueueLength)
                                w[i][c, hkvec] = totalWaitTime
                            }
                        }
                    }
                    SchedStrategy.PS -> {
                        if (ns == 1.0) {
                            var totalQueueLengthKvecC = 0.0
                            for (r in 0 until R) {
                                totalQueueLengthKvecC += L[i][r, hkvec_c]
                            }
                            w[i][c, hkvec] = D[i, c] * (1 + totalQueueLengthKvecC)
                        } else {
                            var totalQueueLengthKvecC = 0.0
                            for (r in 0 until R) {
                                totalQueueLengthKvecC += L[i][r, hkvec_c]
                            }
                            var weightedQueueLength = 0.0
                            for (j in 0 until (ns - 1).toInt()) {
                                weightedQueueLength += ((ns - j - 1) * Pc[i]!![j, hkvec_c])
                            }
                            val totalWaitTime = (D[i, c] / ns) * (1 + totalQueueLengthKvecC + weightedQueueLength)
                            w[i][c, hkvec] = totalWaitTime
                        }
                    }
                    else -> {}
                }
            }
        }

        // Compute throughputs
        for (c in 0 until C) {
            var denom = 0.0
            for (i in 0 until M) {
                denom += (v[i, c] * w[i][c, hkvec])
            }

            for (i in 0 until M) {
                if (denom > 0) {
                    x[i][c, hkvec] = v[i, c] * kvec[c] / denom
                } else {
                    x[i][c, hkvec] = 0.0
                }
            }
        }

        // Update queue lengths and state probabilities
        for (i in 0 until M) {
            for (c in 0 until C) {
                L[i][c, hkvec] = x[i][c, hkvec] * w[i][c, hkvec]
            }
            val ns = S[i].toInt()
            var K_j = 0
            for (r in 0 until R) {
                if (v[i, r] > 0) {
                    K_j += Nc[r].toInt()
                }
            }

            when (sched[i]) {
                SchedStrategy.FCFS -> {
                    var classIndependent = true
                    for (r in 1 until R) {
                        if (D[i, r] != D[i, 0]) {
                            classIndependent = false
                            break
                        }
                    }

                    if (!classIndependent) {
                        var nvec = pprod(kvec)
                        nvec = pprod(nvec, kvec)
                        var sumOfAllProbs = 0.0
                        while (nvec != null && nvec.value() >= 0) {
                            val hnvec = hashpop(nvec, Nc, C, prods)

                            var prob = 0.0
                            for (r in 0 until C) {
                                if (nvec[r] > 0) {
                                    val hnvec_c = hashpop(oner(nvec, r), Nc, C, prods)
                                    val hkvec_c = hashpop(oner(kvec, r), Nc)
                                    val Bcn = getBcnExt(alphas[Pair(i, r)]!!, D, i, r, nvec, C, ns.toDouble())
                                    val capacity_inv = 1 / nvec.elementSum()
                                    val x_ir = x[i][r, hkvec]
                                    val prob_c = Pc[i]!![hnvec_c, hkvec_c]
                                    val classProb = (Bcn * capacity_inv * x_ir * prob_c)
                                    prob += classProb
                                }
                            }
                            Pc[i]!![hnvec, hkvec] = prob
                            sumOfAllProbs += prob
                            nvec = pprod(nvec, kvec)
                        }
                        Pc[i]!![0, hkvec] = max(1e-12, 1.0 - sumOfAllProbs)
                    } else if (ns > 1) {
                        var meanQueueLength = 0.0
                        for (r in 0 until R) {
                            meanQueueLength += L[i][r, hkvec]
                        }
                        for (n in 1 until ns) {
                            val x1 = CombinatoricsUtils.binomialCoefficient(K_j, n).toDouble()
                            val frac = meanQueueLength / K_j
                            val x2 = frac.pow(n.toDouble())
                            val x3 = (1 - frac).pow((K_j - n).toDouble())
                            val prob = x1 * x2 * x3
                            Pc[i]!![n, hkvec] = prob
                        }
                        var sum1 = 0.0
                        var sum2 = 0.0
                        for (r in 0 until R) {
                            sum1 += D[i, r] * x[i][r, hkvec]
                        }
                        for (n in 0 until ns) {
                            sum2 += (ns - n) * Pc[i]!![n, hkvec]
                        }
                        val term = (sum1 + sum2) / ns
                        Pc[i]!![0, hkvec] = max(1e-12, 1 - term)
                    }
                }
                SchedStrategy.PS -> {
                    if (ns > 1) {
                        var meanQueueLength = 0.0
                        for (r in 0 until R) {
                            meanQueueLength += L[i][r, hkvec]
                        }
                        for (n in 1 until ns) {
                            val x1 = CombinatoricsUtils.binomialCoefficient(K_j, n).toDouble()
                            val frac = meanQueueLength / K_j
                            val x2 = frac.pow(n.toDouble())
                            val x3 = (1 - frac).pow((K_j - n).toDouble())
                            val prob = x1 * x2 * x3
                            Pc[i]!![n, hkvec] = prob
                        }
                        var sum1 = 0.0
                        var sum2 = 0.0
                        for (r in 0 until R) {
                            sum1 += D[i, r] * x[i][r, hkvec]
                        }
                        for (n in 0 until ns) {
                            sum2 += (ns - n) * Pc[i]!![n, hkvec]
                        }
                        val term = (sum1 + sum2) / ns
                        Pc[i]!![0, hkvec] = max(1e-12, 1 - term)
                    }
                }
                else -> {}
            }
        }

        kvec = pprod(kvec, Nc)
    }

    val hkvecFinal = hkvec

    // Extract final results
    for (c in 0 until C) {
        var totalResponseTime = 0.0
        for (i in 0 until M) {
            totalResponseTime += w[i][c, hkvecFinal]
        }
        XN[c] = Nc[c] / totalResponseTime
    }

    for (m in 0 until M) {
        for (c in 0 until C) {
            UN[m, c] = (D[m, c] * XN[c]) / S[m]
        }
    }

    for (m in 0 until M) {
        for (c in 0 until C) {
            CN[m, c] = w[m][c, hkvecFinal]
        }
    }

    for (m in 0 until M) {
        val Lmat = L[m]
        for (c in 0 until C) {
            QN[m, c] = Lmat[c, hkvecFinal]
        }
    }

    @Suppress("UNCHECKED_CAST")
    return Ret.pfqnAMVASchmidt(XN, QN, UN, CN, PN as List<MutableMap<Any, Double>>)
}

// Helper functions

private fun allLE(kvec: Matrix, nc: Matrix): Boolean {
    for (i in 0 until kvec.length()) {
        if (kvec[i] > nc[i]) {
            return false
        }
    }
    return true
}

private fun allGE(kvec: Matrix, value: Int): Boolean {
    for (j in 0 until kvec.length()) {
        if (kvec[j] < value) {
            return false
        }
    }
    return true
}

private fun prod(n: Matrix): Double {
    if (n.isEmpty) {
        return 1.0
    }
    var product = 1.0
    for (i in 0 until n.length()) {
        if (n[i].isInfinite()) {
            return Int.MAX_VALUE.toDouble()
        }
        product *= (n[i] + 1)
    }
    return product
}

private fun oner(kvec: Matrix, c: Int): Matrix {
    val result = Matrix(kvec)
    result[c] = result[c] - 1
    return result
}

private fun getBcn(D: Matrix, i: Int, c: Int, nvec: Matrix, C: Int, ns: Double): Double {
    var Bcn = D[i, c]
    if (nvec.elementSum() > 1.0) {
        val eps = 1e-12
        var sum = 0.0
        for (t in 0 until C) {
            sum += nvec[t] * D[i, t]
        }
        Bcn += (max(0.0, nvec.elementSum() - ns) / max(ns * (nvec.elementSum() - 1), eps) * (sum - D[i, c]))
    }
    return Bcn
}

private fun getBcnExt(u: Matrix, D: Matrix, i: Int, c: Int, nvec: Matrix, C: Int, ns: Double): Double {
    var weightedProb = 0.0
    val totalNonPinnedTime = u.sumRows(i) - u[i, C]

    for (s in 0 until C) {
        val prob = u[i, s] / totalNonPinnedTime
        weightedProb += (prob / D[i, s])
    }

    val meanInterdepartureTime = 1.0 / (ns * weightedProb)
    var Bcn = D[i, c]

    if (nvec.elementSum() > 1.0) {
        Bcn += (max(0.0, nvec.elementSum() - ns) * meanInterdepartureTime)
    }

    if (Bcn.isNaN() || Bcn.isInfinite()) {
        Bcn = 0.0
    }

    return Bcn
}
