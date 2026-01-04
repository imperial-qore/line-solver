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
import kotlin.math.floor
import kotlin.math.max
import kotlin.math.min
import kotlin.math.pow

/**
 * Akyildiz-Bolch AMVA method for multi-server BCMP networks.
 *
 * An approximate MVA algorithm using the Akyildiz-Bolch linearizer technique
 * for closed queueing networks with multi-server stations.
 *
 * @param serviceTimes Service time matrix (M x K)
 * @param N Population vector (1 x K)
 * @param V Visit ratio matrix (M x K)
 * @param nservers Number of servers at each station (M x 1)
 * @param schedStrategies Scheduling strategies for each station
 * @param fcfsSchmidt Whether to use Schmidt formula for FCFS stations
 * @param marginalProbMethod Method for marginal probability calculation ("ab" or "scat")
 * @return pfqnAMVAMS result containing Q, U, R, C, X, and iteration count
 */
fun ab_amva(
    serviceTimes: Matrix,
    N: Matrix,
    V: Matrix,
    nservers: Matrix,
    schedStrategies: List<SchedStrategy>,
    fcfsSchmidt: Boolean,
    marginalProbMethod: String
): Ret.pfqnAMVAMS {
    val M = serviceTimes.numRows
    val K = serviceTimes.numCols

    val ret = ab_linearizer(K, M, N, nservers, schedStrategies, V, serviceTimes, fcfsSchmidt, marginalProbMethod)

    return Ret.pfqnAMVAMS(ret.Q, ret.U, ret.R, ret.C, ret.X, ret.totiter)
}

/**
 * Akyildiz-Bolch linearizer method for multi-server BCMP networks.
 *
 * @param K Number of classes
 * @param M Number of stations
 * @param population Population vector
 * @param nservers Number of servers at each station
 * @param type Scheduling discipline at each station
 * @param v Visit matrix
 * @param s Service time matrix
 * @param fcfsSchmidt Whether to use Schmidt formula for FCFS
 * @param marginalProbMethod Method for marginal probability ("ab" or "scat")
 * @return pfqnAMVA result with performance measures
 */
fun ab_linearizer(
    K: Int,
    M: Int,
    population: Matrix,
    nservers: Matrix,
    type: List<SchedStrategy>,
    v: Matrix,
    s: Matrix,
    fcfsSchmidt: Boolean,
    marginalProbMethod: String
): Ret.pfqnAMVA {
    // Initialize queue length matrix
    val L = Matrix(M, K)
    for (i in 0 until M) {
        for (r in 0 until K) {
            L[i, r] = population[r] / M
        }
    }

    // Initialize L_without_r matrices
    val lWithoutR = Array(K) { Matrix(M, K) }

    for (i in 0 until M) {
        for (r in 0 until K) {
            for (t in 0 until K) {
                lWithoutR[r][i, t] = if (r == t) (population[r] - 1) / M else L[i, r]
            }
        }
    }

    // Fractional changes matrix (initialized to 0)
    val D = Array(M) { Array(K) { DoubleArray(K) } }

    // STEP 1: Apply core at full population
    var ret = pfqn_ab_core(K, M, population, nservers, type, v, s, 100, D, L, fcfsSchmidt, marginalProbMethod)
    val LUpdated = ret.Q

    // STEP 2: Apply core at N-e_k populations
    for (r in 0 until K) {
        val populationWithoutC = Matrix(population)
        populationWithoutC[r] = population[r] - 1

        val lWithoutC = Matrix(M, K)
        for (j in 0 until M) {
            for (c in 0 until K) {
                lWithoutC[j, c] = lWithoutR[c][j, c]
            }
        }

        val retWithoutC = pfqn_ab_core(K, M, populationWithoutC, nservers, type, v, s, 100, D, lWithoutC, fcfsSchmidt, marginalProbMethod)

        for (j in 0 until M) {
            for (c in 0 until K) {
                lWithoutR[c][j, r] = retWithoutC.Q[j, c]
            }
        }
    }

    // STEP 3: Compute estimates of F_mk(N) and F_mk(N-e_j)
    for (i in 0 until M) {
        for (r in 0 until K) {
            val F_ir = LUpdated[i, r] / population[r]
            for (t in 0 until K) {
                val divisor = if (r == t) population[r] - 1 else population[r]
                var F_irt = lWithoutR[r][i, t] / divisor
                if (F_irt.isNaN()) F_irt = 0.0
                D[i][r][t] = F_irt - F_ir
            }
        }
    }

    // STEP 4: Apply core at full population using L values from step 1 and D values from step 3
    return pfqn_ab_core(K, M, population, nservers, type, v, s, 100, D, LUpdated, fcfsSchmidt, marginalProbMethod)
}

/**
 * Akyildiz-Bolch core method for multi-server BCMP networks.
 *
 * @param K Number of classes
 * @param M Number of stations
 * @param population Population vector (1D matrix of size K)
 * @param nservers Number of servers at each station
 * @param type Scheduling discipline at each station
 * @param v Visit matrix
 * @param s Service time matrix
 * @param maxiter Maximum number of iterations
 * @param D Fractional changes matrix
 * @param lIn Initial queue length matrix
 * @param fcfsSchmidt Whether to use Schmidt formula for FCFS
 * @param marginalProbMethod Method for marginal probability
 * @return pfqnAMVA result with performance measures
 */
fun pfqn_ab_core(
    K: Int,
    M: Int,
    population: Matrix,
    nservers: Matrix,
    type: List<SchedStrategy>,
    v: Matrix,
    s: Matrix,
    maxiter: Int,
    D: Array<Array<DoubleArray>>,
    lIn: Matrix,
    fcfsSchmidt: Boolean,
    marginalProbMethod: String
): Ret.pfqnAMVA {
    var L = Matrix(lIn)
    val tol = 1 / (4000.0 + (16 * population.sumRows(0)))

    var totalIterations = 0
    val W = Matrix(M, K)

    while (totalIterations < maxiter) {
        val F = Matrix(M, K)
        val lWithoutJ = Array(M) { Array(K) { DoubleArray(K) } }

        for (i in 0 until M) {
            for (r in 0 until K) {
                if (population[r] > 0) {
                    F[i, r] = L[i, r] / population[r]
                } else {
                    F[i, r] = 0.0
                }
            }
        }

        for (i in 0 until M) {
            for (r in 0 until K) {
                for (t in 0 until K) {
                    val scalar = if (r == t) population[r] - 1 else population[r]
                    lWithoutJ[i][r][t] = scalar * (F[i, r] + D[i][r][t])
                }
            }
        }

        for (i in 0 until M) {
            for (r in 0 until K) {
                when {
                    type[i] == SchedStrategy.INF -> {
                        W[i, r] = s[i, r]
                    }
                    nservers[i] == 1.0 -> {
                        var totalQueueLengthKvecC = 0.0
                        for (c in 0 until K) {
                            totalQueueLengthKvecC += lWithoutJ[i][c][r]
                        }
                        W[i, r] = s[i, r] * (1 + totalQueueLengthKvecC)
                    }
                    fcfsSchmidt && type[i] == SchedStrategy.FCFS -> {
                        var waitTime = 0.0
                        var nvec = pprod(population)
                        while (nvec != null && nvec.value() >= 0) {
                            if (nvec[r] > 0) {
                                val Bcn = getBcnForAB(s, i, r, nvec, K, nservers[i])
                                val prob = getMarginalProb(oner(nvec, r), oner(population, r), population, lIn[i, r], K, i)
                                waitTime += (Bcn * prob)
                            }
                            nvec = pprod(nvec, population)
                        }
                        if (waitTime <= 1e-3) {
                            waitTime = 0.0
                        }
                        W[i, r] = waitTime
                    }
                    else -> {
                        var queueLength = 0.0
                        for (j in 0 until K) {
                            queueLength += lWithoutJ[i][j][r]
                        }
                        val numServers = nservers[i].toInt()
                        var multiServerStationWeightedQueueLength = 0.0
                        if (numServers > 1) {
                            val populationWithoutR = Matrix(population)
                            populationWithoutR[r] = populationWithoutR[r] - 1
                            val marginalProbs = findMarginalProbs(queueLength, numServers, populationWithoutR, r, marginalProbMethod)

                            for (j in 1 until numServers) {
                                multiServerStationWeightedQueueLength += (marginalProbs[j - 1] ?: 0.0) * (numServers - j)
                            }
                        }
                        val waitTime = (s[i, r] / numServers) * (1 + queueLength + multiServerStationWeightedQueueLength)
                        W[i, r] = waitTime
                    }
                }
            }
        }

        // Calculate cycle time for each class
        val C = Matrix.zeros(1, K)
        for (r in 0 until K) {
            var cycleTime = 0.0
            for (i in 0 until M) {
                cycleTime += v[i, r] * W[i, r]
            }
            C[0, r] = cycleTime
        }

        // Calculate queue length L_ir = N_r * W_ir / C_r
        val iterationQueueLength = Matrix(M, K)
        for (i in 0 until M) {
            for (r in 0 until K) {
                val queueLength = if (C[r] > 0) {
                    population[r] * (v[i, r] * W[i, r] / C[r])
                } else {
                    0.0
                }
                iterationQueueLength[i, r] = queueLength
            }
        }

        // Check convergence
        var maxDifference = 0.0
        for (i in 0 until M) {
            for (r in 0 until K) {
                val difference = kotlin.math.abs(L[i, r] - iterationQueueLength[i, r]) / population[r]
                maxDifference = max(maxDifference, difference)
            }
        }

        totalIterations++
        L = iterationQueueLength

        if (maxDifference < tol) {
            break
        }
    }

    // Calculate throughput
    val X = Matrix.zeros(1, K)
    for (r in 0 until K) {
        if (W[0, r] > 0) {
            X[0, r] = L[0, r] / W[0, r]
        } else {
            X[0, r] = 0.0
        }
    }

    // Calculate utilization
    val U = Matrix(M, K)
    for (i in 0 until M) {
        for (r in 0 until K) {
            if (s[i, r] > 0) {
                if (type[i] == SchedStrategy.INF) {
                    U[i, r] = X[r] * s[i, r]
                } else {
                    U[i, r] = (X[r] * s[i, r]) / nservers[i]
                }
            } else {
                U[i, r] = 0.0
            }
        }
    }

    val CResult = Matrix.zeros(1, K)

    return Ret.pfqnAMVA(L, U, W, CResult, CResult, X, totalIterations)
}

// Helper functions for AB algorithm

private fun oner(kvec: Matrix, c: Int): Matrix {
    val result = Matrix(kvec)
    result[c] = result[c] - 1
    return result
}

private fun getBcnForAB(D: Matrix, i: Int, c: Int, nvec: Matrix, K: Int, ns: Double): Double {
    var Bcn = D[i, c]
    if (nvec.elementSum() > 1.0) {
        val eps = 1e-12
        var sum = 0.0
        for (t in 0 until K) {
            sum += nvec[t] * D[i, t]
        }
        Bcn += (max(0.0, nvec.elementSum() - ns) / max(ns * (nvec.elementSum() - 1), eps) * (sum - D[i, c]))
    }
    return Bcn
}

private fun getMarginalProb(n: Matrix, k: Matrix, K: Matrix, L_jr: Double, R: Int, j: Int): Double {
    var prob = 1.0
    for (r in 0 until R) {
        val frac = L_jr / K[r]
        if (frac != 0.0) {
            val term1 = CombinatoricsUtils.binomialCoefficient(K[r].toInt(), n[r].toInt()).toDouble()
            val term2 = frac.pow(n[r])
            val term3 = (1 - frac).pow(K[r] - n[r])
            prob *= (term1 * term2 * term3)
        }
    }
    return prob
}

/**
 * Computes weight function for marginal probability calculation.
 */
private fun weightFun(population: Matrix, alpha: Double, beta: Double): Matrix {
    var maxClassPopulation = 0
    for (i in 0 until population.length()) {
        if (population[i] > maxClassPopulation) {
            maxClassPopulation = population[i].toInt()
        }
    }

    // Calculate scaling function PR
    val scalingFun = DoubleArray(maxClassPopulation + 1)
    if (maxClassPopulation >= 1) {
        scalingFun[1] = alpha
        for (n in 2..maxClassPopulation) {
            scalingFun[n] = beta * scalingFun[n - 1]
        }
    }

    // Calculate weight function W
    val weightFun = Array(maxClassPopulation + 1) { DoubleArray(maxClassPopulation + 1) }
    weightFun[0][0] = 1.0

    for (l in 1..maxClassPopulation) {
        for (j in 0 until l) {
            weightFun[l][j] = weightFun[l - 1][j] - (weightFun[l - 1][j] * scalingFun[l]) / 100.0
        }
        var sum = 0.0
        for (j in 0 until l) {
            sum += weightFun[l][j]
        }
        weightFun[l][l] = 1 - sum
    }

    return Matrix(weightFun)
}

/**
 * Finds marginal probabilities using specified method.
 *
 * @param avgJobs Average number of jobs
 * @param numServers Number of servers
 * @param population Population vector
 * @param classIdx Class index
 * @param marginalProbMethod Method ("scat" or "ab")
 * @return Map of job count to probability
 */
fun findMarginalProbs(
    avgJobs: Double,
    numServers: Int,
    population: Matrix,
    classIdx: Int,
    marginalProbMethod: String
): Map<Int, Double> {
    val marginalProbs = mutableMapOf<Int, Double>()

    if (marginalProbMethod == "scat") {
        val floorVal = floor(avgJobs).toInt()
        val ceilVal = floorVal + 1
        marginalProbs[floorVal] = ceilVal - avgJobs
        marginalProbs[ceilVal] = avgJobs - floorVal
        return marginalProbs
    }

    // AB method
    val ALPHA = 45.0
    val BETA = 0.7

    val w = weightFun(population, ALPHA, BETA)

    val floorVal = floor(avgJobs).toInt()
    val ceiling = floorVal + 1
    val maxVal = min((2 * floorVal) + 1, numServers - 2)

    for (j in 0..maxVal) {
        val prob: Double
        if (j <= floorVal) {
            val lDist = floorVal - j
            val lowerVal = floorVal - lDist
            val upperVal = ceiling + lDist
            prob = if (lDist > 25) {
                0.0
            } else {
                if (floorVal < population[classIdx]) {
                    w[floorVal, lDist] * ((upperVal - avgJobs) / (upperVal - lowerVal))
                } else {
                    0.0
                }
            }
            marginalProbs[j] = prob
        } else {
            val uDist = j - ceiling
            if (uDist > 25) {
                marginalProbs[j] = 0.0
            } else if (j > population[classIdx] - 1 && uDist < 25) {
                val existingProb = marginalProbs[(population[classIdx] - 1).toInt()] ?: 0.0
                val newProb = existingProb + (w[floorVal, uDist] - (marginalProbs[floorVal - uDist] ?: 0.0))
                marginalProbs[(population[classIdx] - 1).toInt()] = newProb
            } else {
                val newProb = w[floorVal, uDist] - (marginalProbs[floorVal - uDist] ?: 0.0)
                marginalProbs[j] = newProb
            }
        }
    }

    return marginalProbs
}
