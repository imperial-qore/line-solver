/**
 * @file Akyildiz-Bolch linearizer for load-dependent multi-server queueing networks
 * 
 * Implements the Akyildiz-Bolch linearizer method for analyzing closed product-form queueing networks
 * with load-dependent service rates and multi-server processor sharing stations. Uses iterative
 * marginal probability approximations to handle complex state dependencies in multi-server environments.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.io.Ret
import jline.lang.constant.SchedStrategy
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.CombinatoricsUtils
import kotlin.math.*

/**
 * Akyildiz-Bolch (A/B) linearizer method
 *
 * Reference: A. Akyildiz and G. Bolch, "Mean Value Analysis Approximation for Multiple Server
 * Queueing Networks," Performance Evaluation, vol. 8, no. 2, pp. 77-91, 1988.
 *
 * @param D Service demand matrix (M x R) - service demands per station and class
 * @param N Population vector (R) - number of jobs per class  
 * @param S Server count matrix (M x R) - number of servers per station and class
 * @param sched Scheduling strategy list per station
 * @return pfqnAB result containing performance measures
 */
fun pfqn_ab(D: Matrix, N: Matrix, S: Matrix, sched: List<SchedStrategy>): Ret.pfqnAB {
    val startTime = System.nanoTime()
    
    val M = D.numRows // number of stations
    val R = D.numCols // number of classes
    
    // Initialize performance measures
    var QN = Matrix.zeros(M, R) // Queue lengths
    var UN = Matrix.zeros(M, R) // Utilizations
    var RN = Matrix.zeros(M, R) // Response times
    var TN = Matrix.zeros(M, R) // Throughputs  
    var CN = Matrix.zeros(M, R) // Waiting times
    var XN = Matrix.zeros(1, R) // System throughputs
    
    // Check if network is suitable for A/B method
    var hasMultiServerPS = false
    for (ist in 0 until M) {
        if (sched[ist] == SchedStrategy.PS) {
            for (r in 0 until R) {
                if (S[ist, r] > 1) {
                    hasMultiServerPS = true
                    break
                }
            }
        }
    }
    
    if (!hasMultiServerPS) {
        // Fall back to standard MVA for non-multi-server PS networks
        return fallbackToStandardMVA(D, N, S, sched, startTime)
    }
    
    // A/B algorithm parameters
    val maxIter = 1000
    val tolerance = 1e-6
    var iter = 0
    
    // Initialize mean queue lengths (starting guess)
    var Q_prev = Matrix.zeros(M, R)
    for (ist in 0 until M) {
        for (r in 0 until R) {
            Q_prev[ist, r] = N[r] / M.toDouble() // Uniform distribution
        }
    }
    
    // Iterative A/B linearizer algorithm
    while (iter < maxIter) {
        iter++
        
        // Step 1: Compute marginal probabilities for each station
        val marginalProbs = Array(M) { Array(R) { HashMap<Int, Double>() } }
        
        for (ist in 0 until M) {
            if (sched[ist] == SchedStrategy.PS) {
                for (r in 0 until R) {
                    if (S[ist, r] > 1) {
                        // Use A/B method for multi-server PS stations
                        marginalProbs[ist][r] = findMarginalProbs(
                            Q_prev[ist, r], 
                            S[ist, r].toInt(), 
                            N, 
                            r, 
                            "ab"
                        )
                    } else {
                        // Single server PS - use scatter method
                        marginalProbs[ist][r] = findMarginalProbs(
                            Q_prev[ist, r], 
                            1, 
                            N, 
                            r, 
                            "scat"
                        )
                    }
                }
            }
        }
        
        // Step 2: Solve reduced MVA problems for each class
        for (r in 0 until R) {
            // Create reduced population vector (N-1 for class r)
            val N_reduced = Matrix(N)
            N_reduced[r] = N_reduced[r] - 1
            
            if (N_reduced[r] < 0) continue
            
            // Compute waiting times using marginal probabilities
            for (ist in 0 until M) {
                when (sched[ist]) {
                    SchedStrategy.PS -> {
                        if (S[ist, r] > 1) {
                            // Multi-server PS station - use A/B approximation
                            CN[ist, r] = computeMultiServerPSWaitingTime(
                                ist, r, D[ist, r], S[ist, r].toInt(), marginalProbs[ist][r]
                            )
                        } else {
                            // Single server PS station
                            CN[ist, r] = D[ist, r] * (1 + Q_prev[ist, r])
                        }
                    }
                    SchedStrategy.FCFS -> {
                        // FCFS station (load-independent)
                        CN[ist, r] = D[ist, r] * (1 + Q_prev[ist, r])
                    }
                    SchedStrategy.INF -> {
                        // Infinite server (delay station)
                        CN[ist, r] = D[ist, r]
                    }
                    else -> {
                        // Default to FCFS behavior
                        CN[ist, r] = D[ist, r] * (1 + Q_prev[ist, r])
                    }
                }
            }
            
            // Compute system response time and throughput for class r
            val systemResponseTime = CN.sumCols()[r]
            XN[r] = N[r] / systemResponseTime
        }
        
        // Step 3: Update queue lengths using Little's Law
        val Q_new = Matrix.zeros(M, R)
        for (ist in 0 until M) {
            for (r in 0 until R) {
                Q_new[ist, r] = XN[r] * CN[ist, r]
            }
        }
        
        // Step 4: Check convergence
        var maxChange = 0.0
        for (ist in 0 until M) {
            for (r in 0 until R) {
                val change = abs(Q_new[ist, r] - Q_prev[ist, r])
                maxChange = max(maxChange, change)
            }
        }
        
        Q_prev = Q_new
        
        if (maxChange < tolerance) {
            break
        }
    }
    
    // Final performance measures
    QN = Q_prev
    
    // Compute utilizations
    for (ist in 0 until M) {
        for (r in 0 until R) {
            when (sched[ist]) {
                SchedStrategy.PS -> {
                    if (S[ist, r] > 1) {
                        // Multi-server utilization
                        UN[ist, r] = XN[r] * D[ist, r] / S[ist, r]
                    } else {
                        // Single server utilization  
                        UN[ist, r] = XN[r] * D[ist, r]
                    }
                }
                else -> {
                    UN[ist, r] = XN[r] * D[ist, r]
                }
            }
        }
    }
    
    // Compute response times
    for (ist in 0 until M) {
        for (r in 0 until R) {
            RN[ist, r] = if (XN[r] > 0) QN[ist, r] / XN[r] else 0.0
        }
    }
    
    // Throughputs (per station)
    for (ist in 0 until M) {
        for (r in 0 until R) {
            TN[ist, r] = XN[r]
        }
    }
    
    val runtime = (System.nanoTime() - startTime) / 1e9
    
    return Ret.pfqnAB(
        QN, UN, RN, TN, CN, XN, "ab", iter, runtime
    )
}

/**
 * Weight function calculation for A/B method
 */
private fun weightFun(population: Matrix, alpha: Double, beta: Double): Matrix {
    var maxClassPopulation = 0
    for (i in 0 until population.length()) {
        if (population[i] > maxClassPopulation) {
            maxClassPopulation = population[i].toInt()
        }
    }
    
    // Calculate scaling function PR for all n = 1, ..., max(K_r)
    val scalingFun = DoubleArray(maxClassPopulation + 1)
    if (maxClassPopulation >= 1) {
        scalingFun[1] = alpha // PR[1] = α
        for (n in 2..maxClassPopulation) {
            scalingFun[n] = beta * scalingFun[n - 1] // PR[n] = β * PR[n-1]
        }
    }
    
    // Calculate weight function W for all values l = 1, ..., max(K_r)
    val weightFun = Array(maxClassPopulation + 1) { DoubleArray(maxClassPopulation + 1) }
    
    // Initialize W[0, 0] = 1
    weightFun[0][0] = 1.0
    
    for (l in 1..maxClassPopulation) {
        // Calculate W[l, j] for j = 0, ..., (l-1)
        for (j in 0 until l) {
            weightFun[l][j] = weightFun[l - 1][j] - (weightFun[l - 1][j] * scalingFun[l]) / 100.0
        }
        
        // Calculate W[l, l]
        var sum = 0.0
        for (j in 0 until l) {
            sum += weightFun[l][j]
        }
        weightFun[l][l] = 1.0 - sum
    }
    
    return Matrix(weightFun)
}

/**
 * Find marginal probabilities using A/B or scatter method
 */
private fun findMarginalProbs(
    avgJobs: Double, 
    numServers: Int, 
    population: Matrix, 
    classIdx: Int, 
    method: String
): HashMap<Int, Double> {
    val marginalProbs = HashMap<Int, Double>()
    
    if (method.lowercase() == "scat") {
        // Scatter method for single server
        val floor = floor(avgJobs).toInt()
        val ceiling = ceil(avgJobs).toInt()
        marginalProbs[floor] = ceiling - avgJobs
        marginalProbs[ceiling] = avgJobs - floor
        return marginalProbs
    }
    
    // A/B method parameters
    val ALPHA = 45.0
    val BETA = 0.7
    
    val w = weightFun(population, ALPHA, BETA)
    
    val floor = floor(avgJobs).toInt()
    val ceiling = floor + 1
    val maxVal = min((2 * floor) + 1, numServers - 2)
    
    for (j in 0..maxVal) {
        val prob = if (j <= floor) {
            val lDist = floor - j
            val lowerVal = floor - lDist
            val upperVal = ceiling + lDist
            if (lDist > 25) {
                0.0
            } else {
                if (floor < population[classIdx]) {
                    w[floor, lDist] * ((upperVal - avgJobs) / (upperVal - lowerVal))
                } else {
                    0.0
                }
            }
        } else {
            val uDist = j - ceiling
            val lowerVal = floor - uDist
            val upperVal = ceiling + uDist
            if (uDist > 25) {
                0.0
            } else {
                if (ceiling < population[classIdx]) {
                    w[ceiling, uDist] * ((avgJobs - lowerVal) / (upperVal - lowerVal))
                } else {
                    0.0
                }
            }
        }
        marginalProbs[j] = prob
    }
    
    return marginalProbs
}

/**
 * Compute waiting time for multi-server PS station using A/B approximation
 */
private fun computeMultiServerPSWaitingTime(
    stationIdx: Int,
    classIdx: Int, 
    demand: Double,
    numServers: Int,
    marginalProbs: HashMap<Int, Double>
): Double {
    var expectedWaitingTime = 0.0
    
    for ((jobs, prob) in marginalProbs) {
        val serviceRate = min(jobs.toDouble(), numServers.toDouble())
        if (serviceRate > 0) {
            expectedWaitingTime += prob * demand * (1 + jobs / serviceRate)
        } else {
            expectedWaitingTime += prob * demand
        }
    }
    
    return expectedWaitingTime
}

/**
 * Fallback to standard MVA for non-multi-server networks
 */
private fun fallbackToStandardMVA(
    D: Matrix, 
    N: Matrix, 
    S: Matrix, 
    sched: List<SchedStrategy>, 
    startTime: Long
): Ret.pfqnAB {
    // Simple single-server MVA implementation
    val M = D.numRows
    val R = D.numCols
    
    val QN = Matrix.zeros(M, R)
    val UN = Matrix.zeros(M, R) 
    val RN = Matrix.zeros(M, R)
    val TN = Matrix.zeros(M, R)
    val CN = Matrix.zeros(M, R)
    val XN = Matrix.zeros(1, R)
    
    // Basic MVA for comparison
    for (r in 0 until R) {
        val totalDemand = D.sumCols()[r]
        XN[r] = N[r] / totalDemand
        
        for (ist in 0 until M) {
            CN[ist, r] = D[ist, r]
            QN[ist, r] = XN[r] * CN[ist, r] 
            UN[ist, r] = XN[r] * D[ist, r]
            RN[ist, r] = CN[ist, r]
            TN[ist, r] = XN[r]
        }
    }
    
    val runtime = (System.nanoTime() - startTime) / 1e9
    
    return Ret.pfqnAB(
        QN, UN, RN, TN, CN, XN, "ab-fallback", 1, runtime
    )
}
/**
 * PFQN ab algorithms
 */
@Suppress("unused")
class PfqnAbAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}