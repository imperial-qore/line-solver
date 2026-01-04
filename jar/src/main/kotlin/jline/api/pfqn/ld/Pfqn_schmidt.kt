package jline.api.pfqn.ld

import jline.GlobalConstants.Inf
import jline.io.Ret
import jline.lang.constant.SchedStrategy
import jline.util.matrix.Matrix
import jline.util.PopulationLattice
import jline.util.Utils.isInf
import kotlin.math.*
import org.apache.commons.math3.util.FastMath

/**
 * Schmidt method for load-dependent MVA with multi-server stations
 * 
 * Implementation of the Schmidt population recursion algorithm for analyzing closed 
 * queueing networks with load-dependent service and multi-server stations.
 * Supports INF, PS, and FCFS scheduling strategies with both single and multi-server configurations.
 *
 * Reference: R. Schmidt, "An approximate MVA algorithm for exponential, class-dependent
 * multiple server stations," Performance Evaluation, vol. 29, no. 4, pp. 245-254, 1997.
 *
 * @param D Service demand matrix (M x R) - service demands per station and class
 * @param N Population vector (R) - number of jobs per class
 * @param S Server count matrix (M x R) - number of servers per station and class  
 * @param sched Scheduling strategy list per station
 * @return pfqnSchmidt result containing performance measures and state probabilities
 */
fun pfqn_schmidt(D: Matrix, N: Matrix, S: Matrix, sched: List<SchedStrategy>): Ret.pfqnSchmidt {
    val startTime = System.nanoTime()
    
    val M = D.numRows // number of stations
    val R = D.numCols // number of classes
    
    // Initialize closed classes (all classes are closed)
    val closedClasses = IntArray(R) { it }
    val C = closedClasses.size
    val Nc = Matrix(N) // population vector for closed classes
    
    // Initialize performance measures
    val XN = Matrix.zeros(1, R) // System throughputs
    val UN = Matrix.zeros(M, R) // Utilizations  
    val CN = Matrix.zeros(M, R) // Waiting times
    val QN = Matrix.zeros(M, R) // Queue lengths
    val PN = mutableListOf<Matrix>() // State probabilities per station
    
    // Calculate products for fast hashing: prods[r] = prod(Nc[0:r-1]+1)
    val prods = Matrix.zeros(1, C)
    for (r in 0 until C) {
        var prod = 1.0
        for (i in 0 until r) {
            prod *= (Nc[i] + 1)
        }
        prods[0, r] = prod
    }
    
    // Calculate total number of states to enumerate
    var totalStates = 1.0
    for (r in 0 until R) {
        totalStates *= (Nc[r] + 1)
    }
    val numStates = totalStates.toInt()
    
    // Initialize intermediate data structures
    val L = mutableListOf<Matrix>() // Mean queue lengths
    val Pc = mutableListOf<Matrix?>() // State probabilities
    
    // Setup data structures based on scheduling strategy
    for (ist in 0 until M) {
        L.add(Matrix.zeros(R, numStates))
        
        when (sched[ist]) {
            SchedStrategy.INF -> {
                // Infinite server (delay station) - no state probabilities needed
                Pc.add(null)
                PN.add(Matrix.zeros(1, 1)) // Placeholder
            }
            SchedStrategy.PS -> {
                val isSingleServer = S[ist, 0] == 1.0
                if (isSingleServer) {
                    // Single server PS - no state probabilities needed
                    Pc.add(null)
                    PN.add(Matrix.zeros(1, 1)) // Placeholder
                } else {
                    // Multi-server PS - need scalar state probabilities Pr(j|N)
                    val maxServers = S[ist, 0].toInt()
                    Pc.add(Matrix.zeros(maxServers, numStates))
                    PN.add(Matrix.zeros(maxServers, numStates))
                }
            }
            SchedStrategy.FCFS -> {
                // Check if class-independent
                var classIndependent = true
                for (r in 1 until R) {
                    if (D[ist, r] != D[ist, 0]) {
                        classIndependent = false
                        break
                    }
                }
                
                val isSingleServer = S[ist, 0] == 1.0
                
                if (classIndependent) {
                    if (isSingleServer) {
                        // Single server class-independent FCFS
                        Pc.add(null)
                        PN.add(Matrix.zeros(1, 1)) // Placeholder
                    } else {
                        // Multi-server class-independent FCFS - scalar probabilities
                        val maxServers = S[ist, 0].toInt()
                        Pc.add(Matrix.zeros(maxServers, numStates))
                        PN.add(Matrix.zeros(maxServers, numStates))
                    }
                } else {
                    // Class-dependent FCFS - not product-form, need vector probabilities
                    Pc.add(Matrix.zeros(numStates, numStates))
                    PN.add(Matrix.zeros(numStates, numStates))
                }
            }
            else -> {
                // Default to single server behavior
                Pc.add(null)
                PN.add(Matrix.zeros(1, 1)) // Placeholder
            }
        }
    }
    
    // Initialize recursion variables
    val x = Matrix.zeros(C, numStates) // Throughputs
    val w = Array(M) { Matrix.zeros(C, numStates) } // Waiting times
    
    // Initialize base case: Pc(0|0) = 1 for all stations
    var kvec = PopulationLattice.pprod(Nc) // Start with zero vector
    for (ist in 0 until M) {
        Pc[ist]?.set(0, hashpop(kvec, Nc, prods), 1.0)
    }
    
    // Main population recursion loop
    kvec = PopulationLattice.pprod(kvec, Nc) // Move to next state
    
    while (allGE(kvec, Matrix.zeros(1, R)) && allLE(kvec, Nc)) {
        val hkvec = hashpop(kvec, Nc, prods)
        
        // For each station and class
        for (ist in 0 until M) {
            for (c in 0 until C) {
                if (kvec[c] == 0.0) continue
                
                val ns = S[ist, 0] // Number of servers
                val kvec_c = oner(kvec, c) // Decrement population of class c
                val hkvec_c = hashpop(kvec_c, Nc, prods)
                
                // Compute service rate based on scheduling strategy
                val serviceRate = when (sched[ist]) {
                    SchedStrategy.INF -> {
                        // Infinite server - service rate is always full
                        if (D[ist, c] > 0) 1.0 / D[ist, c] else 0.0
                    }
                    SchedStrategy.PS -> {
                        // Processor sharing
                        if (ns == 1.0) {
                            // Single server PS
                            if (D[ist, c] > 0) 1.0 / D[ist, c] else 0.0
                        } else {
                            // Multi-server PS
                            computeMultiServerPSRate(ist, c, kvec, ns.toInt(), D, Pc[ist])
                        }
                    }
                    SchedStrategy.FCFS -> {
                        // FCFS scheduling
                        if (ns == 1.0) {
                            // Single server FCFS
                            if (D[ist, c] > 0) 1.0 / D[ist, c] else 0.0
                        } else {
                            // Multi-server FCFS
                            computeMultiServerFCFSRate(ist, c, kvec, ns.toInt(), D, Pc[ist])
                        }
                    }
                    else -> {
                        // Default behavior
                        if (D[ist, c] > 0) 1.0 / D[ist, c] else 0.0
                    }
                }
                
                // Update waiting time using recursion
                w[ist][c, hkvec] = if (serviceRate > 0) {
                    (1.0 + L[ist][c, hkvec_c]) / serviceRate
                } else {
                    Inf
                }
            }
        }
        
        // Solve for throughputs using MVA equations
        for (c in 0 until C) {
            if (kvec[c] == 0.0) continue
            
            // Compute total response time for class c
            var totalResponseTime = 0.0
            for (ist in 0 until M) {
                totalResponseTime += w[ist][c, hkvec]
            }
            
            // Compute throughput for class c
            x[c, hkvec] = if (totalResponseTime > 0) {
                kvec[c] / totalResponseTime
            } else {
                0.0
            }
        }
        
        // Update queue lengths using Little's Law
        for (ist in 0 until M) {
            for (c in 0 until C) {
                L[ist][c, hkvec] = x[c, hkvec] * w[ist][c, hkvec]
            }
        }
        
        // Update state probabilities for multi-server stations
        for (ist in 0 until M) {
            if (Pc[ist] != null) {
                updateStateProbabilities(ist, kvec, hkvec, sched[ist], S[ist, 0].toInt(), Pc[ist]!!, L[ist])
            }
        }
        
        // Move to next state
        kvec = PopulationLattice.pprod(kvec, Nc)
    }
    
    // Extract final performance measures
    val finalState = hashpop(Nc, Nc, prods)
    
    // System throughputs
    for (c in 0 until C) {
        XN[0, c] = x[c, finalState]
    }
    
    // Response times and queue lengths
    for (ist in 0 until M) {
        for (c in 0 until C) {
            CN[ist, c] = w[ist][c, finalState]
            QN[ist, c] = L[ist][c, finalState]
        }
    }
    
    // Utilizations
    for (ist in 0 until M) {
        for (c in 0 until C) {
            UN[ist, c] = XN[0, c] * D[ist, c]
        }
    }
    
    // Response times (from queue lengths)
    val RN = Matrix.zeros(M, R)
    for (ist in 0 until M) {
        for (c in 0 until C) {
            RN[ist, c] = if (XN[0, c] > 0) QN[ist, c] / XN[0, c] else 0.0
        }
    }
    
    // Throughputs (per station)
    val TN = Matrix.zeros(M, R)
    for (ist in 0 until M) {
        for (c in 0 until C) {
            TN[ist, c] = XN[0, c]
        }
    }
    
    // Copy final state probabilities
    for (ist in 0 until M) {
        if (Pc[ist] != null) {
            PN[ist] = Matrix(Pc[ist]!!)
        }
    }
    
    val runtime = (System.nanoTime() - startTime) / 1e9
    
    return Ret.pfqnSchmidt(
        QN, UN, RN, TN, CN, XN, PN, "schmidt", numStates, runtime
    )
}

/**
 * Hash function for population state vectors
 */
private fun hashpop(kvec: Matrix, Nc: Matrix, prods: Matrix): Int {
    var hash = 0.0
    for (r in 0 until kvec.length()) {
        hash += kvec[r] * prods[0, r]
    }
    return hash.toInt()
}

/**
 * Decrement population vector for specified class
 */
private fun oner(kvec: Matrix, c: Int): Matrix {
    val result = Matrix(kvec)
    if (result[c] > 0) {
        result[c] = result[c] - 1
    }
    return result
}

/**
 * Check if all elements are greater than or equal to comparison
 */
private fun allGE(a: Matrix, b: Matrix): Boolean {
    for (i in 0 until a.length()) {
        if (a[i] < b[i]) return false
    }
    return true
}

/**
 * Check if all elements are less than or equal to comparison  
 */
private fun allLE(a: Matrix, b: Matrix): Boolean {
    for (i in 0 until a.length()) {
        if (a[i] > b[i]) return false
    }
    return true
}

/**
 * Compute service rate for multi-server PS station
 */
private fun computeMultiServerPSRate(
    stationIdx: Int,
    classIdx: Int, 
    kvec: Matrix,
    numServers: Int,
    D: Matrix,
    Pc: Matrix?
): Double {
    if (D[stationIdx, classIdx] == 0.0) return 0.0
    
    // For multi-server PS, service rate depends on number of active servers
    val totalJobs = kvec.elementSum().toInt()
    val activeServers = min(totalJobs, numServers)
    
    return if (activeServers > 0) {
        activeServers.toDouble() / D[stationIdx, classIdx]
    } else {
        1.0 / D[stationIdx, classIdx]
    }
}

/**
 * Compute service rate for multi-server FCFS station
 */
private fun computeMultiServerFCFSRate(
    stationIdx: Int,
    classIdx: Int,
    kvec: Matrix, 
    numServers: Int,
    D: Matrix,
    Pc: Matrix?
): Double {
    if (D[stationIdx, classIdx] == 0.0) return 0.0
    
    // For multi-server FCFS, service rate depends on number of busy servers
    val totalJobs = kvec.elementSum().toInt()
    val busyServers = min(totalJobs, numServers)
    
    return if (busyServers > 0) {
        busyServers.toDouble() / D[stationIdx, classIdx]
    } else {
        1.0 / D[stationIdx, classIdx]
    }
}

/**
 * Update state probabilities for multi-server stations
 */
private fun updateStateProbabilities(
    stationIdx: Int,
    kvec: Matrix,
    hkvec: Int,
    sched: SchedStrategy,
    numServers: Int,
    Pc: Matrix,
    L: Matrix
) {
    when (sched) {
        SchedStrategy.PS -> {
            if (numServers > 1) {
                // Multi-server PS - update scalar probabilities
                val totalJobs = kvec.elementSum().toInt()
                for (j in 0 until min(totalJobs + 1, Pc.numRows)) {
                    // Simplified probability update - could be enhanced with exact calculations
                    val meanJobs = L.sumRows()[hkvec]
                    if (meanJobs > 0) {
                        Pc[j, hkvec] = exp(-meanJobs) * FastMath.pow(meanJobs, j.toDouble()) / factorial(j)
                    }
                }
            }
        }
        SchedStrategy.FCFS -> {
            if (numServers > 1) {
                // Multi-server FCFS - update probabilities
                val totalJobs = kvec.elementSum().toInt()
                for (j in 0 until min(totalJobs + 1, Pc.numRows)) {
                    // Simplified probability update 
                    val meanJobs = L.sumRows()[hkvec]
                    if (meanJobs > 0) {
                        Pc[j, hkvec] = exp(-meanJobs) * FastMath.pow(meanJobs, j.toDouble()) / factorial(j)
                    }
                }
            }
        }
        else -> {
            // No probability update needed
        }
    }
}

/**
 * Factorial function for probability calculations
 */
private fun factorial(n: Int): Double {
    if (n <= 1) return 1.0
    var result = 1.0
    for (i in 2..n) {
        result *= i.toDouble()
    }
    return result
}
/**
 * PFQN schmidt algorithms
 */
@Suppress("unused")
class PfqnSchmidtAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}