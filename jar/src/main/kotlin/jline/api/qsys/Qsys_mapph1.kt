/**
 * @file MAP/PH/1 queueing system analysis
 *
 * Implements analysis of MAP/PH/1 queues using BUTools MMAPPH1FCFS solver.
 * The queue has a Markovian Arrival Process (MAP) for arrivals and a
 * Phase-Type (PH) distribution for service times.
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys

import jline.api.mc.ctmc_solve
import jline.lib.butools.MMAPPH1FCFS
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Result of MAP/PH type queue analysis.
 *
 * @property meanQueueLength Mean number of customers in the system
 * @property meanWaitingTime Mean waiting time in queue (excluding service)
 * @property meanSojournTime Mean sojourn time (waiting + service)
 * @property utilization Server utilization
 * @property queueLengthDist Queue length distribution P(Q=n) as row vector
 * @property queueLengthMoments Raw moments of queue length [E[Q], E[Q^2], ...]
 * @property sojournTimeMoments Raw moments of sojourn time [E[T], E[T^2], ...]
 * @property analyzer Name of the analyzer used
 */
data class QsysMapPhResult(
    val meanQueueLength: Double,
    val meanWaitingTime: Double,
    val meanSojournTime: Double,
    val utilization: Double,
    val queueLengthDist: Matrix?,
    val queueLengthMoments: Matrix?,
    val sojournTimeMoments: Matrix?,
    val analyzer: String
)

/**
 * Analyzes a MAP/PH/1 queue.
 *
 * Uses BUTools MMAPPH1FCFS which computes performance measures for a
 * single-server FCFS queue with Markovian Arrival Process arrivals
 * and Phase-Type service times.
 *
 * @param D0 MAP hidden transition matrix (n x n)
 * @param D1 MAP arrival transition matrix (n x n)
 * @param sigma PH service initial probability vector (1 x m)
 * @param S PH service generator matrix (m x m)
 * @param numQLMoms Number of queue length moments to compute (default 3)
 * @param numQLProbs Number of queue length probabilities to compute (default 100)
 * @param numSTMoms Number of sojourn time moments to compute (default 3)
 * @return QsysMapPhResult with performance metrics
 */
@JvmOverloads
fun qsys_mapph1(
    D0: Matrix,
    D1: Matrix,
    sigma: Matrix,
    S: Matrix,
    numQLMoms: Int = 3,
    numQLProbs: Int = 100,
    numSTMoms: Int = 3
): QsysMapPhResult {
    // Ensure sigma is a row vector (1 x n)
    // If passed as column vector (n x 1), transpose it
    val sigmaRow = if (sigma.numRows > 1 && sigma.numCols == 1) {
        sigma.transpose()
    } else {
        sigma
    }

    // Build MMAP structure for BUTools (single class)
    val D = MatrixCell(2)
    D[0] = D0
    D[1] = D1

    // Service parameters as maps (single class indexed by 0)
    val sigmaMap = HashMap<Int?, Matrix>()
    sigmaMap[0] = sigmaRow
    val sMap = HashMap<Int?, Matrix>()
    sMap[0] = S

    // Call BUTools solver
    val result = MMAPPH1FCFS(
        D, sigmaMap, sMap,
        numQLMoms, numQLProbs, numSTMoms,
        null, false, false, null, null
    )

    // Compute utilization from arrival and service rates
    val theta = ctmc_solve(D0.add(D1))
    val lambda = theta.mult(D1).elementSum()

    val negSinv = S.scale(-1.0).inv()
    val ones = Matrix.ones(S.numRows, 1)
    val meanService = sigmaRow.mult(negSinv).mult(ones)[0, 0]
    val mu = 1.0 / meanService
    val rho = lambda / mu

    // Extract queue length moments
    val ncMoms = result["ncMoms"]?.get(0)
    val meanQL: Double = ncMoms?.get(0, 0) ?: 0.0

    // Extract sojourn time moments
    val stMoms = result["stNoms"]?.get(0)
    val meanST: Double = stMoms?.get(0, 0) ?: 0.0

    // Waiting time = sojourn time - service time
    val meanWT = maxOf(0.0, meanST - meanService)

    // Extract queue length distribution
    val ncDistr = result["ncDistr"]?.get(0)

    return QsysMapPhResult(
        meanQueueLength = meanQL,
        meanWaitingTime = meanWT,
        meanSojournTime = meanST,
        utilization = rho,
        queueLengthDist = ncDistr,
        queueLengthMoments = ncMoms,
        sojournTimeMoments = stMoms,
        analyzer = "BUTools:MMAPPH1FCFS"
    )
}

/**
 * Simplified MAP/PH/1 analysis using MatrixCell inputs.
 *
 * @param arrival MAP arrival process as MatrixCell [D0, D1]
 * @param service PH service process as MatrixCell [sigma, S] or just [S] with uniform sigma
 * @return QsysMapPhResult with performance metrics
 */
fun qsys_mapph1(arrival: MatrixCell, service: MatrixCell): QsysMapPhResult {
    val D0 = arrival[0]
    val D1 = arrival[1]

    // Extract service PH parameters
    val (sigma, S) = if (service.size() >= 2) {
        Pair(service[0], service[1])
    } else {
        // Only S provided, use uniform initial distribution
        val S = service[0]
        val n = S.numRows
        Pair(Matrix.ones(1, n).scale(1.0 / n), S)
    }

    return qsys_mapph1(D0, D1, sigma, S)
}

/**
 * Queueing system mapph1 algorithms
 */
@Suppress("unused")
class QsysMapph1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}
