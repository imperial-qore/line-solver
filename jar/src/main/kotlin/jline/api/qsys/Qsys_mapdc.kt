/**
 * @file MAP/D/c queueing system analysis
 *
 * Implements analysis of MAP/D/c queues using Q-MAM solver.
 * The queue has a Markovian Arrival Process (MAP) for arrivals,
 * deterministic service times, and c servers.
 *
 * Uses Q-MAM Q_CT_MAP_D_C which implements Non-Skip-Free Markov chain
 * analysis for this queue type.
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys

import jline.api.mc.ctmc_solve
import jline.lib.qmam.MAPDcOptions
import jline.lib.qmam.qCtMapDC
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Result of MAP/D/c queue analysis
 */
data class QsysMapDcResult(
    val meanQueueLength: Double,      // Mean number of customers in system
    val meanWaitingTime: Double,      // Mean waiting time in queue
    val meanSojournTime: Double,      // Mean sojourn time (waiting + service)
    val utilization: Double,          // Server utilization (per server)
    val queueLengthDist: Matrix,      // Queue length distribution P(Q=n)
    val waitingTimeDist: Matrix,      // Waiting time CDF at discrete points
    val analyzer: String              // Analyzer identifier
)

/**
 * Analyzes a MAP/D/c queue (multi-server with deterministic service).
 *
 * Uses Q-MAM qCtMapDC which implements Non-Skip-Free (NSF) Markov chain
 * analysis for MAP/D/c/FCFS queues.
 *
 * @param D0 MAP hidden transition matrix (n x n)
 * @param D1 MAP arrival transition matrix (n x n)
 * @param s Deterministic service time (positive scalar)
 * @param c Number of servers
 * @param maxNumComp Maximum number of queue length components (default 1000)
 * @param numSteps Number of waiting time distribution points per service interval (default 1)
 * @param verbose Verbosity level (default 0)
 * @return QsysMapDcResult with performance metrics
 */
@JvmOverloads
fun qsys_mapdc(
    D0: Matrix,
    D1: Matrix,
    s: Double,
    c: Int,
    maxNumComp: Int = 1000,
    numSteps: Int = 1,
    verbose: Int = 0
): QsysMapDcResult {
    // Call Q-MAM solver
    val result = qCtMapDC(D0, D1, s, c, MAPDcOptions(maxNumComp, verbose, numSteps))

    // Compute utilization
    val theta = ctmc_solve(D0.add(D1))
    val lambda = theta.mult(D1).elementSum()
    val rho = lambda * s / c

    // Compute mean queue length from distribution
    val ql = result.queueLength
    var meanQL = 0.0
    for (i in 0 until ql.numCols) {
        meanQL += i * ql[0, i]
    }

    // Compute mean waiting time from CDF using survival function integration
    // W(t) = CDF, so P(W > t) = 1 - W(t)
    // E[W] = integral of P(W > t) dt = sum of (1 - W(k*s)) * s (for basic numSteps=1 case)
    val w = result.waitingTime
    var meanWT = 0.0
    val stepSize = s / numSteps
    for (i in 0 until w.numCols - 1) {
        meanWT += (1.0 - w[0, i]) * stepSize
    }

    // Mean sojourn time = waiting + service
    val meanST = meanWT + s

    return QsysMapDcResult(
        meanQueueLength = meanQL,
        meanWaitingTime = meanWT,
        meanSojournTime = meanST,
        utilization = rho,
        queueLengthDist = ql,
        waitingTimeDist = w,
        analyzer = "Q-MAM:MAP/D/$c"
    )
}

/**
 * Simplified MAP/D/c analysis using MatrixCell input for arrival.
 *
 * @param arrival MAP arrival process as MatrixCell [D0, D1]
 * @param s Deterministic service time
 * @param c Number of servers
 * @return QsysMapDcResult with performance metrics
 */
fun qsys_mapdc(arrival: MatrixCell, s: Double, c: Int): QsysMapDcResult {
    require(arrival.size() >= 2) { "Arrival MAP must have at least 2 matrices [D0, D1]" }
    return qsys_mapdc(arrival[0], arrival[1], s, c)
}

/**
 * Analyzes a MAP/D/1 queue (single server convenience function).
 *
 * @param D0 MAP hidden transition matrix
 * @param D1 MAP arrival transition matrix
 * @param s Deterministic service time
 * @param maxNumComp Maximum number of queue length components (default 1000)
 * @param numSteps Number of waiting time distribution points per service interval (default 1)
 * @return QsysMapDcResult with performance metrics
 */
@JvmOverloads
fun qsys_mapd1(
    D0: Matrix,
    D1: Matrix,
    s: Double,
    maxNumComp: Int = 1000,
    numSteps: Int = 1
): QsysMapDcResult {
    return qsys_mapdc(D0, D1, s, 1, maxNumComp, numSteps)
}

/**
 * Simplified MAP/D/1 analysis using MatrixCell input for arrival.
 *
 * @param arrival MAP arrival process as MatrixCell [D0, D1]
 * @param s Deterministic service time
 * @return QsysMapDcResult with performance metrics
 */
fun qsys_mapd1(arrival: MatrixCell, s: Double): QsysMapDcResult {
    require(arrival.size() >= 2) { "Arrival MAP must have at least 2 matrices [D0, D1]" }
    return qsys_mapdc(arrival[0], arrival[1], s, 1)
}

/**
 * Analyzes a PH/D/c queue.
 *
 * Converts PH arrival process to equivalent MAP and uses Q-MAM.
 *
 * @param alpha Arrival PH initial probability vector
 * @param T Arrival PH generator matrix
 * @param s Deterministic service time
 * @param c Number of servers
 * @return QsysMapDcResult with performance metrics
 */
fun qsys_phdc(alpha: Matrix, T: Matrix, s: Double, c: Int): QsysMapDcResult {
    // Convert PH to MAP
    // D0 = T
    // D1 = (-T)*e * alpha
    val ones = Matrix.ones(T.numRows, 1)
    val exitRates = T.mult(ones).scale(-1.0)
    val D1 = exitRates.mult(alpha)

    return qsys_mapdc(T, D1, s, c)
}

/**
 * Queueing system mapdc algorithms
 */
@Suppress("unused")
class QsysMapdcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
