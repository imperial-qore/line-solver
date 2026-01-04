/**
 * @file MAP/M/c queueing system analysis
 *
 * Implements analysis of MAP/M/c queues using Q-MAM solver.
 * The queue has a Markovian Arrival Process (MAP) for arrivals,
 * exponential service times, and c servers.
 *
 * Uses Q-MAM since BUTools MMAPPH1FCFS only supports single-server queues.
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys

import jline.api.mc.ctmc_solve
import jline.lib.qmam.MAPMcOptions
import jline.lib.qmam.qCtMapMC
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Analyzes a MAP/M/c queue (multi-server with exponential service).
 *
 * Uses Q-MAM qCtMapMC which implements the Gaver-Jacobs-Latouche
 * level-dependent QBD approach for MAP/M/c/FCFS queues.
 *
 * @param D0 MAP hidden transition matrix (n x n)
 * @param D1 MAP arrival transition matrix (n x n)
 * @param mu Exponential service rate
 * @param c Number of servers
 * @param mode Solver mode: "SylvesCR" (default), "DirectCR", "SylvesFI", "DirectFI"
 * @param maxNumComp Maximum number of queue length components (default 1000)
 * @param verbose Verbosity level (default 0)
 * @return QsysMapPhResult with performance metrics
 */
@JvmOverloads
fun qsys_mapmc(
    D0: Matrix,
    D1: Matrix,
    mu: Double,
    c: Int,
    mode: String = "SylvesCR",
    maxNumComp: Int = 1000,
    verbose: Int = 0
): QsysMapPhResult {
    // Call Q-MAM solver
    val result = qCtMapMC(D0, D1, mu, c, MAPMcOptions(mode, maxNumComp, verbose))

    // Compute utilization
    val theta = ctmc_solve(D0.add(D1))
    val lambda = theta.mult(D1).elementSum()
    val rho = lambda / (mu * c)

    // Compute mean queue length from distribution
    val ql = result.queueLength
    var meanQL = 0.0
    for (i in 0 until ql.numCols) {
        meanQL += i * ql[0, i]
    }

    // Compute mean waiting time from PH representation
    var meanWT = 0.0
    if (result.waitAlpha != null && result.Smat != null) {
        val negSinv = result.Smat.scale(-1.0).inv()
        val ones = Matrix.ones(result.Smat.numRows, 1)
        meanWT = result.waitAlpha.mult(negSinv).mult(ones)[0, 0]
    }

    // Mean service time
    val meanService = 1.0 / mu

    // Mean sojourn time = waiting + service
    val meanST = meanWT + meanService

    return QsysMapPhResult(
        meanQueueLength = meanQL,
        meanWaitingTime = meanWT,
        meanSojournTime = meanST,
        utilization = rho,
        queueLengthDist = ql,
        queueLengthMoments = null,  // Not directly computed by Q-MAM
        sojournTimeMoments = null,
        analyzer = "Q-MAM:MAP/M/$c"
    )
}

/**
 * Simplified MAP/M/c analysis using MatrixCell input for arrival.
 *
 * @param arrival MAP arrival process as MatrixCell [D0, D1]
 * @param mu Exponential service rate
 * @param c Number of servers
 * @return QsysMapPhResult with performance metrics
 */
fun qsys_mapmc(arrival: MatrixCell, mu: Double, c: Int): QsysMapPhResult {
    require(arrival.size() >= 2) { "Arrival MAP must have at least 2 matrices [D0, D1]" }
    return qsys_mapmc(arrival[0], arrival[1], mu, c)
}

/**
 * Analyzes a MAP/M/1 queue (single server convenience function).
 *
 * For single-server queues, this delegates to BUTools via qsys_mapph1
 * for better performance metrics including sojourn time moments.
 *
 * @param D0 MAP hidden transition matrix
 * @param D1 MAP arrival transition matrix
 * @param mu Exponential service rate
 * @return QsysMapPhResult with performance metrics
 */
fun qsys_mapm1(D0: Matrix, D1: Matrix, mu: Double): QsysMapPhResult {
    // Create exponential service PH
    val sigma = Matrix(1, 1)
    sigma[0, 0] = 1.0
    val S = Matrix(1, 1)
    S[0, 0] = -mu

    // Use BUTools for single-server case
    return qsys_mapph1(D0, D1, sigma, S)
}

/**
 * Analyzes a PH/M/c queue.
 *
 * Converts PH arrival process to equivalent MAP and uses Q-MAM.
 *
 * @param alpha Arrival PH initial probability vector
 * @param T Arrival PH generator matrix
 * @param mu Exponential service rate
 * @param c Number of servers
 * @return QsysMapPhResult with performance metrics
 */
fun qsys_phmc(alpha: Matrix, T: Matrix, mu: Double, c: Int): QsysMapPhResult {
    // Convert PH to MAP
    // D0 = T
    // D1 = (-T)*e * alpha
    val ones = Matrix.ones(T.numRows, 1)
    val exitRates = T.mult(ones).scale(-1.0)
    val D1 = exitRates.mult(alpha)

    return qsys_mapmc(T, D1, mu, c)
}

/**
 * Queueing system mapmc algorithms
 */
@Suppress("unused")
class QsysMapmcAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
