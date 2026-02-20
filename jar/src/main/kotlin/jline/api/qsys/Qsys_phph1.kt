/**
 * @file PH/PH/1 queueing system analysis
 *
 * Implements analysis of PH/PH/1 queues using BUTools MMAPPH1FCFS solver.
 * The queue has Phase-Type distributions for both arrivals and service.
 * The arrival PH is converted to equivalent MAP representation for the solver.
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys

import jline.lib.butools.MMAPPH1FCFS
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Analyzes a PH/PH/1 queue.
 *
 * Uses BUTools MMAPPH1FCFS by converting the arrival PH to an equivalent
 * Markovian Arrival Process. The PH renewal process is embedded as a MAP
 * where arrivals occur at PH absorption epochs.
 *
 * @param alpha Arrival PH initial probability vector (1 x n)
 * @param T Arrival PH generator matrix (n x n)
 * @param beta Service PH initial probability vector (1 x m)
 * @param S Service PH generator matrix (m x m)
 * @param numQLMoms Number of queue length moments to compute (default 3)
 * @param numQLProbs Number of queue length probabilities to compute (default 100)
 * @param numSTMoms Number of sojourn time moments to compute (default 3)
 * @return QsysMapPhResult with performance metrics
 */
@JvmOverloads
fun qsys_phph1(
    alpha: Matrix,
    T: Matrix,
    beta: Matrix,
    S: Matrix,
    numQLMoms: Int = 3,
    numQLProbs: Int = 100,
    numSTMoms: Int = 3
): QsysMapPhResult {
    // Convert arrival PH to MAP representation
    // For a PH renewal process:
    // D0 = T (hidden transitions within the PH)
    // D1 = (-T)*e * alpha (absorption followed by restart)
    val (D0, D1) = phToMap(alpha, T)

    // Build arrival MMAP structure for BUTools (single class)
    val D = MatrixCell(2)
    D[0] = D0
    D[1] = D1

    // Service parameters as maps (single class indexed by 0)
    val sigmaMap = HashMap<Int?, Matrix>()
    sigmaMap[0] = beta
    val sMap = HashMap<Int?, Matrix>()
    sMap[0] = S

    // Call BUTools solver
    val result = MMAPPH1FCFS(
        D, sigmaMap, sMap,
        numQLMoms, numQLProbs, numSTMoms,
        null, false, false, null, null
    )

    // Compute rates
    val negTinv = T.scale(-1.0).inv()
    val ones = Matrix.ones(T.numRows, 1)
    val meanInterarrival = alpha.mult(negTinv).mult(ones)[0, 0]
    val lambda = 1.0 / meanInterarrival

    val negSinv = S.scale(-1.0).inv()
    val onesS = Matrix.ones(S.numRows, 1)
    val meanService = beta.mult(negSinv).mult(onesS)[0, 0]
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
 * Converts a PH distribution to its equivalent MAP representation.
 *
 * For a PH renewal process, the MAP has:
 * - D0 = T (transitions within the PH, no arrival)
 * - D1 = t * alpha where t = -T*e (exit rates times restart distribution)
 *
 * @param alpha PH initial probability vector (1 x n)
 * @param T PH generator matrix (n x n)
 * @return Pair of (D0, D1) matrices representing the MAP
 */
private fun phToMap(alpha: Matrix, T: Matrix): Pair<Matrix, Matrix> {
    // D0 = T (hidden transitions)
    val D0 = T.copy()

    // t = -T * e (exit rate vector, column)
    val ones = Matrix.ones(T.numRows, 1)
    val exitRates = T.mult(ones).scale(-1.0)

    // D1 = t * alpha (restart to initial distribution)
    val D1 = exitRates.mult(alpha)

    return Pair(D0, D1)
}

/**
 * Simplified PH/PH/1 analysis using MatrixCell inputs.
 *
 * @param arrival PH arrival process as MatrixCell [alpha, T] or [T] with uniform alpha
 * @param service PH service process as MatrixCell [beta, S] or [S] with uniform beta
 * @return QsysMapPhResult with performance metrics
 */
fun qsys_phph1(arrival: MatrixCell, service: MatrixCell): QsysMapPhResult {
    // Extract arrival PH parameters
    val (alpha, T) = extractPH(arrival, "arrival")

    // Extract service PH parameters
    val (beta, S) = extractPH(service, "service")

    return qsys_phph1(alpha, T, beta, S)
}

/**
 * Extracts PH parameters from MatrixCell.
 *
 * @param ph MatrixCell containing [alpha, T] or just [T]
 * @param name Name for error messages
 * @return Pair of (alpha, T)
 */
private fun extractPH(ph: MatrixCell, name: String): Pair<Matrix, Matrix> {
    require(ph.size() >= 1) { "$name PH must have at least 1 matrix" }

    return if (ph.size() >= 2) {
        // Explicit alpha provided
        val alpha = if (ph[0].numRows == 1) ph[0] else ph[0].transpose()
        val T = ph[1]
        Pair(alpha, T)
    } else {
        // Only T provided, use uniform initial distribution
        val T = ph[0]
        val n = T.numRows
        val alpha = Matrix.ones(1, n).scale(1.0 / n)
        Pair(alpha, T)
    }
}

/**
 * Analyzes a PH/PH/1 queue with exponential service (simplified E/M/1).
 *
 * @param alpha Arrival PH initial probability vector
 * @param T Arrival PH generator matrix
 * @param mu Exponential service rate
 * @return QsysMapPhResult with performance metrics
 */
fun qsys_phm1(alpha: Matrix, T: Matrix, mu: Double): QsysMapPhResult {
    // Create exponential service PH
    val beta = Matrix(1, 1)
    beta[0, 0] = 1.0
    val S = Matrix(1, 1)
    S[0, 0] = -mu

    return qsys_phph1(alpha, T, beta, S)
}

/**
 * Queueing system phph1 algorithms
 */
@Suppress("unused")
class QsysPhph1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}
