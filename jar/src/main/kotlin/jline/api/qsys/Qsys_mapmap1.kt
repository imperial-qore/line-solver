/**
 * @file MAP/MAP/1 queueing system analysis
 *
 * Implements analysis of MAP/MAP/1 queues using BUTools MMAPPH1FCFS solver.
 * The queue has Markovian Arrival Processes for both arrivals and service.
 * Service MAP is converted to equivalent PH representation for the solver.
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys

import jline.api.mc.ctmc_solve
import jline.lib.butools.MMAPPH1FCFS
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Analyzes a MAP/MAP/1 queue.
 *
 * Uses BUTools MMAPPH1FCFS by converting the service MAP to an equivalent
 * Phase-Type representation. The service time distribution is extracted
 * from the MAP embedded at arrival epochs.
 *
 * @param C0 Arrival MAP hidden transition matrix (n x n)
 * @param C1 Arrival MAP arrival transition matrix (n x n)
 * @param D0 Service MAP hidden transition matrix (m x m)
 * @param D1 Service MAP observable transition matrix (m x m)
 * @param numQLMoms Number of queue length moments to compute (default 3)
 * @param numQLProbs Number of queue length probabilities to compute (default 100)
 * @param numSTMoms Number of sojourn time moments to compute (default 3)
 * @return QsysMapPhResult with performance metrics
 */
@JvmOverloads
fun qsys_mapmap1(
    C0: Matrix,
    C1: Matrix,
    D0: Matrix,
    D1: Matrix,
    numQLMoms: Int = 3,
    numQLProbs: Int = 100,
    numSTMoms: Int = 3
): QsysMapPhResult {
    // Build arrival MMAP structure for BUTools (single class)
    val D = MatrixCell(2)
    D[0] = C0
    D[1] = C1

    // Convert service MAP to PH representation
    // The service time PH is derived from the MAP with:
    // - Generator: D0 (hidden transitions)
    // - Initial distribution: proportional to exit rates from arrival MAP
    val (sigma, S) = mapToPh(D0, D1)

    // Service parameters as maps (single class indexed by 0)
    val sigmaMap = HashMap<Int?, Matrix>()
    sigmaMap[0] = sigma
    val sMap = HashMap<Int?, Matrix>()
    sMap[0] = S

    // Call BUTools solver
    val result = MMAPPH1FCFS(
        D, sigmaMap, sMap,
        numQLMoms, numQLProbs, numSTMoms,
        null, false, false, null, null
    )

    // Compute utilization from arrival and service rates
    val thetaArr = ctmc_solve(C0.add(C1))
    val lambda = thetaArr.mult(C1).elementSum()

    val thetaSvc = ctmc_solve(D0.add(D1))
    val mu = thetaSvc.mult(D1).elementSum()
    val rho = lambda / mu

    // Extract queue length moments
    val ncMoms = result["ncMoms"]?.get(0)
    val meanQL: Double = ncMoms?.get(0, 0) ?: 0.0

    // Extract sojourn time moments
    val stMoms = result["stNoms"]?.get(0)
    val meanST: Double = stMoms?.get(0, 0) ?: 0.0

    // Mean service time from PH
    val negSinv = S.scale(-1.0).inv()
    val ones = Matrix.ones(S.numRows, 1)
    val meanService = sigma.mult(negSinv).mult(ones)[0, 0]

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
 * Converts a MAP to its equivalent PH representation.
 *
 * The PH representation captures the service time distribution
 * of the MAP, with initial distribution derived from the
 * stationary distribution weighted by observable transition rates.
 *
 * @param D0 MAP hidden transition matrix
 * @param D1 MAP observable transition matrix
 * @return Pair of (sigma, S) where sigma is initial distribution and S is generator
 */
private fun mapToPh(D0: Matrix, D1: Matrix): Pair<Matrix, Matrix> {
    // Stationary distribution of the MAP
    val theta = ctmc_solve(D0.add(D1))

    // Initial distribution for PH is proportional to
    // how customers enter the service process
    val sigma = theta.mult(D1)
    val total = sigma.elementSum()
    if (total > 0) {
        sigma.scaleEq(1.0 / total)
    } else {
        // Fallback to uniform if no arrivals
        for (i in 0 until sigma.numCols) {
            sigma[0, i] = 1.0 / sigma.numCols
        }
    }

    // The PH generator is the hidden transition matrix
    return Pair(sigma, D0)
}

/**
 * Simplified MAP/MAP/1 analysis using MatrixCell inputs.
 *
 * @param arrival MAP arrival process as MatrixCell [C0, C1]
 * @param service MAP service process as MatrixCell [D0, D1]
 * @return QsysMapPhResult with performance metrics
 */
fun qsys_mapmap1(arrival: MatrixCell, service: MatrixCell): QsysMapPhResult {
    require(arrival.size() >= 2) { "Arrival MAP must have at least 2 matrices [D0, D1]" }
    require(service.size() >= 2) { "Service MAP must have at least 2 matrices [D0, D1]" }

    return qsys_mapmap1(arrival[0], arrival[1], service[0], service[1])
}

/**
 * Queueing system mapmap1 algorithms
 */
@Suppress("unused")
class QsysMapmap1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}
