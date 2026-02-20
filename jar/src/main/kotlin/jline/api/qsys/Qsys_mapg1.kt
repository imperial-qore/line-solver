/**
 * @file MAP/G/1 queueing system analysis
 *
 * Implements analysis of MAP/G/1 queues using BUTools MMAPPH1FCFS solver.
 * The queue has a Markovian Arrival Process (MAP) for arrivals and a
 * General service time distribution specified by its moments.
 *
 * The general service time is fitted to a Phase-Type distribution using
 * moment matching before analysis.
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys

import jline.api.mc.ctmc_solve
import jline.lib.butools.APHFrom3Moments
import jline.lib.butools.MMAPPH1FCFS
import jline.lib.butools.ph.ph2From3Moments
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Analyzes a MAP/G/1 queue.
 *
 * Uses BUTools MMAPPH1FCFS by fitting the general service time distribution
 * to a Phase-Type distribution via moment matching. Supports 2 or 3 moment
 * matching depending on the number of moments provided.
 *
 * @param D0 MAP hidden transition matrix (n x n)
 * @param D1 MAP arrival transition matrix (n x n)
 * @param serviceMoments First k raw moments of service time [E[S], E[S^2], ...] (k = 2 or 3)
 * @param numQLMoms Number of queue length moments to compute (default 3)
 * @param numQLProbs Number of queue length probabilities to compute (default 100)
 * @param numSTMoms Number of sojourn time moments to compute (default 3)
 * @return QsysMapPhResult with performance metrics
 */
@JvmOverloads
fun qsys_mapg1(
    D0: Matrix,
    D1: Matrix,
    serviceMoments: DoubleArray,
    numQLMoms: Int = 3,
    numQLProbs: Int = 100,
    numSTMoms: Int = 3
): QsysMapPhResult {
    require(serviceMoments.size >= 2) { "At least 2 service moments are required" }

    // Fit service distribution to PH using moment matching
    val (sigma, S) = fitServiceToPH(serviceMoments)

    // Build arrival MMAP structure for BUTools (single class)
    val D = MatrixCell(2)
    D[0] = D0
    D[1] = D1

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
    val theta = ctmc_solve(D0.add(D1))
    val lambda = theta.mult(D1).elementSum()

    val meanService = serviceMoments[0]
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
 * Fits a general service time distribution to a Phase-Type representation
 * using moment matching.
 *
 * @param moments Raw moments [E[X], E[X^2], E[X^3], ...]
 * @return Pair of (alpha, A) where alpha is initial distribution and A is generator
 */
private fun fitServiceToPH(moments: DoubleArray): Pair<Matrix, Matrix> {
    val m1 = moments[0]

    // Check if we have enough moments for PH fitting
    return if (moments.size >= 3) {
        // Use APH from 3 moments for best accuracy
        try {
            val aph = APHFrom3Moments(moments)
            val alpha = aph.initProb
            val T = aph.getParam(3).value as Matrix
            Pair(alpha, T)
        } catch (e: Exception) {
            // Fallback to 2-phase PH
            try {
                val ph2 = ph2From3Moments(moments)
                Pair(ph2.alpha, ph2.A)
            } catch (e2: Exception) {
                // Ultimate fallback to exponential
                createExponentialPH(m1)
            }
        }
    } else if (moments.size >= 2) {
        // Use 2 moments - create PH(2) with matching mean and variance
        val m2 = moments[1]
        val cv2 = m2 / (m1 * m1) - 1.0  // Squared coefficient of variation

        when {
            cv2 <= 0.0 -> {
                // Deterministic or near-deterministic: use Erlang approximation
                createErlangPH(m1, maxOf(1, (1.0 / maxOf(cv2, 0.01)).toInt()))
            }
            cv2 < 1.0 -> {
                // Hypoexponential: use Erlang-2 approximation scaled
                createErlang2PH(m1, cv2)
            }
            cv2 == 1.0 -> {
                // Exponential
                createExponentialPH(m1)
            }
            else -> {
                // Hyperexponential: use 2-phase hyperexponential
                createHyperexp2PH(m1, cv2)
            }
        }
    } else {
        // Only mean provided - use exponential
        createExponentialPH(m1)
    }
}

/**
 * Creates an exponential PH distribution.
 */
private fun createExponentialPH(mean: Double): Pair<Matrix, Matrix> {
    val alpha = Matrix(1, 1)
    alpha[0, 0] = 1.0
    val A = Matrix(1, 1)
    A[0, 0] = -1.0 / mean
    return Pair(alpha, A)
}

/**
 * Creates an Erlang-k PH distribution.
 */
private fun createErlangPH(mean: Double, k: Int): Pair<Matrix, Matrix> {
    val mu = k.toDouble() / mean
    val alpha = Matrix(1, k)
    alpha[0, 0] = 1.0

    val A = Matrix(k, k)
    for (i in 0 until k) {
        A[i, i] = -mu
        if (i < k - 1) {
            A[i, i + 1] = mu
        }
    }
    return Pair(alpha, A)
}

/**
 * Creates an Erlang-2 based hypoexponential with matched mean and cv2.
 */
private fun createErlang2PH(mean: Double, cv2: Double): Pair<Matrix, Matrix> {
    // For cv2 < 1, use generalized Erlang
    // Solve for parameters that match mean and cv2
    val k = maxOf(2, (1.0 / cv2).toInt())
    return createErlangPH(mean, k)
}

/**
 * Creates a 2-phase hyperexponential with matched mean and cv2.
 */
private fun createHyperexp2PH(mean: Double, cv2: Double): Pair<Matrix, Matrix> {
    // For cv2 > 1, use balanced hyperexponential
    // H2 with parameters (p, lambda1, lambda2) matching mean and cv2

    // Balanced means approach
    val cv = kotlin.math.sqrt(cv2)
    val p = 0.5 * (1.0 + kotlin.math.sqrt((cv2 - 1.0) / (cv2 + 1.0)))
    val lambda1 = 2.0 * p / mean
    val lambda2 = 2.0 * (1.0 - p) / mean

    val alpha = Matrix(1, 2)
    alpha[0, 0] = p
    alpha[0, 1] = 1.0 - p

    val A = Matrix(2, 2)
    A[0, 0] = -lambda1
    A[1, 1] = -lambda2

    return Pair(alpha, A)
}

/**
 * Analyzes a MAP/G/1 queue with service time specified as mean and coefficient of variation.
 *
 * @param D0 MAP hidden transition matrix
 * @param D1 MAP arrival transition matrix
 * @param meanService Mean service time
 * @param cvService Coefficient of variation of service time
 * @return QsysMapPhResult with performance metrics
 */
fun qsys_mapg1(
    D0: Matrix,
    D1: Matrix,
    meanService: Double,
    cvService: Double
): QsysMapPhResult {
    // Compute first two moments from mean and CV
    val m1 = meanService
    val m2 = m1 * m1 * (1.0 + cvService * cvService)

    return qsys_mapg1(D0, D1, doubleArrayOf(m1, m2))
}

/**
 * Analyzes a MAP/G/1 queue using MatrixCell input for arrival.
 *
 * @param arrival MAP arrival process as MatrixCell [D0, D1]
 * @param serviceMoments First k raw moments of service time
 * @return QsysMapPhResult with performance metrics
 */
fun qsys_mapg1(arrival: MatrixCell, serviceMoments: DoubleArray): QsysMapPhResult {
    require(arrival.size() >= 2) { "Arrival MAP must have at least 2 matrices [D0, D1]" }
    return qsys_mapg1(arrival[0], arrival[1], serviceMoments)
}

/**
 * Queueing system mapg1 algorithms
 */
@Suppress("unused")
class QsysMapg1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}
