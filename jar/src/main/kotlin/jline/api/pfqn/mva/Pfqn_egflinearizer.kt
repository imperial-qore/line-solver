/**
 * @file Extended General-Form linearizer approximate MVA with class-specific parameters
 * 
 * Implements the Extended General-Form linearizer approximation for closed queueing networks
 * with class-specific linearization parameters (alpha matrix). Provides enhanced accuracy
 * over standard linearizer methods through flexible per-class parameter tuning.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.io.Ret
import jline.lang.constant.SchedStrategy
import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Extended general form linearizer approximate mean value analysis algorithm.
 * This method performs approximate mean value analysis using an extended general
 * form linearizer approach.
 *
 * @param L       - the service demand matrix, where rows represent stations and columns represent classes
 * @param N       - the population vector indicating the number of requests for each class
 * @param Z       - the think times, representing the average time between completion of one service request and the beginning of the next
 * @param type    - the types of scheduling disciplines (e.g., FCFS, PS) used at each station
 * @param tol     - the maximum tolerance admitted between successive iterations, used for convergence checking
 * @param maxiter - the maximum number of iterations allowed in the algorithm
 * @param alpha   - matrix of alphas, which provides weightings or adjustments for each class in the linearizer
 * @return - a pfqnAMVA object containing the computed performance measures: queue lengths, utilizations, response times, and throughput rates
 */

fun pfqn_egflinearizer(L: Matrix,
                       N: Matrix,
                       Z: Matrix = Matrix(1, L.numCols),
                       type: Array<SchedStrategy>,
                       tol: Double,
                       maxiter: Int,
                       alpha: Matrix): Ret.pfqnAMVA {
    var Z = Z
    val M = L.numRows
    val R = L.numCols

    Z = if (Z.isEmpty) Matrix(1, R) else Z
    Z = Z.sumCols()
    var Lmaxcols = true
    for (j in 0..<L.numCols) {
        var maxcol = 0.0
        for (i in 0..<L.numRows) {
            if (L[i, j] > maxcol) maxcol = L[i, j]
        }
        if (maxcol != 0.0) {
            Lmaxcols = false
            break
        }
    }
    if (L.isEmpty || Lmaxcols) {
        val X = Matrix(N.numRows, N.numCols)
        for (i in 0..<X.numRows) {
            for (j in 0..<X.numCols) {
                X[i, j] = N[i, j] / Z[i, j]
            }
        }
        val Q = Matrix(M, R)
        val U = Matrix(M, R)
        val W = Matrix(M, R)
        val T = Matrix(M, R)
        val C = Matrix(1, R)
        for (r in 0..<R) {
            for (i in 0..<M) {
                U[i, r] = X[r] * L[i, r]
            }
        }
        for (r in 0..<R) {
            for (i in 0..<M) {
                T[i, r] = Q[i, r] / W[i, r]
            }
        }
        val totiter = 0
        return Ret.pfqnAMVA(Q, U, W, T, C, X, totiter)
    }
    // Initialise
    val Q = arrayOfNulls<Matrix>(M)
    val Delta = arrayOfNulls<Matrix>(M)
    for (i in 0..<M) {
        Q[i] = Matrix(R, 1 + R)
        Delta[i] = Matrix(R, R)
    }

    for (s in -1..<R) {
        val N_1 = Matrix.oner(N, ArrayList(listOf(s)))
        val init0 = pfqn_bs(L, N_1, Z)
        for (i in 0..<M) {
            for (r in 0..<R) {
                Q[i]!![r, 1 + s] = init0.Q[i, r]
            }
        }
    }

    var totiter = 0

    // Main loop
    for (I in 0..2) {
        for (s in -1..<R) {
            val N_1 = Matrix.oner(N, ArrayList(listOf(s)))
            val Q1 = Matrix(M, R)
            for (i in 0..<M) {
                for (j in 0..<R) {
                    Q1[i, j] = Q[i]!![j, 1 + s]
                }
            }
            val ret1 = egflinearizer_core(L, M, R, N_1, Z, Q1, Delta, type, tol, maxiter - totiter, alpha)
            for (i in 0..<M) {
                for (j in 0..<R) {
                    Q[i]!![j, 1 + s] = ret1.Q[i, j]
                }
            }

            totiter += ret1.iter
        }
        // Upgrade delta
        for (i in 0..<M) {
            for (r in 0..<R) {
                if (N[r] == 1.0) {
                    for (j in 0..<R) {
                        Q[i]!![j, 1 + r] = 0
                    }
                }
                for (s in 0..<R) {
                    if (N[s] > 1) {
                        val Ns = Matrix.oner(N, ArrayList(listOf(s)))
                        Delta[i]!![r, s] =
                            Q[i]!![r, 1 + s] / FastMath.pow(Ns[r], alpha[r]) - Q[i]!![r, 0] / FastMath.pow(N[r],
                                alpha[r])
                    }
                }
            }
        }
    }

    // Core(N)
    val Q1 = Matrix(M, R)
    for (i in 0..<M) {
        for (j in 0..<R) {
            Q1[i, j] = Q[i]!![j, 0]
        }
    }
    val ret1 = egflinearizer_core(L, M, R, N, Z, Q1, Delta, type, tol, maxiter - totiter, alpha)
    val newQ = ret1.Q
    val W = ret1.W
    val X = ret1.T
    totiter += ret1.iter
    // Compute performance metrics
    val U = Matrix(M, R)
    for (i in 0..<M) {
        for (r in 0..<R) {
            U[i, r] = X[r] * L[i, r]
        }
    }
    val C = Matrix(1, R)
    for (i in 0..<R) {
        C[0, i] = N[i] / X[i] - Z[i]
    }
    return Ret.pfqnAMVA(newQ, U, W, null, C, X, totiter)
}

/**
 * Core method for the extended general form linearizer. This method is responsible for
 * iterating the linearizer steps, adjusting the queue length estimates, and checking for convergence.
 *
 * @param L       - the service demand matrix
 * @param M       - the number of stations
 * @param R       - the number of classes
 * @param N_1     - adjusted population vector
 * @param Z       - the think times
 * @param Q       - current queue length estimates
 * @param Delta   - array of matrices representing adjustments for the linearizer
 * @param type    - scheduling disciplines at each station
 * @param tol     - convergence tolerance
 * @param maxiter - maximum number of iterations
 * @param alpha   - matrix of alpha values for adjustments
 * @return - a pfqnLinearizerCore object containing the updated queue lengths, response times, throughputs, and iteration count
 */
internal fun egflinearizer_core(L: Matrix,
                                M: Int,
                                R: Int,
                                N_1: Matrix,
                                Z: Matrix,
                                Q: Matrix,
                                Delta: Array<Matrix?>,
                                type: Array<SchedStrategy>,
                                tol: Double,
                                maxiter: Int,
                                alpha: Matrix): Ret.LinearizerResult {
    var Q = Q
    var hasConverged = false
    var W = Matrix(L)
    var iter = 0
    var T: Matrix? = null
    while (!hasConverged) {
        val Qlast = Matrix(Q)
        // Estimate population
        val ret1 = egflinearizer_estimate(L, M, R, N_1, Z, Q, Delta, W, alpha)
        val Q_1 = ret1.Q_1
        // Forward MVA
        val ret2 = egflinearizer_forwardMVA(L, M, R, type, N_1, Z, Q_1)
        Q = ret2.Q
        W = ret2.W
        T = ret2.T
        if (Q.sub(Qlast).norm() < tol || iter > maxiter) {
            hasConverged = true
        }
        iter++
    }
    return Ret.LinearizerResult(Q, W, T, iter)
}

/**
 * Estimate method for the extended general form linearizer. This method computes
 * the intermediate estimates of the queue lengths and response times.
 *
 * @param L     - the service demand matrix
 * @param M     - the number of stations
 * @param R     - the number of classes
 * @param N_1   - adjusted population vector
 * @param Z     - the think times
 * @param Q     - current queue length estimates
 * @param Delta - array of matrices representing adjustments for the linearizer
 * @param W     - response time matrix
 * @param alpha - matrix of alpha values for adjustments
 * @return - a pfqnLinearizerEstimate object containing the intermediate queue length estimates
 */
internal fun egflinearizer_estimate(L: Matrix?,
                                    M: Int,
                                    R: Int,
                                    N_1: Matrix,
                                    Z: Matrix?,
                                    Q: Matrix,
                                    Delta: Array<Matrix?>,
                                    W: Matrix?,
                                    alpha: Matrix): Ret.pfqnLinearizerEstimate {
    val Q_1 = arrayOfNulls<Matrix>(M)
    for (i in 0..<M) {
        Q_1[i] = Matrix(R, 1 + R)
    }
    val T_1 = Matrix(R, 1 + R)
    for (i in 0..<M) {
        for (r in 0..<R) {
            for (s in 0..<R) {
                val Ns = Matrix.oner(N_1, ArrayList(listOf(s)))
                Q_1[i]!![r, 1 + s] =
                    FastMath.pow(Ns[r], alpha[r]) * (Q[i, r] / FastMath.pow(N_1[r], alpha[r]) + Delta[i]!![r, s])
            }
        }
    }
    return Ret.pfqnLinearizerEstimate(Q_1, T_1)
}

/**
 * Forward Mean Value Analysis method for the extended general form linearizer.
 * This method computes the response times, throughputs, and queue lengths
 * using the mean value analysis technique.
 *
 * @param L    - the service demand matrix
 * @param M    - the number of stations
 * @param R    - the number of classes
 * @param type - scheduling disciplines at each station
 * @param N_1  - adjusted population vector
 * @param Z    - the think times
 * @param Q_1  - intermediate queue length estimates
 * @return - a LinearizerResult object containing the computed queue lengths, response times, and throughputs
 */
internal fun egflinearizer_forwardMVA(L: Matrix,
                                      M: Int,
                                      R: Int,
                                      type: Array<SchedStrategy>,
                                      N_1: Matrix,
                                      Z: Matrix,
                                      Q_1: Array<Matrix>): Ret.LinearizerResult {
    val W = Matrix(M, R)
    val T = Matrix(1, R)
    val Q = Matrix(M, R)

    // Compute residence time
    for (i in 0..<M) {
        for (r in 0..<R) {
            if (type[i] == SchedStrategy.FCFS) {
                W[i, r] = L[i, r]
                if (L[i, r] != 0.0) {
                    for (s in 0..<R) {
                        W[i, r] = W[i, r] + L[i, s] * Q_1[i][s, 1 + r]
                    }
                }
            } else {
                var sumQ = 0.0
                for (k in 0..<Q_1[i].numRows) {
                    sumQ += Q_1[i][k, 1 + r]
                }
                W[i, r] = L[i, r] * (1 + sumQ)
            }
        }
    }

    // Compute throughputs and queue lengths
    for (r in 0..<R) {
        T[r] = N_1[r] / (Z[r] + Matrix.extractColumn(W, r, null).elementSum())
        for (i in 0..<M) {
            Q[i, r] = T[r] * W[i, r]
        }
    }
    return Ret.LinearizerResult(Q, W, T)
}
/**
 * PFQN egflinearizer algorithms
 */
@Suppress("unused")
class PfqnEgflinearizerAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}