/**
 * @file Multi-server Krzesinski linearizer approximate MVA
 * 
 * Implements the multi-server version of Krzesinski's linearizer approximation for closed
 * queueing networks with multi-server stations. Uses Conway's algorithm description with
 * De Souza-Muntz adjustments for handling complex multi-server dependencies.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.io.Ret
import jline.GlobalConstants
import jline.lang.constant.SchedStrategy
import jline.util.matrix.Matrix

/**
 * Mean Value Analysis (MVA) method for closed networks with load dependent service
 *
 * @param L         - service demand matrix
 * @param N         - population vector
 * @param Z         - think times
 * @param mu        - load dependent service rate. mu[i][j] - load dependent service rate of the ith node when there are
 * j jobs in it
 * @param stabilize - whether to stabilize the probabilities or not (ensures that probabilities do not become negative)
 * @return - performance measures for the closed network.
 *
 * Multiserver version of Krzesinski's Linearizer
 */

fun pfqn_linearizerms(L: Matrix, N: Matrix, Z: Matrix = Matrix(1, L.numCols), nservers: Matrix): Ret.pfqnAMVAMS {
    L.numRows
    val type: MutableList<SchedStrategy> = ArrayList()
    type.add(SchedStrategy.PS)
    return pfqn_linearizerms(L, N, Z, nservers, type)
}

/**
 * Multiserver version of Krzesinski's Linearizer as described in Conway 1989, Fast Approximate Solution of
 * Queueing Networks with Multi-Server Chain- Dependent FCFS Queues. Minor adjustments based on De Souza-Muntz's
 * description of the algorithm.
 *
 * @param L        - service demand matrix
 * @param N        - population vector
 * @param Z        - think times
 * @param nservers - number of servers at each station
 * @param type     - scheduling discipline at each station
 * @param tol      - max tolerance admitted between successive iterations
 * @param maxiter  - maximum number of iterations
 * @return - the performance measures of the network.
 *
 *
 * Multiserver version of Krzesinski's Linearizer
 */
fun pfqn_linearizerms(L: Matrix,
                      N: Matrix,
                      Z: Matrix = Matrix(1, L.numCols),
                      nservers: Matrix,
                      type: List<SchedStrategy>,
                      tol: Double = GlobalConstants.FineTol,
                      maxiter: Int = 1000): Ret.pfqnAMVAMS {
    var Z = Z
    val M = L.numRows
    val R = L.numCols

    Z = if (Z.isEmpty) Matrix(1, R) else Z

    // Initialise
    val Q = arrayOfNulls<Matrix>(M)
    val PB = Matrix(M, 1 + R)
    val P = arrayOfNulls<Matrix>(M)
    val Delta = arrayOfNulls<Matrix>(M)
    for (i in 0..<M) {
        Q[i] = Matrix(R, 1 + R)
        P[i] = Matrix(nservers.elementMax().toInt(), 1 + R)
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

    for (i in 0..<M) {
        for (r in 0..<R) {
            for (s in -1..<R) {
                val N_1 = Matrix.oner(N, ArrayList(listOf(s)))
                val pop = N_1.elementSum()
                if (nservers[i] > 1) {
                    var sumQ = 0.0
                    for (k in 0..<R) {
                        sumQ += Q[i]!![k, 1 + s]
                    }
                    var j = 0
                    while (j < nservers[i] - 1) {
                        P[i]!![1 + j, 1 + s] = 2 * sumQ / (pop * (pop + 1))
                        j++
                    }
                    PB[i, 1 + s] = 2 * sumQ / (pop + 1 - nservers[i]) / (pop * (pop + 1))
                    var sumP = 0.0
                    var k = 0
                    while (k < nservers[i] - 1) {
                        sumP += P[i]!![k + 1, 1 + s]
                        k++
                    }
                    P[i]!![0, 1 + s] = 1 - PB[i, 1 + s] - sumP
                }
            }
        }
    }

    var totiter = 0

    // Main loop
    for (I in 0..1) {
        for (s in -1..<R) {
            val N_1 = Matrix.oner(N, ArrayList(listOf(s)))
            val Q1 = Matrix(M, R)
            val P1 = Matrix(M, nservers.elementMax().toInt())
            val PB1 = Matrix(M, 1)
            for (i in 0..<M) {
                for (j in 0..<R) {
                    Q1[i, j] = Q[i]!![j, 1 + s]
                }
                for (j in 0..<P1.numCols) {
                    P1[i, j] = P[i]!![j, 1 + s]
                }
                PB1[i, 0] = PB[i, 1 + s]
            }
            val ret1 = linearizerms_core(L, M, R, N_1, Z, nservers, Q1, P1, PB1, Delta, type, tol, maxiter - totiter)
            for (i in 0..<M) {
                for (j in 0..<R) {
                    Q[i]!![j, 1 + s] = ret1.Q[i, j]
                }
                var j = 0
                while (j < nservers.elementMax()) {
                    P[i]!![j, 1 + s] = ret1.P[i, j]
                    j++
                }
                PB[i, 1 + s] = ret1.PB[i]
            }
            totiter += ret1.iter
        }
        // Upgrade delta
        for (i in 0..<M) {
            for (r in 0..<R) {
                for (s in 0..<R) {
                    val Ns = Matrix.oner(N, ArrayList(listOf(s)))
                    Delta[i]!![r, s] = Q[i]!![r, 1 + s] / Ns[r] - Q[i]!![r, 0] / N[r]
                }
            }
        }
    }

    val Q1 = Matrix(M, R)
    val P1 = Matrix(M, nservers.elementMax().toInt())
    val PB1 = Matrix(M, 1)
    for (i in 0..<M) {
        for (j in 0..<R) {
            Q1[i, j] = Q[i]!![j, 0]
        }
        for (j in 0..<P1.numCols) {
            P1[i, j] = P[i]!![j, 0]
        }
        PB1[i, 0] = PB[i, 0]
    }
    val ret1 = linearizerms_core(L, M, R, N, Z, nservers, Q1, P1, PB1, Delta, type, tol, maxiter - totiter)
    totiter += ret1.iter
    val newQ = ret1.Q
    val W = ret1.W
    val X = ret1.T
    val U = Matrix(M, R)

    for (i in 0..<M) {
        for (r in 0..<R) {
            if (nservers[i] == 1.0) {
                U[i, r] = X[r] * L[i, r]
            } else {
                U[i, r] = X[r] * L[i, r] / nservers[i]
            }
        }
    }
    val C = Matrix(1, R)
    for (i in 0..<R) {
        C[0, i] = N[i] / X[i] - Z[i]
    }
    return Ret.pfqnAMVAMS(newQ, U, W, C, X, totiter)
}

internal fun linearizerms_core(L: Matrix,
                               M: Int,
                               R: Int,
                               N_1: Matrix,
                               Z: Matrix,
                               nservers: Matrix,
                               Q: Matrix,
                               P: Matrix,
                               PB: Matrix,
                               Delta: Array<Matrix?>,
                               type: List<SchedStrategy>,
                               tol: Double,
                               maxiter: Int): Ret.LinearizerResult {
    var Q = Q
    var P = P
    var PB = PB
    var iter = 0
    var hasConverged = false
    var W: Matrix? = null
    var T: Matrix? = null
    while (!hasConverged) {
        iter++
        val Qlast = Matrix(Q)
        // Estimate population
        val ret1 = linearizerms_estimate(M, R, N_1, nservers, Q, P, PB, Delta)
        // Forward MVA
        val ret2 = linearizerms_forwardMVA(L, M, R, N_1, Z, nservers, type, ret1.Q_1, ret1.P_1, ret1.PB_1)
        Q = ret2.Q
        W = ret2.W
        T = ret2.T
        P = ret2.P
        PB = ret2.PB
        if (Q.sub(Qlast).norm() < tol || iter > maxiter) {
            hasConverged = true
        }
    }
    return Ret.LinearizerResult(Q, W, T, P, PB, iter)
}

internal fun linearizerms_estimate(M: Int,
                                   R: Int,
                                   N_1: Matrix,
                                   nservers: Matrix,
                                   Q: Matrix,
                                   P: Matrix,
                                   PB: Matrix,
                                   Delta: Array<Matrix?>): Ret.pfqnLinearizerMSEstimate {
    val P_1 = arrayOfNulls<Matrix>(M)
    val Q_1 = arrayOfNulls<Matrix>(M)
    for (i in 0..<M) {
        P_1[i] = Matrix(nservers.elementMax().toInt(), 1 + R)
        Q_1[i] = Matrix(R, 1 + R)
    }
    val PB_1 = Matrix(M, 1 + R)
    for (i in 0..<M) {
        if (nservers[i] > 1) {
            var j = -1
            while (j < nservers[i] - 1) {
                for (s in -1..<R) {
                    P_1[i]!![1 + j, 1 + s] = P[i, 1 + j]
                }
                j++
            }
            for (s in -1..<R) {
                PB_1[i, 1 + s] = PB[i, 0]
            }
        }
        for (r in 0..<R) {
            for (s in 0..<R) {
                val Ns = Matrix.oner(N_1, ArrayList(listOf(s)))
                Q_1[i]!![r, 1 + s] = Ns[r] * (Q[i, r] / N_1[r] + Delta[i]!![r, s])
            }
        }
    }
    return Ret.pfqnLinearizerMSEstimate(Q_1, P_1, PB_1)
}

internal fun linearizerms_forwardMVA(L: Matrix,
                                     M: Int,
                                     R: Int,
                                     N_1: Matrix,
                                     Z: Matrix,
                                     nservers: Matrix,
                                     type: List<SchedStrategy>,
                                     Q_1: Array<Matrix>,
                                     P_1: Array<Matrix>,
                                     PB_1: Matrix): Ret.LinearizerResult {
    val W = Matrix(M, R)
    val T = Matrix(1, R)
    val Q = Matrix(M, R)
    val P = Matrix(M, nservers.elementMax().toInt())
    val PB = Matrix(M, 1)
    for (i in 0..<M) {
        for (r in 0..<R) {
            W[i, r] = L[i, r] / nservers[i]
            if (L[i, r] == 0.0) {
                // 0 service demand at this station => this class does not visit the current node
                continue
            }
            var flag = true // flag = whether all nodes have the FCFS scheduling strategy
            for (k in 0..<M) {
                if (type[k] == SchedStrategy.FCFS) {
                    flag = false
                    break
                }
            }
            if (flag) {
                for (s in 0..<R) {
                    W[i, r] = W[i, r] + (L[i, s] / nservers[i]) * Q_1[i][s, 1 + r]
                }
            } else {
                for (s in 0..<R) {
                    W[i, r] = W[i, r] + (L[i, r] / nservers[i]) * Q_1[i][s, 1 + r]
                }
            }
            if (nservers[i] > 1) {
                var j = 0
                while (j <= nservers[i] - 2) {
                    if (flag) {
                        for (s in 0..<R) {
                            W[i, r] = W[i, r] + L[i, s] * (nservers[i] - 1 - j) * P_1[i][j, 1 + r]
                        }
                    } else {
                        for (s in 0..<R) {
                            W[i, r] = W[i, r] + L[i, r] * (nservers[i] - 1 - j) * P_1[i][j, 1 + r]
                        }
                    }
                    j++
                }
            }
        }
    }
    for (r in 0..<R) {
        T[r] = N_1[r] / (Z[r] + Matrix.extractColumn(W, r, null).elementSum())
        for (i in 0..<M) {
            Q[i, r] = T[r] * W[i, r]
        }
    }
    for (i in 0..<M) {
        if (nservers[i] > 1) {
            for (k in 0..<P.numCols) {
                P[i, k] = 0
            }
            var j = 1
            while (j <= nservers[i] - 1) {
                for (s in 0..<R) {
                    P[i, j] = P[i, j] + L[i, s] * T[s] * P_1[i][j - 1, 1 + s] / j
                }
                j++
            }
        }
    }
    for (i in 0..<M) {
        if (nservers[i] > 1) {
            PB[i, 0] = 0
            for (s in 0..<R) {
                PB[i, 0] =
                    PB[i, 0] + L[i, s] * T[s] * (PB_1[i, 1 + s] + P_1[i][nservers[i].toInt() - 1, 1 + s]) / nservers[i]
            }
        }
    }
    for (i in 0..<M) {
        if (nservers[i] > 1) {
            P[i, 0] = 1 - PB[i]
            var j = 0
            while (j < nservers[i] - 1) {
                P[i, 0] = P[i, 0] - P[i, 1 + j]
                j++
            }
        }
    }
    return Ret.LinearizerResult(Q, W, T, P, PB)
}
/**
 * PFQN linearizerms algorithms
 */
@Suppress("unused")
class PfqnLinearizermsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}