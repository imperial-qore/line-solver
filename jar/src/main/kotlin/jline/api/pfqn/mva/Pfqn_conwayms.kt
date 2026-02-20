/**
 * @file Conway-Maxwell approximate MVA for multi-server queueing networks
 * 
 * Implements the Conway-Maxwell approximation method for analyzing closed queueing networks
 * with multi-server stations. Uses iterative estimation and forward MVA techniques to
 * handle complex multi-server dependencies with various scheduling disciplines.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.io.Ret.*
import jline.GlobalConstants
import jline.lang.constant.SchedStrategy
import jline.util.Maths
import jline.util.PopulationLattice
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath
import java.util.*

fun pfqn_conwayms(L: Matrix, N: Matrix, Z: Matrix, nservers: IntArray): pfqnAMVAMS {
    val M = L.numRows
    val sched = arrayOfNulls<SchedStrategy>(M)
    Arrays.fill(sched, SchedStrategy.FCFS)
    val tol = GlobalConstants.FineTol
    val iter = 1000
    return pfqn_conwayms(L, N, Z, nservers, sched, tol, iter)
}

fun pfqn_conwayms(L: Matrix, N: Matrix, Z: Matrix): pfqnAMVAMS {
    val M = L.numRows
    val nservers = Matrix.ones(1, M).toIntArray1D()
    val sched = arrayOfNulls<SchedStrategy>(M)
    Arrays.fill(sched, SchedStrategy.FCFS)
    val tol = GlobalConstants.FineTol
    val iter = 1000
    return pfqn_conwayms(L, N, Z, nservers, sched, tol, iter)
}

fun pfqn_conwayms(L: Matrix,
                  N: Matrix,
                  Z: Matrix,
                  nservers: IntArray,
                  type: Array<SchedStrategy?>,
                  tol: Double,
                  maxiter: Int): pfqnAMVAMS {
    var Z = Z
    val M = L.numRows
    val R = L.numCols

    Z = Z.sumCols()

    val Q = MatrixCell(1 + R)
    val PB = MatrixCell(1 + R)
    val P = MatrixCell(1 + R)
    val Delta = MatrixCell(1 + R)
    val maxNumServers = Arrays.stream(nservers).max().orElse(1)
    for (s in 0..<1 + R) {
        Q[s] = Matrix(M, R)
        P[s] = Matrix(M, maxNumServers)
        PB[s] = Matrix(M, 1 + R)
        Delta[s] = Matrix(R, R)
    }

    for (i in 0..<M) {
        for (r in 0..<R) {
            for (s in 0..R) {
                val N_1 = if (s == 0) N else Matrix.oner(N, s - 1)
                Q[s][i, r] = N_1[r] / M
            }
        }
    }

    for (i in 0..<M) {
        for (r in 0..<R) {
            for (s in 0..R) {
                val N_1 = if (s == 0) N else Matrix.oner(N, s - 1)
                val pop = N_1.elementSum()

                if (nservers[i] > 1) {
                    for (j in 1..<nservers[i]) {
                        P[s][i, j] = 2 * Q[s].getRow(i).elementSum() / (pop * (pop + 1))
                    }
                    PB[s][i] = 2 * Q[s].getRow(i).elementSum() / (pop + 1 - nservers[i]) / (pop * (pop + 1))
                    P[s][i, 0] = 1 - PB[s][i] - P[s].sumSubMatrix(i, i + 1, 1, nservers[i])
                }
            }
        }
    }

    var totiter = 0

    // Main loop
    for (I in 0..1) {
        for (s in 0..R) {
            val N_1 = if (s == 0) N else Matrix.oner(N, s - 1)
            // Core(N_1)
            val coreResult =
                pfqn_conwayms_core(L, M, R, N_1, Z, nservers, Q[s], P[s], PB[s], Delta, type, tol, maxiter - totiter)
            Q[s].setTo(coreResult.Q)
            P[s].setTo(coreResult.P)
            PB[s].setTo(coreResult.PB)
            totiter += coreResult.iter
        }
        // Update_Delta
        for (i in 0..<M) {
            for (r in 0..<R) {
                for (s in 0..<R) {
                    val Ns = Matrix.oner(N, s)

                    if (N[s] > 2) {
                        Delta[s][i, r] = Q[1 + s][i, r] / Ns[r] - Q[0][i, r] / N[r]
                    }
                }
            }
        }
    }

    val finalCoreResult = pfqn_conwayms_core(L, M, R, N, Z, nservers, Q[0], P[0], PB[0], Delta, type, tol, maxiter)
    val retQ = finalCoreResult.Q
    val retW = finalCoreResult.W
    val retX = finalCoreResult.X
    totiter += finalCoreResult.iter

    val retU = Matrix(M, R)
    for (i in 0..<M) {
        for (r in 0..<R) {
            if (nservers[i] == 1) {
                retU[i, r] = retX[r] * L[i, r]
            } else {
                retU[i, r] = retX[r] * L[i, r] / nservers[i]
            }
        }
    }

    val retC = N.copy().elementDiv(retX).sub(Z)

    return pfqnAMVAMS(retQ, retU, retW, retC, retX, totiter)
}

fun pfqn_conwayms_core(L: Matrix,
                       M: Int,
                       R: Int,
                       N_1: Matrix,
                       Z: Matrix,
                       nservers: IntArray,
                       Q: Matrix,
                       P: Matrix,
                       PB: Matrix,
                       Delta: MatrixCell,
                       type: Array<SchedStrategy?>,
                       tol: Double,
                       maxiter: Int): LinearizerResult {
    var Q = Q
    var P = P
    var PB = PB
    var hasConverged = false
    var W = L.copy()
    var T = Matrix.createLike(L)
    var iter = 1

    while (!hasConverged) {
        val Qlast = Q.copy()
        val estimateResult = pfqn_conwayms_estimate(M, R, N_1, nservers, Q, P, PB, Delta, W)

        val forwardMVAResult = pfqn_conwayms_forwardmva(L,
            M,
            R,
            N_1,
            Z,
            nservers,
            type,
            estimateResult.Q_1,
            estimateResult.P_1,
            estimateResult.PB_1,
            estimateResult.T_1)
        Q = forwardMVAResult.Q
        W = forwardMVAResult.W
        T = forwardMVAResult.T
        P = forwardMVAResult.P
        PB = forwardMVAResult.PB

        if (Q.copy().sub(Qlast).norm() < tol || iter > maxiter) {
            hasConverged = true
        }
        iter++
    }
    return LinearizerResult.withX(Q, W, T, P, PB, iter)
}


fun pfqn_conwayms_estimate(M: Int,
                           R: Int,
                           N_1: Matrix,
                           nservers: IntArray,
                           Q: Matrix,
                           P: Matrix,
                           PB: Matrix,
                           Delta: MatrixCell,
                           W: Matrix): pfqnEstimate {
    val maxNumServers = Arrays.stream(nservers).max().orElse(1)
    val P_1 = MatrixCell(1 + R)
    val Q_1 = MatrixCell(1 + R)
    for (r in 0..<R + 1) {
        Q_1[r] = Matrix(M, R)
        P_1[r] = Matrix(M, maxNumServers)
    }
    val PB_1 = Matrix(M, 1 + R)
    val T_1 = Matrix(R, 1 + R)

    for (i in 0..<M) {
        if (nservers[i] > 1) {
            for (j in 0..<nservers[i]) {
                for (s in 0..R) {
                    P_1[s][i, j] = P[i, j]
                }
            }
            for (s in 0..R) {
                PB_1[i, s] = PB[i, 0]
            }
        }
        for (r in 0..<R) {
            for (s in 0..<R) {
                val Ns = Matrix.oner(N_1, s)
                Q_1[1 + s][i, r] = Ns[r] * (Q[i, r] / N_1[r] + Delta[s][i, r])
            }
        }
    }
    for (r in 0..<R) {
        for (s in 0..<R) {
            val Nr = Matrix.oner(N_1, r)
            for (i in 0..<M) {
                if (W[i, s] > 0) {
                    T_1[s, 1 + r] = Nr[s] * (Q[i, s] / N_1[s] + Delta[s][i, r]) / W[i, s]
                    break
                }
            }
        }
    }
    return pfqnEstimate(Q_1, P_1, PB_1, T_1)
}


fun pfqn_conwayms_forwardmva(L: Matrix,
                             M: Int,
                             R: Int,
                             N_1: Matrix,
                             Z: Matrix,
                             nservers: IntArray,
                             type: Array<SchedStrategy?>,
                             Q_1: MatrixCell,
                             P_1: MatrixCell,
                             PB_1: Matrix,
                             T_1: Matrix): LinearizerResult {
    val maxNumServers = Arrays.stream(nservers).max().orElse(1)
    val W = Matrix(M, R)
    val T = Matrix(1, R)
    val Q = Matrix(M, R)
    val P = Matrix(M, maxNumServers)
    val PB = Matrix(M, 1)
    val XR = Matrix(M, R)
    val XE = MatrixCell(R)
    for (r in 0..<R) {
        XE[r] = Matrix(M, R)
    }

    val F = MatrixCell(R)
    for (r in 0..<R) {
        F[r] = Matrix(M, R)
        for (i in 0..<M) {
            val den = L.getRow(i).mult(T_1.getColumn(1 + r)).value()
            for (s in 0..<R) {
                F[r][i, s] = L[i, s] * T_1[s, 1 + r] / den
            }
        }
    }

    val mu = L.reciprocal()
    val C = Matrix(M, R + 1)
    for (i in 0..<M) {
        for (r in 0..<R) {
            if (nservers[i] > 1) {
                XR[i, r] = 0
                C[i, 1 + r] = 0
                var sprodRes = PopulationLattice.sprod(R, Matrix.singleton(nservers[i].toDouble()))
                while (sprodRes.s.value() >= 0) {
                    if (Matrix.compare(sprodRes.n.transpose(), Matrix.oner(N_1, r), "lte")) {
                        val n = sprodRes.n.copy().transpose()
                        val Ai = FastMath.exp(Maths.multinomialln(n) + n.mult(F[r].getRow(i).log().transpose()).value())
                        C[i, 1 + r] = C[i, 1 + r] + Ai
                        XR[i, r] = XR[i, r] + Ai / mu.getRow(i).mult(n.transpose()).value()
                    }
                    sprodRes = PopulationLattice.sprod(sprodRes.s, sprodRes.S, sprodRes.D)
                }
                XR[i, r] = XR[i, r] / C[i, 1 + r]
            }
        }
    }

    val Cx = Matrix(M, 1 + R)
    for (i in 0..<M) {
        for (r in 0..<R) {
            if (nservers[i] > 1) {
                for (c in 0..<R) {
                    XE[c][i, r] = 0
                    Cx[i, 1 + r] = 0
                    var sprodResult = PopulationLattice.sprod(R, Matrix.singleton(nservers[i].toDouble()))
                    while (sprodResult.s.value() >= 0) {
                        if (Matrix.compare(sprodResult.n.transpose(),
                                Matrix.oner(N_1, r),
                                "lte") && sprodResult.n[c] >= 1) {
                            val n = sprodResult.n.copy().transpose()
                            val Aix =
                                FastMath.exp(Maths.multinomialln(n) + n.mult(F[r].getRow(i).log().transpose()).value())
                            Cx[i, 1 + r] = Cx[i, 1 + r] + Aix
                            XE[c][i, r] = XE[c][i, r] + Aix / mu.getRow(i).mult(n.transpose()).value()
                        }
                        sprodResult = PopulationLattice.sprod(sprodResult.s, sprodResult.S, sprodResult.D)
                    }
                    XE[c][i, r] = XE[c][i, r] / Cx[i, 1 + r]
                }
            }
        }
    }

    for (i in 0..<M) {
        for (r in 0..<R) {
            if (nservers[i] == 1) {
                if (type[i] == SchedStrategy.FCFS) {
                    W[i, r] = L[i, r]
                    for (c in 0..<R) {
                        W[i, r] = W[i, r] + L[i, c] * Q_1[1 + r][i, c]
                    }
                } else {
                    W[i, r] = L[i, r]
                    for (c in 0..<R) {
                        W[i, r] = W[i, r] + L[i, r] * Q_1[1 + r][i, c]
                    }
                }
            } else {
                W[i, r] = L[i, r] + PB_1[i, 1 + r] * XR[i, r]
                for (c in 0..<R) {
                    W[i, r] = W[i, r] + XE[c][i, r] * (Q_1[1 + r][i, c] - L[i, c] * T_1[c, 1 + r])
                }
            }
        }
    }

    for (r in 0..<R) {
        T[r] = N_1[r] / (Z[r] + W.getColumn(r).elementSum())
        for (i in 0..<M) {
            Q[i, r] = T[r] * W[i, r]
        }
    }

    for (i in 0..<M) {
        if (nservers[i] > 1) {
            for (j in 0..<nservers[i] - 1) {
                for (c in 0..<R) {
                    P[i, 1 + j] = P[i, 1 + j] + L[i, c] * T[c] * P_1[1 + c][i, 1 + (j - 1)] / (j + 1)
                }
            }
        }
    }

    for (i in 0..<M) {
        if (nservers[i] > 1) {
            PB[i] = 0.0
            for (c in 0..<R) {
                PB[i] = PB[i] + L[i, c] * T[c] * (PB_1[i, 1 + c] + P_1[1 + c][i, nservers[i] - 1]) / nservers[i]
            }
        }
    }

    for (i in 0..<M) {
        if (nservers[i] > 1) {
            P[i, 0] = FastMath.max(0.0, 1 - PB[i])
            for (j in 1..<nservers[i]) {
                P[i, 0] = FastMath.max(0.0, P[i, 0] - P[i, j])
            }
        }
    }

    return LinearizerResult(Q, W, T, P, PB)
}
/**
 * PFQN conwayms algorithms
 */
@Suppress("unused")
class PfqnConwaymsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}