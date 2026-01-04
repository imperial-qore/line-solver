/**
 * @file Linearizer++ approximate MVA with higher-order moment corrections
 * 
 * Implements the Linearizer++ algorithm for closed queueing networks with enhanced accuracy
 * through higher-order moment corrections. Supports multiple approximation levels with
 * configurable precision-performance trade-offs for complex multi-class systems.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.io.Ret.pfqnAMVA
import jline.util.matrix.Matrix
import kotlin.math.abs
import kotlin.math.max

/**
 * Linearizer approximate mean value analysis algorithm
 *
 * @param L       - the service demand matrix
 * @param N       - the population vector
 * @param Z       - the think times
 * @param type    - the types of the scheduling disciplines at each station
 * @param tol     - max tolerance admitted between successive iterations
 * @param maxiter - maximum number of iterations
 * @return - the performance measures for the given network
 *
 *
 * Linearizer approximate mean value analysis algorithm
 */

/* Linearizer++ algorithm for closed networks without think times */
fun pfqn_linearizerpp(L: Matrix,
                      N: Matrix,
                      Z: Matrix = Matrix(1, L.numCols),
                      level: Int,
                      tol: Double = 1e-4,
                      maxiter: Int = 1000,
                      flag: Int = 0): pfqnAMVA {
    var Z = Z
    val M = L.numRows
    val R = L.numCols

    Z = if (Z.isEmpty) Matrix(1, R) else Z

    val Nvec = N.toArray1D()
    val Zvec = Z.toArray1D()
    val Lmat = L.toArray2D()
    val J = if (level == 1) {
        1
    } else if (level == 2) {
        1 + R
    } else {
        (R + 1) * (R + 2) / 2
    }

    val X = Array(R) { Array(M) { DoubleArray(J) } }
    val Xn = Array(R) { Array(M) { DoubleArray(J) } }
    val Y = Array(R) { Array(J) { DoubleArray(J) } }
    val W = Array(R) { Array(J) { DoubleArray(J) } }

    var Ntot = 0.0
    for (v in Nvec) {
        Ntot += v
    }

    for (r in 0..<R) {
        for (i in 0..<M) {
            for (v in 0..<J) {
                X[r][i][v] = nvr(Nvec, v + 1, r + 1, R) / M
            }
        }
    }

    for (r in 0..<R) {
        for (u in 0..<J) {
            for (v in 0..<J) {
                Y[r][u][v] = ff(Nvec, r + 1, u + 1, v + 1, R)
            }
        }
    }

    for (r in 0..<R) {
        val Ymat = Array(J) { DoubleArray(J) }
        for (i in 0..<J) {
            for (j in 0..<J) {
                Ymat[i][j] = Y[r][i][j]
            }
        }
        val Ymatrix = Matrix(Ymat)
        val Winv = Matrix.inv(Ymatrix)
        for (i in 0..<J) {
            for (j in 0..<J) {
                W[r][i][j] = Winv[i, j]
            }
        }
    }

    var err = tol + 1
    var iter = 0
    while (err > tol && iter < maxiter) {
        for (r in 0..<R) {
            val N1 = Nvec.clone()
            N1[r] -= 1.0
            for (v in 0..<J) {
                var den = Zvec[r]
                for (j in 0..<M) {
                    var tmp = 1.0
                    for (s in 0..<R) {
                        val vettX = X[s][j]
                        for (u in 0..<J) {
                            val f = ff(N1, s + 1, u + 1, v + 1, R)
                            val vettW = DoubleArray(J)
                            for (k in 0..<J) {
                                vettW[k] = W[s][k][u]
                            }
                            var inner = 0.0
                            for (k in 0..<J) {
                                inner += vettW[k] * vettX[k]
                            }
                            tmp += inner * f
                        }
                    }
                    den += tmp * Lmat[j][r]
                }
                for (i in 0..<M) {
                    var tmp = 1.0
                    for (s in 0..<R) {
                        val vettX = X[s][i]
                        for (u in 0..<J) {
                            val f = ff(N1, s + 1, u + 1, v + 1, R)
                            val vettW = DoubleArray(J)
                            for (k in 0..<J) {
                                vettW[k] = W[s][k][u]
                            }
                            var inner = 0.0
                            for (k in 0..<J) {
                                inner += vettW[k] * vettX[k]
                            }
                            tmp += inner * f
                        }
                    }
                    Xn[r][i][v] = tmp * Lmat[i][r] * nvr(Nvec, v + 1, r + 1, R) / den
                }
            }
        }

        err = 0.0
        for (r in 0..<R) {
            for (i in 0..<M) {
                err = max(err, abs(X[r][i][0] - Xn[r][i][0]))
            }
        }
        err /= Ntot
        iter++
        for (r in 0..<R) {
            for (i in 0..<M) {
                System.arraycopy(Xn[r][i], 0, X[r][i], 0, J)
            }
        }
    }

    val Q = Matrix(M, R)
    for (r in 0..<R) {
        for (i in 0..<M) {
            Q[i, r] = X[r][i][0]
        }
    }

    var Zsum = 0.0
    for (z in Zvec) {
        Zsum += z
    }
    val Xres = Matrix(1, R)
    if (Zsum == 0.0) {
        var Nsum = 0.0
        for (`val` in Nvec) {
            Nsum += `val`
        }
        for (r in 0..<R) {
            var sumQdivL = 0.0
            for (i in 0..<M) {
                sumQdivL += Q[i, r] / Lmat[i][r]
            }
            Xres[0, r] = sumQdivL / (Nsum + M - 1)
        }
    } else {
        throw RuntimeException("pfqn_linearizerpp does not support think times")
    }

    val U = Matrix(M, R)
    val Rmat = Matrix(M, R)
    for (i in 0..<M) {
        for (r in 0..<R) {
            U[i, r] = Xres[0, r] * L[i, r]
            if (Xres[0, r] > 0) {
                Rmat[i, r] = Q[i, r] / Xres[0, r]
            }
        }
    }

    return pfqnAMVA(Q, U, Rmat, null, null, Xres, iter)
}

private fun nv(N: DoubleArray, v: Int, R: Int): DoubleArray {
    var v = v
    val res = N.clone()
    if (v - 1 != 0) {
        v = v - 1
        if (v <= R) {
            res[v - 1] -= 1.0
        } else {
            v = v - R
            outer@ for (i in 1..R) {
                for (j in i..R) {
                    v -= 1
                    if (v == 0) {
                        res[i - 1] -= 1.0
                        res[j - 1] -= 1.0
                        break@outer
                    }
                }
            }
        }
    }
    return res
}

private fun nvr(N: DoubleArray, v: Int, r: Int, R: Int): Double {
    return nv(N, v, R)[r - 1]
}

private fun ff(N: DoubleArray, r: Int, u: Int, v: Int, R: Int): Double {
    var u = u
    val Nv = nv(N, v, R)
    if (u == 1) {
        return Nv[r - 1]
    } else {
        u -= 1
        if (u <= R) {
            return Nv[r - 1] * Nv[u - 1]
        } else {
            u -= R
            for (i in 1..R) {
                for (j in i..R) {
                    u -= 1
                    if (u == 0) {
                        return Nv[r - 1] * Nv[i - 1] * Nv[j - 1]
                    }
                }
            }
        }
    }
    return 0.0
}
/**
 * PFQN linearizerpp algorithms
 */
@Suppress("unused")
class PfqnLinearizerppAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}